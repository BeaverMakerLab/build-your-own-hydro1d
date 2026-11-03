"""Minimal 1-D nonlinear shallow-water solver, for hydro education only.

Governing equations on a horizontal, rectangular channel, per unit width:
 U_t + F(U)_x = 0
 U = [h, q=hu], F(U) = [q, q^2/h + g*h^2/2].
Finite-volume first order + Rusanov flux + forward Euler.
No bed slope, friction, lateral flow, dry states, structures or networks.
"""
from __future__ import annotations
import csv
from pathlib import Path
from dataclasses import dataclass
from typing import Iterable, Sequence
import numpy as np

GRAVITY = 9.81


def physical_flux(U: np.ndarray, g: float = GRAVITY) -> np.ndarray:
    h = U[..., 0]
    q = U[..., 1]
    if np.any(h <= 0.0):
        raise ValueError("This minimal wet-bed solver requires h > 0 everywhere.")
    return np.stack((q, q * q / h + 0.5 * g * h * h), axis=-1)


def rusanov_flux(left: np.ndarray, right: np.ndarray, g: float = GRAVITY) -> np.ndarray:
    """Conservative Rusanov (local Lax-Friedrichs) interface flux."""
    hL, hR = left[..., 0], right[..., 0]
    qL, qR = left[..., 1], right[..., 1]
    if np.any(hL <= 0) or np.any(hR <= 0):
        raise ValueError("Dry states are outside this introductory model.")
    smax = np.maximum(np.abs(qL / hL) + np.sqrt(g * hL),
                      np.abs(qR / hR) + np.sqrt(g * hR))
    return 0.5 * (physical_flux(left, g) + physical_flux(right, g)) - 0.5 * smax[..., None] * (right - left)


def make_dambreak(length: float = 1000.0, nx: int = 200,
                  left_depth: float = 2.0, right_depth: float = 1.0,
                  split: float = 500.0):
    """Cell-center coordinates and initial wet-dam-break state."""
    dx = length / nx
    x = (np.arange(nx, dtype=float) + 0.5) * dx
    U = np.column_stack((np.where(x < split, left_depth, right_depth),
                         np.zeros(nx, dtype=float)))
    return x, dx, U


@dataclass
class Simulation:
    x: np.ndarray
    dx: float
    state: np.ndarray
    time: float
    snapshots: dict[float, np.ndarray]
    steps: int
    mass_initial: float
    mass_final: float
    mass_boundary_accounting: float
    mass_balance_error: float
    min_depth: float


def simulate(U0: np.ndarray, dx: float, t_final: float,
             snapshot_times: Sequence[float] = (),
             snapshot_interval: float | None = None,
             cfl: float = 0.4, g: float = GRAVITY) -> Simulation:
    """Run with transmissive ghost cells. Waves must not reach boundaries.

    The mass conservation budget explicitly includes the numerical boundary flux.
    """
    if not (0 < cfl <= 0.5):
        raise ValueError("For this teaching solver use 0 < cfl <= 0.5.")
    U = np.asarray(U0, dtype=float).copy()
    if U.ndim != 2 or U.shape[1] != 2 or np.any(U[:, 0] <= 0):
        raise ValueError("State must be [nx, 2] and strictly wet.")
    nx = len(U)
    x = (np.arange(nx) + 0.5) * dx
    mass0 = dx * np.sum(U[:, 0])
    total_boundary_mass = 0.0

    target_times = set([0.0, float(t_final), *[float(t) for t in snapshot_times]])
    if snapshot_interval is not None and snapshot_interval > 0:
        n_interval_steps = int(np.floor(t_final / snapshot_interval))
        interval_times = [float(i * snapshot_interval) for i in range(n_interval_steps + 1)]
        target_times.update(interval_times)

    times = sorted(target_times)
    if times[0] < 0 or times[-1] > t_final + 1e-12:
        raise ValueError("Invalid snapshot times")

    snapshots = {0.0: U.copy()}
    t, steps = 0.0, 0
    for target in times[1:]:
        while t < target - 1e-11:
            h, q = U[:, 0], U[:, 1]
            smax = np.max(np.abs(q / h) + np.sqrt(g * h))
            dt = min(cfl * dx / smax, target - t)
            padded = np.concatenate((U[0:1], U, U[-1:]), axis=0)
            fluxes = rusanov_flux(padded[:-1], padded[1:], g)
            U -= (dt / dx) * (fluxes[1:] - fluxes[:-1])
            if not np.all(np.isfinite(U)) or np.min(U[:, 0]) <= 0:
                raise FloatingPointError("Negative or invalid water depth; stop simulation.")
            total_boundary_mass += dt * (fluxes[0, 0] - fluxes[-1, 0])
            t += dt
            steps += 1
        snapshots[target] = U.copy()

    mass1 = dx * np.sum(U[:, 0])
    return Simulation(x=x, dx=dx, state=U, time=t, snapshots=snapshots,
                      steps=steps, mass_initial=mass0, mass_final=mass1,
                      mass_boundary_accounting=total_boundary_mass,
                      mass_balance_error=(mass1 - mass0) - total_boundary_mass,
                      min_depth=float(np.min(U[:, 0])))


def export_water_level_matrix(sim: Simulation, filepath: str | Path) -> Path:
    """Export water level as a 2D matrix: rows are time steps, columns are cross-sections (x).

    Note: In this horizontal flat-bed benchmark (z_b = 0), water level equals water depth h.
    """
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ['time_s'] + [f'x_{xi:.2f}m' for xi in sim.x]
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for t in sorted(sim.snapshots.keys()):
            h_row = sim.snapshots[t][:, 0]
            writer.writerow([f'{t:.3f}'] + [f'{val:.6f}' for val in h_row])
    return path


def export_water_level_series(sim: Simulation, filepath: str | Path) -> Path:
    """Export detailed cross-section water levels and hydraulic variables in long table format."""
    path = Path(filepath)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ['time_s', 'x_m', 'water_level_m', 'unit_discharge_q_m2s', 'velocity_u_ms']
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        for t in sorted(sim.snapshots.keys()):
            snap = sim.snapshots[t]
            h = snap[:, 0]
            q = snap[:, 1]
            u = q / h
            for xi, hi, qi, ui in zip(sim.x, h, q, u):
                writer.writerow([f'{t:.3f}', f'{xi:.2f}', f'{hi:.6f}', f'{qi:.6f}', f'{ui:.6f}'])
    return path


if __name__ == '__main__':
    print("=== 一维圣维南方程湿床溃坝求解器 ===")
    x, dx, U0 = make_dambreak(length=1000.0, nx=200, left_depth=2.0, right_depth=1.0, split=500.0)
    sim = simulate(U0, dx, t_final=60.0, snapshot_interval=1.0)
    print(f"计算完成: 断面数 = {len(x)}, 计算总步数 = {sim.steps}, 质量守恒误差 = {sim.mass_balance_error:.2e} m^2")
    matrix_path = export_water_level_matrix(sim, "water_level_matrix.csv")
    print(f"水位矩阵数据已导出至: {matrix_path} (时刻数: {len(sim.snapshots)}, 断面数: {len(x)})")
