"""Run 1D Saint-Venant wet-bed dam-break simulation and export water level data.

No plotting or GUI dependencies required. Outputs cross-sectional water level
data at each time snapshot to CSV files.
"""
import json
from pathlib import Path
import numpy as np
from solver import (
    make_dambreak,
    simulate,
    export_water_level_matrix,
    export_water_level_series,
)
from exact_solution import intermediate_state

CHAPTER_ROOT = Path(__file__).resolve().parent
DATA_DIR = CHAPTER_ROOT / 'data'
DATA_DIR.mkdir(parents=True, exist_ok=True)


def run_benchmark(
    length: float = 1000.0,
    nx: int = 200,
    left_depth: float = 2.0,
    right_depth: float = 1.0,
    split: float = 500.0,
    t_final: float = 60.0,
    dt_snapshot: float = 1.0,
    cfl: float = 0.4,
):
    """Run the dam-break case and export water level time-series."""
    print("=" * 60)
    print("一维圣维南方程湿床溃坝算例求解器 (无绘图依赖版本)")
    print("=" * 60)
    print(f"网格与物理参数:")
    print(f"  - 河道长度: {length:.1f} m, 单元格数: {nx} (断面间距 dx = {length/nx:.2f} m)")
    print(f"  - 初始水深: 左侧 h_L = {left_depth:.2f} m, 右侧 h_R = {right_depth:.2f} m (位置 x = {split:.1f} m)")
    print(f"  - 模拟总时长: {t_final:.1f} s, 记录时间间隔 dt = {dt_snapshot:.2f} s, CFL = {cfl}")
    print(f"  - 水平平底河道 (z_b = 0 m)，因此各断面水位严格等于水深 h")

    x, dx, U0 = make_dambreak(
        length=length,
        nx=nx,
        left_depth=left_depth,
        right_depth=right_depth,
        split=split,
    )

    sim = simulate(
        U0,
        dx=dx,
        t_final=t_final,
        snapshot_interval=dt_snapshot,
        cfl=cfl,
    )

    n_snapshots = len(sim.snapshots)
    print(f"\n求解计算完成:")
    print(f"  - 总计算步数: {sim.steps} 步")
    print(f"  - 实际记录时刻数: {n_snapshots} 个时刻 (t = 0.0 s 到 {sim.time:.1f} s)")
    print(f"  - 最小水深: {sim.min_depth:.4f} m (无负水深)")
    print(f"  - 质量守恒误差: {sim.mass_balance_error:.2e} m^2 (相对误差: {abs(sim.mass_balance_error)/sim.mass_initial:.2e})")

    matrix_file = DATA_DIR / 'water_level_matrix.csv'
    export_water_level_matrix(sim, matrix_file)
    print(f"\n数据导出成功:")
    print(f"  [1] 断面-时刻水位矩阵表: {matrix_file.name}")
    print(f"      格式: 每行为一个记录时刻，每列对应一个断面坐标 x_im，值为该时刻该断面的水位(m)")

    series_file = DATA_DIR / 'water_level_series.csv'
    export_water_level_series(sim, series_file)
    print(f"  [2] 逐断面逐时刻明细长表: {series_file.name}")
    print(f"      格式: [time_s, x_m, water_level_m, unit_discharge_q_m2s, velocity_u_ms]")

    summary = {
        'scenario': '1D wet-bed dam break (Stoker benchmark)',
        'channel_length_m': length,
        'nx': nx,
        'dx_m': dx,
        't_final_s': t_final,
        'dt_snapshot_s': dt_snapshot,
        'n_snapshots': n_snapshots,
        'n_steps': sim.steps,
        'min_depth_m': sim.min_depth,
        'mass_initial_m2': sim.mass_initial,
        'mass_final_m2': sim.mass_final,
        'mass_balance_error_m2': sim.mass_balance_error,
        'stoker_intermediate_plateau': intermediate_state(left_depth, right_depth),
        'files': {
            'matrix_csv': str(matrix_file.name),
            'series_csv': str(series_file.name),
        },
    }
    summary_file = DATA_DIR / 'simulation_summary.json'
    summary_file.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f"  [3] 算例元数据配置摘要: {summary_file.name}")
    print("=" * 60)
    return sim


if __name__ == '__main__':
    run_benchmark()
