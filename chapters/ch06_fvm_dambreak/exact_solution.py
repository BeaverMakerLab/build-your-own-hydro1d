"""Independent exact wet-bed dam-break (Stoker) Riemann reference solution.

Assumptions: left depth > right depth > 0, zero velocity on both sides,
flat bed, no friction, no boundary reflections. Reference is not computed
by calling or reusing the finite-volume numerical solver.
"""
import numpy as np
from scipy.optimize import brentq


def intermediate_state(h_left=2., h_right=1., g=9.81):
    if not (h_left > h_right > 0):
        raise ValueError('Requires h_left > h_right > 0')
    cl = np.sqrt(g*h_left)
    def root(h):
        return 2*(np.sqrt(g*h)-cl) + (h-h_right)*np.sqrt(0.5*g*(1/h + 1/h_right))
    hm = brentq(root, h_right, h_left, xtol=1e-13)
    um = 2*(cl - np.sqrt(g*hm))
    shock_speed = hm*um/(hm-h_right)
    return {'h_middle':float(hm), 'u_middle':float(um),
            'left_fan_head':float(-cl),
            'left_fan_tail':float(um-np.sqrt(g*hm)),
            'right_shock_speed':float(shock_speed)}


def exact_dambreak(x, time, split=500., h_left=2., h_right=1., g=9.81):
    x = np.asarray(x, dtype=float)
    if time < 0:
        raise ValueError('time must be nonnegative')
    if time == 0:
        return np.where(x < split, h_left, h_right), np.zeros_like(x)
    S = intermediate_state(h_left,h_right,g)
    xi = (x-split)/time
    cl = np.sqrt(g*h_left)
    head, tail, shock = S['left_fan_head'],S['left_fan_tail'],S['right_shock_speed']
    h = np.full_like(x, h_right)
    u = np.zeros_like(x)
    left = xi < head
    fan = (xi >= head) & (xi < tail)
    middle = (xi >= tail) & (xi < shock)
    h[left] = h_left
    c_fan = (2*cl - xi[fan])/3.
    h[fan] = c_fan*c_fan/g
    u[fan] = (2./3.)*(cl + xi[fan])
    h[middle] = S['h_middle']
    u[middle] = S['u_middle']
    return h, h*u
