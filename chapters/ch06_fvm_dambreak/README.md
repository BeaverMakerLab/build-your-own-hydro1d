# 第 06 章：从零构建一维有限体积水动力求解器（湿床溃坝）

本章为【水鉴Lab 一维水动力模型教程】第 06 篇配套代码，实现了一个紧凑、自包含、数学严谨的一维圣维南方程（浅水方程）有限体积求解器。

## 物理模型与数值方案

- **控制方程**：水平矩形明渠、等宽、无底坡、无糙率摩阻的非线性浅水方程（Saint-Venant Equations）：
  $$\frac{\partial U}{\partial t} + \frac{\partial F(U)}{\partial x} = 0$$
  状态量 $U = [h, q]^T$，通量 $F(U) = [q, q^2/h + \frac{1}{2}gh^2]^T$（单位宽流量 $q = hu$）。
- **空间离散**：一阶有限体积法（Cell-Centered FVM）+ Rusanov（Local Lax-Friedrichs）数值界面通量。
- **时间推进**：显式一阶 Euler 格式，动态满足 CFL 条件（$\text{CFL} \le 0.4$）。
- **边界条件**：零梯度（Transmissive / Extrapolation）边界条件。
- **验证基准**：一维湿床溃坝（Stoker Riemann 问题），附带完全独立的解析参考解。

## 文件结构

- `solver.py`：核心求解器（包含控制方程通量计算、Rusanov 格式、显式推进及水位导出函数）。
- `exact_solution.py`：Stoker 湿床溃坝精确 Riemann 解析解（通过 `scipy.optimize.brentq` 独立求根）。
- `run_experiment.py`：算例执行入口，计算 $0 \sim 60\text{ s}$ 溃坝演进，自动导出全断面逐时刻水位。
- `test_solver.py`：5 项单元测试（静水、均匀动水、通量相容、质量守恒闭环及空间收敛性）。
- `data/`：
  - `water_level_matrix.csv`：**断面-时刻水位矩阵表**（行=时间，列=断面坐标，单元格=水位值）。
  - `water_level_series.csv`：**时空长明细表**（时间、断面、水位、单宽流量、流速）。
  - `simulation_summary.json`：算例参数与守恒误差总结。

## 快速运行

进入本章目录：

```bash
cd chapters/ch06_fvm_dambreak

# 1. 运行算例并导出全断面水位
python run_experiment.py

# 2. 运行单元测试
python -m unittest discover -s . -p "test_*.py" -v
```
