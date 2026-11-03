# Build Your Own Hydro1D · 从零构建一维水动力学模型

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python Version](https://img.shields.io/badge/Python-3.10%2B-brightgreen.svg)](requirements.txt)
[![Tests](https://github.com/BeaverMakerLab/build-your-own-hydro1d/actions/workflows/test.yml/badge.svg)](https://github.com/BeaverMakerLab/build-your-own-hydro1d/actions)

本项目是 **【水鉴Lab】** 一维水动力模型系列教程的开源代码库。

旨在帮助水利工程、土木环境、力学与计算机背景的学习者与工程师，**抛开商业软件黑盒，用纯粹、简洁、严谨的 Python 代码从零实现圣维南方程（Saint-Venant Equations / 1D 浅水方程）求解器**。

---

## 教程全章节索引与源码导航

各章节采用**独立自治结构**，各章代码即下即跑、互不耦合。

| 章节 | 核心主题 | 内容类型 / 数值方法 | 源码目录 |
| :---: | :--- | :--- | :---: |
| **01** | 一维水动力计算，走过了怎样一条路 | 水动力学前世今生与发展脉络（理论背景） | - |
| **02** | 拆解圣维南方程：物理意义与特征波 | 质量守恒、动量守恒与重力长波（物理推导） | - |
| **03** | 计算机如何求解方程？连续到离散 | 离散化基本思想与时空网格（概念入门） | - |
| **04** | 有限差分 vs 有限体积：本质区别 | 通量守恒性、激波捕捉对比（数值对比） | - |
| **05** | 水动力模型为什么会算崩？CFL 条件 | 双曲型方程特征线与稳定性（稳定理论） | - |
| **06** | **从零构建一维求解器：让圣维南方程真正跑起来** | **一阶 FVM + Rusanov 通量 + 显式推进** | [**`ch06_fvm_dambreak`**](chapters/ch06_fvm_dambreak/) |

---

## 快速上手（Quick Start）

### 1. 克隆代码仓库

```bash
git clone https://github.com/BeaverMakerLab/build-your-own-hydro1d.git
cd build-your-own-hydro1d
```

### 2. 安装基础依赖

本项目追求极简依赖，仅需基础科学计算库：

```bash
pip install -r requirements.txt
```

### 3. 运行第 06 章湿床溃坝基准算例

```bash
cd chapters/ch06_fvm_dambreak
python run_experiment.py
```

### 4. 运行全套单元测试

```bash
python -m unittest discover -s chapters/ch06_fvm_dambreak -p "test_*.py" -v
```

---

## 第 06 章算例输出说明

第 06 章执行后，将在 `chapters/ch06_fvm_dambreak/data/` 目录中生成全时空高精度水位数据：

1. **`water_level_matrix.csv`（断面-时刻矩阵表，推荐）**：
   - 每行为一个时间步（`0.0s, 1.0s, ..., 60.0s`）；
   - 每列为一个空间断面中心坐标（`x_2.50m, x_7.50m, ..., x_997.50m`）；
   - 单元格为对应断面的实时水位高度（平底渠道 $z_b=0\text{ m}$，水位严格等于水深 $h$）。
   - Excel 双击即可直接打开，方便快速透视。
2. **`water_level_series.csv`（时空长明细表）**：
   - 包含 `[time_s, x_m, water_level_m, unit_discharge_q_m2s, velocity_u_ms]` 五列长数据。
3. **`simulation_summary.json`（算例参数与守恒误差）**：
   - 记录时间步数、CFL、Stoker 理论激波平台理论值、质量守恒误差等指标。

---

## 许可证 (License)

本项目遵循 [MIT License](LICENSE) 开源协议。
欢迎学习、引用与交流。商业使用及衍生工程请保留原始著作权声明。
