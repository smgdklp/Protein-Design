# 多肽链关节级运动学与软约束模型

---

## 记号约定

| 记号 | 含义 |
|---|---|
| $n$ | 残基数，链方向 N→C，$i=1,\dots,n$ |
| $B_i$ | 第 $i$ 个残基对应的二自由度复合节点（$\phi_i,\psi_i$ 两条 1R 副串联） |
| $J_i$ | 第 $i$ 个肽键对应的 1R 转动副（$\omega_i$），$i=1,\dots,n-1$ |
| $q\in\mathbb T^{3n-1}$ | 关节位形向量 |
| $S_k\in\mathbb R^6$ | 零位构型下的空间旋量 |
| $M=T(0)$ | 零位构型下末端位姿 |
| $f(q)$ | 正运动学映射 |
| $p(\phi,\psi\mid m_{i-1},m_i,m_{i+1})$ | 由经验数据拟合的条件概率密度（软约束来源） |
| $\omega^0(m_i,m_{i+1})\in\{0,\pi\}$ | 肽键角的离散映射值（cis / trans） |
| $\beta$ | 统计约束强度（软约束的"刚度"，$\beta=1$ 为纯统计） |

---

# 一、模型部分

## 1. 模型定义与适用范围

将含 $n$ 个残基的多肽链抽象为一条沿 **N→C** 方向串联的开链机构。模型仅在**关节（joint）层级**描述运动传递，刚体内部几何视为常量。

- 第 $i$ 个残基抽象为二自由度复合节点 $B_i$，串联两个相邻的单自由度转动副，即 $\phi_i$ 与 $\psi_i$；
- 第 $i$ 个肽键抽象为单自由度转动副 $J_i$，自由度为 $\omega_i$；
- 各刚体间的连接仅通过上述转动副实现，不含移动副、球副等其他关节型式。

## 2. 符号与索引

链方向为 N→C，残基索引 $i=1,\dots,n$。

全部关节按链序统一编号为 $k=1,\dots,3n-1$，映射关系为

$$
\phi_i \leftrightarrow k=3i-2,\qquad
\psi_i \leftrightarrow k=3i-1,\qquad
\omega_i \leftrightarrow k=3i\ \ (i=1,\dots,n-1).
$$

关节序列显式写为

$$
\phi_1,\ \psi_1,\ \omega_1,\ \phi_2,\ \psi_2,\ \omega_2,\ \dots,\ \omega_{n-1},\ \phi_n,\ \psi_n .
$$

位形变量取

$$
q=(q_1,\dots,q_{3n-1})\in Q,\qquad
q_{3i-2}=\phi_i,\quad q_{3i-1}=\psi_i,\quad q_{3i}=\omega_i .
$$

无约束条件下 $Q=\mathbb T^{\,3n-1}=(S^1)^{3n-1}$，各角量以 $2\pi$ 为周期。

## 3. 关节旋量

对第 $k$ 个转动副，定义**零位构型**（$q=0$）下的空间旋量

$$
S_k=\begin{bmatrix}\boldsymbol a_k\\[2pt] \boldsymbol v_k\end{bmatrix}\in\mathbb R^6,
\qquad
\boldsymbol v_k=-\boldsymbol a_k\times \boldsymbol r_k,
\qquad
\|\boldsymbol a_k\|=1,
$$

其中 $\boldsymbol a_k\in\mathbb R^3$ 为关节轴单位向量，$\boldsymbol r_k\in\mathbb R^3$ 为轴上任一点。该旋量螺距为零，对应纯转动

$$
T_k(q_k)=\exp\!\big([S_k]\,q_k\big)\in SE(3),
$$

$[S_k]$ 为旋量的 $4\times4$ 矩阵表示。三组关节轴均在零位构型下定义，正转方向按右手定则取轴正向：

| 关节 | 轴单位向量 $\boldsymbol a_k$ | 轴上点 $\boldsymbol r_k$ | 转动时保持不动的刚体 | 随动段 |
|---|---|---|---|---|
| $\phi_i$ | $\dfrac{C_{\alpha i}-N_i}{\lVert C_{\alpha i}-N_i\rVert}$ | $N_i$ | $N_i,\ C_{\alpha i}$ | $C_i$ 及下游 |
| $\psi_i$ | $\dfrac{C_i-C_{\alpha i}}{\lVert C_i-C_{\alpha i}\rVert}$ | $C_{\alpha i}$ | $C_{\alpha i},\ C_i$ | $N_{i+1}$ 及下游 |
| $\omega_i$ | $\dfrac{N_{i+1}-C_i}{\lVert N_{i+1}-C_i\rVert}$ | $C_i$ | $C_i,\ N_{i+1}$ | $C_{\alpha i+1}$ 及下游 |

**注记.** $B_i$ 标称二自由度，实为两条在 $C_{\alpha i}$ 处相交的独立转动副串联而成，即两个 1R 关节的复合，**不能合并为一条单一二自由度轴**。因子分解

$$
T_{B_i}(\phi_i,\psi_i)=\exp\!\big([S_{i,\phi}]\,\phi_i\big)\,\exp\!\big([S_{i,\psi}]\,\psi_i\big)
$$

只是链式乘积中的**记号分组**；两个旋量必须统一表示在零位空间坐标系下，乘积才成立。

## 4. 链式结构与自由度

结构序列为

$$
B_1\to J_1\to B_2\to\cdots\to J_{n-1}\to B_n ,
$$

等价关节序列见 §2。自由度计数为

$$
\underbrace{2n}_{\phi_i,\ \psi_i}+\underbrace{(n-1)}_{\omega_i}=3n-1 .
$$

因此该机构共有 $3n-1$ 个独立关节变量，位形空间维数为 $3n-1$。

## 5. 基于 PoE 的正运动学

设 $\{S\}$ 为固定于基座的空间坐标系。取

$$
M\equiv T(0)\in SE(3)
$$

为**零位构型**下末端坐标系相对 $\{S\}$ 的齐次变换，由零位链路几何完全确定。

空间形式正运动学为

$$
T(q)=\left[\prod_{k=1}^{3n-1}\exp\!\big([S_k]\,q_k\big)\right]M
=\exp\!\big([S_1]q_1\big)\cdots\exp\!\big([S_{3n-1}]q_{3n-1}\big)\,M,
$$

其中各 $S_k$ 均取零位构型下的空间旋量，乘积严格按 N→C 链序排列。

等价地可采用物体形式：

$$
T(q)=M\prod_{k=1}^{3n-1}\exp\!\big([B_k]\,q_k\big),
\qquad
B_k=\mathrm{Ad}_{M^{-1}}\,S_k,
$$

其中 $\mathrm{Ad}$ 为伴随变换。正运动学映射记为

$$
f:Q\to SE(3),\qquad q\mapsto f(q)=T(q).
$$

**零位参考说明.** $M$ 应取整条链的零位末端位姿 $T(0)$；若将末端坐标系原点约定在 $B_1$ 的零位，则 $M$ 即由该约定与其余刚体常量共同确定。其选择属坐标系规范问题，不影响自由度与运动学等价性。

## 6. 关节与输出函数

单元变换（作为链式乘积中的因子）：

$$
T_{B_i}(\phi_i,\psi_i)=\exp\!\big([S_{i,\phi}]\phi_i\big)\exp\!\big([S_{i,\psi}]\psi_i\big),
\qquad
T_{J_i}(\omega_i)=\exp\!\big([S_{i,\omega}]\omega_i\big).
$$

整链输出函数为 $T(q)=f(q)\in SE(3)$。

空间雅可比由标准 PoE 给出：

$$
J_s(q)=\begin{bmatrix}S_1'&S_2'&\cdots&S_{3n-1}'\end{bmatrix}\in\mathbb R^{6\times(3n-1)},
$$

$$
S_k'=\mathrm{Ad}_{\,e^{[S_1]q_1}\cdots e^{[S_{k-1}]q_{k-1}}}\,S_k,
\qquad S_1'=S_1 .
$$

末端速度关系为 $\mathcal V_s=J_s(q)\,\dot q$。

## 7. 结构常量（零位构型）

以下常量用于确定 $M$ 与各旋量，属刚体几何参数，不随关节运动改变：

| 参数 | 典型值 | 说明 |
|---|---:|---|
| $N-C_\alpha$ 键长 | $\approx 1.46\ \text{\AA}$ | 主链键长 |
| $C_\alpha-C$ 键长 | $\approx 1.52\ \text{\AA}$ | 主链键长 |
| $C-N$ 肽键长 | $\approx 1.33\ \text{\AA}$ | 部分双键 |
| $\angle N-C_\alpha-C$ | $\approx 111^\circ$ | 主链键角 |
| $\angle C_\alpha-C-N$ | $\approx 116^\circ$ | 主链键角 |
| $\angle C-N-C_\alpha$ | $\approx 121^\circ$ | 肽键附近键角 |
| $L_i$ | $L_1,\dots,L_{20}$ | 残基标签，映射零位坐标与几何参数，不计入连续自由度 |

标签映射记为

$$
L_i\ \mapsto\ \bigl(\text{原子类型},\ T_{B_i}(0,0),\ \text{零位几何参数}\bigr).
$$

## 8. 小结

关节级正运动学可概括为

$$
\boxed{
T(q)=\left[\prod_{k=1}^{3n-1}\exp\!\big([S_k]\,q_k\big)\right]M,
\qquad
M=T(0),
\qquad
q\in\mathbb T^{\,3n-1}
}
$$

其中 $q_k$ 按 N→C 链序取 $\phi_i,\psi_i,\omega_i$；$S_k$ 为零位空间旋量；标签 $L_i$ 为离散几何填充变量，不增加连续自由度。

---

# 二、约束部分

本节给出模型所受的两条约束的**声明**。

## 9. 约束声明

模型受且仅受以下两条约束。

### 9.1 约束一：$\phi/\psi$ 的条件统计软约束

对每个残基 $i=1,\dots,n$，$(\phi_i,\psi_i)$ 受给定相邻标签下的条件概率密度约束：

$$
\text{(C1)}\qquad
\rho_{\phi\psi}\ \propto\ \prod_{i=1}^{n}
p\big(\phi_i,\psi_i\mid m_{i-1},m_i,m_{i+1}\big)^{\beta},
\qquad \beta>0 .
$$

- 性质：**软约束**，不设不可达区；
- 禁区为高势能（低概率）区域，其密度由正地板限定，不取零；
- 端点 $i=1,n$ 的三元组条件不成立，按降阶条件处理，或令其势能为常数。

### 9.2 约束二：$\omega$ 的离散映射

对每个肽键 $i=1,\dots,n-1$，其转角不保留连续自由度，由相邻残基标签唯一确定：

$$
\text{(C2)}\qquad
\omega_i=\omega^{0}(m_i,m_{i+1})\in\{0,\ \pi\},
\qquad i=1,\dots,n-1 .
$$

- 性质：**离散赋值**；
- **不含活动范围**，不设 $\pm$ 容差带；
- $0$ 对应 cis，$\pi$ 对应 trans；映射表由标签对给定，不随位形变化。

### 9.3 约束范围

除上述两条外不引入其他约束：键长与键角取 §7 的零位常量；碰撞条件不在本模型内讨论；关节型式仅限转动副。

| 编号 | 类型 | 作用变量 | 形式 |
|---|---|---|---|
| C1 | 条件统计软约束 | $q_{3i-2}=\phi_i,\ q_{3i-1}=\psi_i$ | 权重 $\prod_{i}p_i^{\beta}$ |
| C2 | 离散映射 | $q_{3i}=\omega_i$ | $\omega^{0}(m_i,m_{i+1})\in\{0,\pi\}$ |

---

# 三、最终简化模型

将 §9 的两条约束代入 §一 的机构，得到降自由度后的简化模型。

## 10. 简化位形与自由度

记连续坐标为

$$
\boxed{
\tilde q=(\phi_1,\psi_1,\ \phi_2,\psi_2,\ \dots,\ \phi_n,\psi_n)\in\mathbb T^{\,2n}
}
$$

肽键角取离散映射值

$$
\omega_i=\omega^{0}(m_i,m_{i+1})\in\{0,\pi\},\qquad i=1,\dots,n-1 .
$$

| 量 | 简化前 | 简化后 |
|---|---|---|
| 连续坐标 | $q\in\mathbb T^{3n-1}$ | $\tilde q\in\mathbb T^{2n}$ |
| 离散赋值 | 无 | $\omega_i=\omega^{0}(m_i,m_{i+1})$，共 $n-1$ 个 |
| 连续自由度 | $3n-1$ | $2n$ |

## 11. 简化后的正运动学

按 N→C 链序代入 $\omega_i=\omega^{0}(m_i,m_{i+1})$：

$$
\boxed{
T(\tilde q)=\left[\prod_{i=1}^{n-1}T_{B_i}(\phi_i,\psi_i)\,T_{J_i}\big(\omega^{0}(m_i,m_{i+1})\big)\right]
T_{B_n}(\phi_n,\psi_n)\,M
}
$$

其中

$$
T_{B_i}(\phi_i,\psi_i)=\exp\!\big([S_{i,\phi}]\phi_i\big)\exp\!\big([S_{i,\psi}]\psi_i\big),
\qquad
T_{J_i}(\omega_i)=\exp\!\big([S_{i,\omega}]\omega_i\big).
$$

等价地，保留原关节编号，令

$$
q_{3i-2}=\phi_i,\qquad q_{3i-1}=\psi_i,\qquad q_{3i}=\omega^{0}(m_i,m_{i+1}),
$$

则

$$
T(\tilde q)=\left[\prod_{k=1}^{3n-1}\exp\!\big([S_k]\,q_k\big)\right]M .
$$

## 12. 简化后的位形权重

$$
\boxed{
\rho(\tilde q)\ \propto\ \exp\!\big[-\beta\,U(\tilde q)\big],
\qquad
U(\tilde q)=\sum_{i=1}^{n}u_i(\phi_i,\psi_i),
\qquad
u_i=-\log p\big(\phi_i,\psi_i\mid m_{i-1},m_i,m_{i+1}\big)
}
$$

## 13. 简化后的雅可比

$\omega_i$ 固定后不再参与变分，$\mathrm d q_{3i}=0$，故速度关系只保留连续列：

$$
\boxed{
\mathcal V_s=J_{\rm red}(\tilde q)\,\dot{\tilde q},
\qquad
J_{\rm red}=\begin{bmatrix}S'_{1,\phi}&S'_{1,\psi}&S'_{2,\phi}&S'_{2,\psi}&\cdots&S'_{n,\phi}&S'_{n,\psi}\end{bmatrix}
\in\mathbb R^{6\times 2n}
}
$$

**注记.** $\omega$ 固定只去掉了对应的变分自由度，**并未从乘积中删去其指数因子**；$J_{\rm red}$ 中的旋量仍按包含 $\exp([S_{i,\omega}]\omega_i^{0})$ 的完整链计算。

## 14. 最终模型小结

$$
\boxed{
\begin{aligned}
&\text{位形：}\quad \tilde q=(\phi_1,\psi_1,\dots,\phi_n,\psi_n)\in\mathbb T^{2n},
\qquad \omega_i=\omega^{0}(m_i,m_{i+1})\in\{0,\pi\};\\[4pt]
&\text{正运动学：}\quad T(\tilde q)=\Big[\prod_{i=1}^{n-1}T_{B_i}(\phi_i,\psi_i)\,T_{J_i}(\omega_i^{0})\Big]T_{B_n}(\phi_n,\psi_n)\,M;\\[4pt]
&\text{权重：}\quad \rho(\tilde q)\propto\exp\!\Big[-\beta\sum_{i=1}^{n}u_i(\phi_i,\psi_i)\Big];\\[4pt]
&\text{速度：}\quad \mathcal V_s=J_{\rm red}(\tilde q)\,\dot{\tilde q},\qquad J_{\rm red}\in\mathbb R^{6\times 2n}.
\end{aligned}
}
$$

简化模型相对原机构的变化只有三处：连续自由度由 $3n-1$ 降为 $2n$；肽键角由关节变量变为标签函数 $\omega^{0}\in\{0,\pi\}$；位形权重由条件统计软约束给出。$S_k$、$M$ 与各单元变换的定义均不变。
