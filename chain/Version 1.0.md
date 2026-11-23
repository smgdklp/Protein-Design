# 肽键链模型

---

## 1. 简化为主碳链的模型

### 1.1 最简表述

模型对象为残基 \(0,\dots,n+1\) 的主链。残基 \(1..n\) 为待定窗口，残基 \(0\) 与 \(n+1\) 为上下文残基，其绝对位姿 \(T_0^*,T_{n+1}^*\) 已知。

主链约化为 \(n+2\) 个主碳节点，每个节点携带一个局部标架；相邻节点之间由一段几何固定的双铰链连杆连接。于是整条链的构型由

$$
\big(\ T_0^*,\ \Theta\ \big),\qquad
\Theta=(\psi_0,\varphi_1,\psi_1,\dots,\varphi_n,\psi_n,\varphi_{n+1})\in\mathbb T^{2n+2}
$$

唯一确定。其中序列、常数与肽键角均视为已知量；自由度数目为

$$
1+2n+1=2n+2 .
$$

### 1.2 统一代数

将残基 \(0..n+1\) 的主链原子按链序记为

$$
Y_{3i}=N_i,\qquad Y_{3i+1}=C\alpha_i,\qquad Y_{3i+2}=C_i,\qquad i=0,\dots,n+1,
$$

即

$$
(Y_0,\dots,Y_{3n+5})=(N_0,C\alpha_0,C_0,\;N_1,C\alpha_1,C_1,\;\dots,\;N_{n+1},C\alpha_{n+1},C_{n+1}),
$$

共 \(3n+6\) 个原子、\(3n+5\) 条键、\(3n+4\) 个内角。下文统一采用该链序。

$$
\ell_k=|Y_k-Y_{k-1}|\quad(k=1..3n+5),\qquad
\beta_k=\angle(Y_{k-1},Y_k,Y_{k+1})\quad(k=1..3n+4),
$$

$$
\tau_k=\operatorname{Dih}(Y_{k-2},Y_{k-1},Y_k,Y_{k+1})\quad(k=2..3n+4).
$$

三个扭转与链上位置的对应为

$$
\tau_{3i+2}=\psi_i,\qquad
\tau_{3i+3}=\omega_i,\qquad
\tau_{3i+4}=\varphi_{i+1},\qquad i=0..n.
$$

### 1.3 常量：键长与键角

键长：

$$
\ell^{C\alpha C}_i=|C\alpha_i-C_i|,\qquad
\ell^{CN}_i=|C_i-N_{i+1}|,\qquad
\ell^{N C\alpha}_i=|N_{i+1}-C\alpha_{i+1}|,\qquad i=0..n,
$$

$$
\ell^{N C\alpha}_{-1}=|N_0-C\alpha_0|,\qquad
\ell^{C\alpha C}_{n+1}=|C\alpha_{n+1}-C_{n+1}| .
$$

合计 \(3(n+1)+2=3n+5\) 条。

键角：

$$
\kappa_i=\angle(N_i,C\alpha_i,C_i)\quad(i=0..n+1),\qquad
\eta_i=\angle(C\alpha_i,C_i,N_{i+1}),\qquad
\theta_i=\angle(C_i,N_{i+1},C\alpha_{i+1})\quad(i=0..n).
$$

合计 \((n+2)+(n+1)+(n+1)=3n+4\) 个。

与旧记号之对照为：\(a_i=\kappa_i\)，\(b_i=\theta_i\)，\(c_i=\kappa_{i+1}\)。其中 \(c_i\) 与 \(a_{i+1}\) 重复，旧记号未含 \(\eta_i\)。

### 1.4 方向与定义

**（一）键角的定义与方向。** 对任意三点 \(A,B,C\)，定义

$$
\angle(A,B,C)=\arccos\frac{(A-B)\cdot(C-B)}{|A-B|\,|C-B|}\in(0,\pi),
$$

其顶点为 \(B\)，两条边分别指向 \(A\) 与 \(C\)，沿链方向 \(Y_{k-1}\to Y_k\to Y_{k+1}\) 取值。递推中采用其补角

$$
\alpha_k=\pi-\beta_k,
$$

即入键 \(Y_{k-1}\to Y_k\) 与出键 \(Y_k\to Y_{k+1}\) 之间的夹角。

**（二）二面角的定义。** 对四点 \(A,B,C,D\)，定义

$$
\operatorname{Dih}(A,B,C,D)=\operatorname{atan2}\!\big((u\times w)\cdot\hat e,\ u\cdot w\big),
$$

$$
\hat e=\frac{C-B}{|C-B|},\qquad u=(A-B)_\perp,\qquad w=(D-C)_\perp,\qquad x_\perp=x-(x\cdot\hat e)\hat e .
$$

该式取 IUPAC 正方向：反式肽键对应 \(\omega=\pm\pi\)，顺式对应 \(\omega=0\)。式中 \(u\) 取 \(A-B\)；若取 \(B-A\)，则全部二面角整体平移 \(\pi\)。

**（三）主链二面角。** 依上定义，

$$
\varphi_i=\operatorname{Dih}(C_{i-1},N_i,C\alpha_i,C_i),\qquad
\psi_i=\operatorname{Dih}(N_i,C\alpha_i,C_i,N_{i+1}),
$$

$$
\omega_i=\operatorname{Dih}(C\alpha_i,C_i,N_{i+1},C\alpha_{i+1}),
$$

其存在范围如下表。

| 角 | 存在范围 | 个数 | 说明 |
|---|---|---|---|
| \(\varphi\) | \(i=1..n+1\) | \(n+1\) | \(\varphi_0\) 需 \(C_{-1}\)，不在模型内 |
| \(\psi\) | \(i=0..n\) | \(n+1\) | \(\psi_{n+1}\) 需 \(N_{n+2}\)，不在模型内 |
| \(\omega\) | \(i=0..n\) | \(n+1\) | 由序列确定，见 §4.2 |

---

## 2. 基于肽键为单位的片层

### 2.1 第 \(i\) 位的原子与链式方向

以肽键为基本单位，第 \(i\) 个片层定义为

$$
S_i=(C\alpha_i,\ C_i,\ N_{i+1},\ C\alpha_{i+1}),\qquad i=0..n,
$$

共 \(n+1\) 个，其中 \(S_0\) 与 \(S_n\) 为两个接入片层。其链式方向为

$$
C\alpha_i\ \longrightarrow\ C_i\ \longrightarrow\ N_{i+1}\ \longrightarrow\ C\alpha_{i+1}.
$$

记 \(T_{n,i}=N_{i+1}\)，\(T_{c,i}=C\alpha_{i+1}\)，则上述链式方向即

$$
C_i\ -\ T_{n,i}\ -\ T_{c,i}\ \longrightarrow\ C_{i+1}.
$$

片层内部唯一的扭转是 \(\omega_i\)；\(\omega_i\) 一经确定，\(S_i\) 的六个内部坐标全部确定，故 \(S_i\) 为刚体。片层的内部几何由三条键长 \(\ell^{C\alpha C}_i,\ell^{CN}_i,\ell^{N C\alpha}_i\) 与两个内角 \(\eta_i,\theta_i\) 给定。

### 2.2 二面角 \(\varphi(i)\to\psi(i+1)\)

片层 \(S_i\) 自身不转动，其前后各有一个铰链，如下表。

| 位置 | 铰链 | 绕行键 | 归属 |
|---|---|---|---|
| 左端 | \(\psi_i\) | \(C\alpha_i\to C_i\) | 残基 \(i\) |
| 内部 | \(\omega_i\) | \(C_i\to N_{i+1}\) | 肽键，已知 |
| 右端 | \(\varphi_{i+1}\) | \(N_{i+1}\to C\alpha_{i+1}\) | 残基 \(i+1\) |

因此，片层 \(S_i\) 所联系的一对 Ramachandran 角为 \((\varphi_{i+1},\psi_{i+1})\)；记作 \(\varphi(i)\to\psi(i+1)\)。

### 2.3 索引约定

本文中片层 \(S_i\) 与残基 \(i\) 共用同一下标，\(S_i\) 的首原子为 \(C\alpha_i\)。片层 \(i\) 自 \(C_i\) 起始、至 \(C_{i+1}\) 终止，故

$$
\varphi(i)=\varphi_{i+1},\qquad \psi(i+1)=\psi_{i+1},
$$

两者同属残基 \(i+1\)。相应地，

$$
P\big(\varphi(i),\psi(i+1)\mid m_{i-1},m_i,m_{i+1}\big)
=
P\big(\varphi_{i+1},\psi_{i+1}\mid m_{i-1},m_i,m_{i+1}\big).
$$

下文统一采用残基编号，即记作 \(P(\varphi_i,\psi_i\mid m_{i-1},m_i,m_{i+1})\)；两种写法仅差一个整体平移。

---

## 3. 运动学递推

### 3.1 主碳标架

对每个主碳定义规范标架

$$
p_i=C\alpha_i,\qquad
\hat z_i=\frac{C_i-C\alpha_i}{|C_i-C\alpha_i|},
$$

$$
\hat x_i=\frac{(N_i-C\alpha_i)_\perp}{\|(N_i-C\alpha_i)_\perp\|},
\qquad
\hat y_i=\hat z_i\times\hat x_i,
\qquad
T_i=[\,\hat x_i\ \ \hat y_i\ \ \hat z_i\ ;\ p_i\,]\in SE(3).
$$

其中 \(\hat z\) 为出键 \(C\alpha_i\to C_i\) 的方向，\(\hat x\) 为入键 \(N_i\to C\alpha_i\) 在垂直于 \(\hat z\) 的平面上的投影方向。端点 \(i=0,n+1\) 分别取 \(N_0\)、\(C_{n+1}\)，由 \(T_0^*,T_{n+1}^*\) 提供，无需另行选定方向。

由构造可知，\(T_i\) 与 \(\psi_i,\omega_i,\varphi_{i+1}\) 无关，仅依赖上游（含 \(\varphi_i\)）。

### 3.2 单步连杆

$$
T_{i+1}=T_i\,L_i(\psi_i,\varphi_{i+1}),\qquad i=0..n,
$$

$$
L_i=B_i^{0}\,\mathrm{Rz}(\psi_i)\,B_i^{1}\,\mathrm{Rz}(\varphi_{i+1})\,B_i^{2},
$$

其中 \(\omega_i\) 已并入 \(B_i^{1}\)，三个常数矩阵 \(B_i^{0},B_i^{1},B_i^{2}\in SE(3)\) 每个残基预计算一次。其显式形式为

$$
L_i=
B_i^{0}\,
\mathrm{Rz}(\psi_i)\,\mathrm{Tz}(\ell^{C\alpha C}_i)\,\mathrm{Rx}(\pi-\eta_i)\,
\mathrm{Rz}(\omega_i)\,\mathrm{Tz}(\ell^{CN}_i)\,\mathrm{Rx}(\pi-\theta_i)\,
\mathrm{Rz}(\varphi_{i+1})\,\mathrm{Tz}(\ell^{N C\alpha}_i)\,
B_i^{3},
$$

其中 \(B_i^{0}\) 把主碳标架变换至 \(\psi\) 铰链标架，\(B_i^{3}\) 把 \(\varphi\) 铰链标架变换至主碳标架（含 \(\kappa_{i+1}\) 所定的一转）。

铰链标架的 \(x\) 轴取在 IUPAC 参考方向 \((A-B)_\perp\) 上，故 \(\mathrm{Rz}\) 的转角即为 §1.4 所定义的二面角。

### 3.3 最终式

$$
F(\Theta)=T_0^*\prod_{i=0}^{n}L_i(\psi_i,\varphi_{i+1})=T_{n+1}\in SE(3).
$$

### 3.4 端点条件

入口位姿取 \(T_0=T_0^*\)；出口位姿要求 \(T_{n+1}=T_{n+1}^*\)。端点仅以上述两式进入模型，不引入额外变量。

### 3.5 计算次序

依次预计算 \(B_i^{0},B_i^{1},B_i^{2}\)，随后按 \(i=0..n\) 递推 \(T_{i+1}=T_iL_i\)，最终输出 \(T_{n+1}\)。内层循环仅含 \(n+1\) 次 \(4\times4\) 矩阵乘法，不涉及反解二面角。

---

## 4. 主要自由度约束

### 4.1 二面角概率约束（逐残基二维联合）

概率约束的对象是**同一主碳上的两个二面角**。对残基 \(i\)，\(\varphi_i\) 绕 \(N_i\!-\!C\alpha_i\)，\(\psi_i\) 绕 \(C\alpha_i\!-\!C_i\)，二者共享顶点 \(C\alpha_i\)；其统计规律由相邻三个残基的标签共同决定，记为

$$
P(\varphi_i,\psi_i\mid m_{i-1},m_i,m_{i+1}).
$$

该记号是 \(\mathbb T^2\) 上的**联合**密度：逗号表示联合，不是 \(P(\varphi_i)\,P(\psi_i)\)；它也不是肽键 \(\omega_i\) 的分布（\(\omega_i\) 由 §4.2 确定，不进入概率）。它描述的是一个主碳上的两个二面角受前后残基标签影响而呈现的概率规律。

现以完成概率分布函数
已知m_{i-1},m_i,m_{i+1}可获得任意φ_i,ψ_i对的概率

先验密度取

$$
p(\Theta)=P_0(\psi_0\mid m_0,m_1)\prod_{i=1}^{n}P(\varphi_i,\psi_i\mid m_{i-1},m_i,m_{i+1})\,P_{n+1}(\varphi_{n+1}\mid m_n,m_{n+1}),
$$

上式含一条**跨残基条件独立假设**：给定标签序列，各残基的角对相互独立。若要保留相邻残基角之间的耦合，应把乘积替换为马尔可夫链

$$
P(\varphi_1,\psi_1\mid m_0,m_1,m_2)\prod_{i=2}^{n}P(\varphi_i,\psi_i\mid \varphi_{i-1},\psi_{i-1},m_{i-1},m_i,m_{i+1});
$$

此时因子仍相乘（链式法则），改变的只是条件项。

边界项由三标签表对缺失一侧的标签按背景频率 \(\pi\) 边缘化、再对另一角取边缘得到：

$$
P_0(\psi_0\mid m_0,m_1)=\sum_{a}\pi(a)\int P(\varphi,\psi_0\mid a,m_0,m_1)\,d\varphi,
$$

$$
P_{n+1}(\varphi_{n+1}\mid m_n,m_{n+1})=\sum_{b}\pi(b)\int P(\varphi_{n+1},\psi\mid m_n,m_{n+1},b)\,d\psi .
$$

端点的标签窗口不完整，故边界表是截断与边缘的结果，与内部表并非同一对象；若直接估计两标签表，则以上两式取直接估计值。

各因子在自身变量上归一，故 \(p\) 在 \(\mathbb T^{2n+2}\) 上已归一，无需额外的归一化常数。若仅有单变量表，可按链式法则 \(P(\varphi,\psi)=P(\varphi)P(\psi\mid\varphi)\) 构造联合密度；这是精确分解，不引入独立性假设。因子数与未知量数相等：

$$
1+2n+1=2n+2 .
$$

### 4.2 肽键转角的映射约束

$$
\omega_i=\Omega(m_i,m_{i+1})\in\{0,\pi\},\qquad i=0..n,
$$

其中 \(|\omega_i|=\pi\) 对应反式，\(\omega_i=0\) 对应顺式。由于整条序列已知，\(\omega_i\) 全部为已知量：它们既不是未知量，也不进入积分，无需引入 \(Q\) 因子。此约束属参数消去，与 §4.1 的概率约束性质不同。

若需保留顺反不确定性，可将 \(\omega_i\in\{0,\pi\}\) 视为离散隐变量，对 \(2^{n+1}\) 种模式求和；模型其余部分不变。

### 4.3 约束的分类

| 类型 | 对象 | 数学形式 |
|---|---|---|
| 参数消去 | \(\omega_0..\omega_n\) | \(\omega_i=\Omega(m_i,m_{i+1})\) |
| 概率约束 | \(\psi_0,\ \{\varphi_i,\psi_i\}_{i=1}^{n},\ \varphi_{n+1}\) | 每个残基一个 \(\mathbb T^2\) 联合密度，总密度 \(p(\Theta)\) |
| 硬约束 | 末端位姿及可选的端点角 | \(\mathbf 1_\Omega\) 或 \(\delta(F-T^*)\) |

---

## 5. 场景描述

设前段序列及其结构已知，后段序列及其结构已知，两端点在空间中的位姿分别为 \(T_0^*\) 与 \(T_{n+1}^*\)；中间序列 \(m_1,\dots,m_n\) 仅已知标签顺序。要求将中间序列接入两端，即确定满足 \(F(\Theta)=T_{n+1}^*\) 的构型 \(\Theta\)。

在此设定下，期望为氨基酸编码器提供终止依据的上下文

简化目标为问题

其一，可达性判定。对给定的 \(T_0^*,T_{n+1}^*\) 与序列，判定是否存在 \(\Theta\in\mathbb T^{2n+2}\) 使 \(F(\Theta)=T_{n+1}^*\)，即判定闭链约束集是否非空；该判定在何种条件下可以给出确定的结论。

其二，方向相容性。端点携带上下文，其半二面角 \(\psi_0\) 与 \(\varphi_{n+1}\) 的参考已由上下文确定；当这些端角另受数值或分布约束时，如何判定存在与端点二面角相容的接入构型，即附加端角条件后的闭链约束集是否非空。只做可接入可行性判断,不获得确切二面角序列.

其三，能量指标。能否以概率密度代替能量，为已接入的构型给出一个标量评价，用以刻画其紧绷程度或能量高低。

其四，统一指标。能否构造一个（或一组）指标，同时反映序列对两端的可达性，以及可达之后构型的紧绷程度或能量大小。

---

## 6. 问题简化

### 6.1 \(\omega\) 已知

序列全知，故 \(\omega_i=\Omega(m_i,m_{i+1})\) 全知。每个片层因而是刚体，模型只需确定各主链二面角；未知量数目由 \(3n+3\) 降为 \(2n+2\)。

### 6.2 端点肽键的补充

两端各有一个半二面角，即 \(\psi_0\) 与 \(\varphi_{n+1}\)。它们并非新的对象，而是端键上的普通二面角：

- 其参考原子由端点帽 \(T_0^*,T_{n+1}^*\) 给出；
- 若上下文结构给出其数值，则作为硬约束；
- 若仅给出分布，则采用 §4.1 的 \(P_0\) 与 \(P_{n+1}\)；
- 端部肽键 \(\omega_0,\omega_n\) 本身已由 §4.2 确定。

因此，接入问题是在两条闭链硬约束与两个端角的概率约束之下，确定中间序列的二面角。

### 6.3 变量集

$$
\Theta=\big(\underbrace{\psi_0}_{\text{N 端半二面角}},\
\underbrace{\{\varphi_i,\psi_i\}_{i=1}^{n}}_{\text{中间序列}},\
\underbrace{\varphi_{n+1}}_{\text{C 端半二面角}}\big)\in\mathbb T^{2n+2}.
$$

注：C 端边界角按残基编号为 \(\varphi_{n+1}\)；\(\psi_{n+1}\) 需 \(N_{n+2}\)，不在模型范围内。

给定 \(\Theta\) 与 \(T_0^*\)，整条序列的空间构型（含全部非主碳原子）完全确定。

---

## 附录 A　计数与维度

| 量 | 值 |
|---|---|
| 原子 | \(3n+6\) |
| 键 / 内角 / 扭转 | \(3n+5\) / \(3n+4\) / \(3n+3\) |
| 片层数 | \(n+1\) |
| 肽键角（由序列确定） | \(n+1\) |
| 未知量 | \(2n+2\) |
| 闭链约束 | \(6\) |
| \(\dim\Omega\) | \(2n-4\) |

| \(n\) | \(2n+2\) | \(\dim\Omega\) | 情形 |
|---|---|---|---|
| 1 | 4 | \(-2\) | 一般无解 |
| 2 | 6 | \(0\) | 离散解 |
| 3 | 8 | \(2\) | 解流形 |
| 4 | 10 | \(4\) | 解流形 |

## 附录 B　计算实现

1. 常数矩阵 \(B_i\) 预计算；\(T_i\) 为递推缓存，不作为自由变量。
2. 采样在 \(\Omega_\varepsilon\) 上进行，以 \(p(\Theta)\) 为提议密度。
3. 梯度可由串链 Jacobian 解析给出，无需数值差分。
4. \(T_0..T_{n+1}\) 可作为区间分支定界的检查点。
5. 二面角符号只需标定一次：取一已知构象核对 \(\operatorname{Dih}\) 与递推转角同号，此后固定。
6. 手性由 \(\hat y=\hat z\times\hat x\) 固定；若改为 \(\hat x\times\hat z\)，则全部 \(\psi,\varphi\) 反号。
