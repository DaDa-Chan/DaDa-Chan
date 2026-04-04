---
title: "Ch5 Hypothesis Testing"
date: 2025-07-11
category: "Mathematical Statistics"
excerpt: "Neyman-Pearson framework, normal population tests, exponential family tests, likelihood ratio test, and p-values."
---

## §5.1 假设检验的若干基本概念

### 一、假设检验相关概念

#### 假设检验的做法:
1.  **建立假设**:
    * **原假设 (或称零假设) $H_0$**: 所要检验的假设，通常是保守的、不易被否定的命题。
    * **备择假设 (或称对立假设) $H_1$**: 与$H_0$不相容的假设，当$H_0$被拒绝时所选择的假设。

2.  **构造检验统计量**:
    * 先假定$H_0$为真，然后构造合适的统计量（检验统计量）来判断其真伪。

3.  **确定拒绝域与接受域**:
    * 将样本空间$\Omega$划分为两个互斥的集合 $D$ 和 $\overline{D}$。
    * **拒绝域 $D$**: 若样本观测值 $\tilde{x} \in D$，则拒绝$H_0$。
    * **接受域 $\overline{D}$**: 若样本观测值 $\tilde{x} \in \overline{D}$，则接受$H_0$。

4.  **据给定的显著性水平来确定临界点**:
    * **第一类错误**: 假设$H_0$为真，但我们拒绝了它。这个错误的概率是我们希望控制的风险。
    * **显著性水平 $\alpha$**: 我们要求犯第一类错误的概率不超过$\alpha$。

5.  **做判断**:
    * **实际推断原理**: 小概率事件在一次实验中是不会发生的。
    * 检验准则：若样本观测值导致了一个在$H_0$为真前提下的小概率事件发生，我们就有理由怀疑$H_0$的真实性，从而拒绝它。

### 二、假设

**定义 (假设)**: 指的是关于总体分布的命题。
* **原假设 (Null Hypothesis)**: $H_0$
* **备择假设 (Alternative Hypothesis)**: $H_1$

**定义 (参数假设)**: 在参数分布族 $\mathcal{F}=\{p_{\theta}(x):\theta\in\Theta\}$ 场合所作的关于其未知参数$\theta$的假设。
$$H_{0}:\theta\in\Theta_{0} \quad \text{vs} \quad H_{1}:\theta\in\Theta_{1}$$
其中，$\Theta_{0}$ 与 $\Theta_{1}$ 是参数空间 $\Theta$ 的两个互不相交的非空子集。
* **简单假设**: 假设中只含一个元素。
* **复杂假设**: 假设中包含多个元素。

### 三、检验

**定义 (检验)**: 给出一个规则，凭此规则，在有了样本观测值之后，就可以作出接受还是拒绝原假设$H_0$的判断。
* **检验函数 $\varphi(\tilde{x})$**:
    $$\varphi(\tilde{x})=\begin{cases}1, & \text{若 }\tilde{x}\in D \\ 0, & \text{若 }\tilde{x}\in\overline{D} \end{cases}$$
* **检验统计量**: 能从样本空间中划分出拒绝域的统计量。

### 四、两类错误与势函数

* **第一类错误 (Type I Error)**: 原假设$H_0$为真，但被拒绝（弃真）。其发生的概率为：
    $$\alpha(\theta)=P_{\theta}(\tilde{X}\in D), \quad \theta\in\Theta_{0}$$
* **第二类错误 (Type II Error)**: 原假设$H_0$为假，但被接受（取伪）。其发生的概率为：
    $$\beta(\theta)=P_{\theta}(\tilde{X}\in\overline{D}), \quad \theta\in\Theta_{1}$$
* **注**: 当样本量n固定时，$\alpha$和$\beta$此消彼长，不能同时减小。

**定义 5.1.2 (势函数)**: 设检验的拒绝域为D。称
$$g(\theta)=P_{\theta}(\tilde{X}\in D), \quad \theta\in\Theta$$
为此检验的**势函数** (power function)，也称功效函数。
$$g(\theta) = \begin{cases} \text{犯第一类错误的概率}, & \text{当} \theta \in \Theta_0 \\ 1-\text{犯第二类错误的概率}, & \text{当} \theta \in \Theta_1 \end{cases}$$

### 五、Neyman-Pearson原则与显著性水平
**Neyman-Pearson (NP)原则**: 控制犯第一类错误的概率，即选定一个数$\alpha, (0<\alpha<1)$，要求检验犯第一类错误的概率不超过$\alpha$：
$$sup_{\theta\in\Theta_{0}}\alpha(\theta)\le\alpha$$
* **显著性水平为$\alpha$的检验**: 犯第一类错误的概率不超过$\alpha$的检验。
* 通常取检验临界点使得 $sup_{\theta\in\Theta_{0}}\alpha(\theta)=\alpha$。

### 六、处理假设检验问题的一般步骤
1.  根据问题提出原假设$H_0$和备择假设$H_1$。
2.  确定检验统计量$T(X_1, \dots, X_n)$，并根据备择假设的形式确定拒绝域D的形式。
    * **单侧拒绝域**: $D=\{\tilde{x}:T(\tilde{x})<C\}$ 或 $D=\{\tilde{x}:T(\tilde{x})>C\}$
    * **双侧拒绝域**: $D=\{\tilde{x}:T(\tilde{x})<C_1 \text{ 或 } T(\tilde{x})>C_2\}$ 或 $D=\{\tilde{x}:|T(\tilde{x})|>C\}$
3.  选取适当的显著性水平$\alpha$，并求出临界值，使得 $sup P(\tilde{X}\in D|H_{0} \text{为真})\le\alpha$ 且尽可能地接近$\alpha$。
4.  由样本的观察值算出检验统计量的值 $t=T(\tilde{x})$，并与临界点进行比较，若观察值落入拒绝域D，则拒绝原假设$H_0$，否则接受原假设。

---

## §5.2 正态分布总体参数的假设检验

### 一、单个正态总体

设总体 $X\sim N(\mu,\sigma^{2})$，$\tilde{X}=(X_{1},X_{2},\cdots,X_{n})$ 是取自该总体的样本，并记
$$\overline{X}=\frac{1}{n}\sum_{i=1}^{n}X_{i}, \quad S^{2}=\frac{1}{n-1}\sum_{i=1}^{n}(X_{i}-\overline{X})^{2}$$

#### 1、正态总体均值的假设检验

#### (A) $\sigma^{2}$ 已知时
对于假设 $H_{0}:\mu=\mu_{0} \leftrightarrow H_{1}:\mu>\mu_{0}$，检验统计量为：
$$u=u(\tilde{X})=\frac{\overline{X}-\mu_{0}}{\sigma/\sqrt{n}}$$
在 $H_0$ 为真时，$u \sim N(0,1)$。
* 对于 $H_{0}:\mu=\mu_{0} \leftrightarrow H_{1}:\mu>\mu_{0}$ 或 $H_{0}:\mu\le\mu_{0} \leftrightarrow H_{1}:\mu>\mu_{0}$，拒绝域为 $D=\{u>u_{\alpha}\}$。
* 对于 $H_{0}:\mu=\mu_{0} \leftrightarrow H_{1}:\mu<\mu_{0}$ 或 $H_{0}:\mu\ge\mu_{0} \leftrightarrow H_{1}:\mu<\mu_{0}$，拒绝域为 $D=\{u<-u_{\alpha}\}$。
* 对于 $H_{0}:\mu=\mu_{0} \leftrightarrow H_{1}:\mu\ne\mu_{0}$，拒绝域为 $D=\{|u|>u_{\alpha/2}\}$。
这些检验称为 **u-检验法**。

#### (B) $\sigma^{2}$ 未知时
用$\sigma^2$的无偏估计 $S^2$ 代替 $\sigma^2$，得检验统计量为：
$$t=t(\tilde{X})=\frac{\overline{X}-\mu_{0}}{S/\sqrt{n}}$$
**引理**: 设随机变量T服从自由度为n，非中心参数为$\delta$的非中心t分布 $t(n,\delta)$，则对任意常数C，$P_{\delta}(T<C)$ 是$\delta$的严格单调减函数。

**结论**: $t(\tilde{X})=\frac{\sqrt{n}(\overline{X}-\mu_{0})/\sigma}{S/\sigma} \sim t(n-1,\frac{\mu-\mu_{0}}{\sigma/\sqrt{n}})$。
* 对于 $H_{0}:\mu\ge\mu_{0}\leftrightarrow H_{1}:\mu<\mu_{0}$，拒绝域为 $D=\{t<-t_{\alpha}(n-1)\}$。
* 对于 $H_{0}:\mu\le\mu_{0}\leftrightarrow H_{1}:\mu>\mu_{0}$，拒绝域为 $D=\{t>t_{\alpha}(n-1)\}$。
* 对于 $H_{0}:\mu=\mu_{0}\leftrightarrow H_{1}:\mu\ne\mu_{0}$，拒绝域为 $D=\{|t|>t_{\alpha/2}(n-1)\}$。
这些检验称为 **t-检验**。

#### 2、正态总体方差 $\sigma^2$ 的假设检验

##### (A) $\mu$ 已知时
取检验统计量为：
$$\chi^{2}=\chi^{2}(\tilde{X})=\frac{\sum_{i=1}^{n}(X_{i}-\mu)^{2}}{\sigma_{0}^{2}}$$
在 $H_0: \sigma^2 = \sigma_0^2$ 为真时，$\chi^2 \sim \chi^2(n)$。
* 对于 $H_{0}:\sigma^{2}\ge\sigma_{0}^{2}\leftrightarrow H_{1}:\sigma^{2}<\sigma_{0}^{2}$，拒绝域为 $D=\{\chi^{2}<\chi_{1-\alpha}^{2}(n)\}$。
* 对于 $H_{0}:\sigma^{2}\le\sigma_{0}^{2}\leftrightarrow H_{1}:\sigma^{2}>\sigma_{0}^{2}$，拒绝域为 $D=\{\chi^{2}>\chi_{\alpha}^{2}(n)\}$。
* 对于 $H_{0}:\sigma^{2}=\sigma_{0}^{2}\leftrightarrow H_{1}:\sigma^{2}\ne\sigma_{0}^{2}$，拒绝域为 $D=\{\chi^{2}<\chi_{1-\alpha/2}^{2}(n) \text{ 或 } \chi^{2}>\chi_{\alpha/2}^{2}(n)\}$。

##### (B) $\mu$ 未知时
取检验统计量为：
$$\chi^{2}=\chi^{2}(\tilde{X})=\frac{\sum_{i=1}^{n}(X_{i}-\overline{X})^{2}}{\sigma_{0}^{2}}=\frac{(n-1)S^{2}}{\sigma_{0}^{2}}$$
在 $H_0: \sigma^2 = \sigma_0^2$ 为真时，$\chi^2 \sim \chi^2(n-1)$。
* 对于 $H_{0}:\sigma^{2}\ge\sigma_{0}^{2}\leftrightarrow H_{1}:\sigma^{2}<\sigma_{0}^{2}$，拒绝域为 $D=\{\chi^{2}<\chi_{1-\alpha}^{2}(n-1)\}$。
* 对于 $H_{0}:\sigma^{2}\le\sigma_{0}^{2}\leftrightarrow H_{1}:\sigma^{2}>\sigma_{0}^{2}$，拒绝域为 $D=\{\chi^{2}>\chi_{\alpha}^{2}(n-1)\}$。
* 对于 $H_{0}:\sigma^{2}=\sigma_{0}^{2}\leftrightarrow H_{1}:\sigma^{2}\ne\sigma_{0}^{2}$，拒绝域为 $D=\{\chi^{2}<\chi_{1-\alpha/2}^{2}(n-1) \text{ 或 } \chi^{2}>\chi_{\alpha/2}^{2}(n-1)\}$。

---

### 二、两个正态总体

设 $X_1, \dots, X_m$ 来自 $N(\mu_X, \sigma_X^2)$，$Y_1, \dots, Y_n$ 来自 $N(\mu_Y, \sigma_Y^2)$，两样本相互独立。

#### 1、比较均值的假设检验

##### (A) $\sigma_X^2$ 和 $\sigma_Y^2$ 均已知时
检验统计量为：
$$U=U(\tilde{X},\tilde{Y})=\frac{\overline{X}-\overline{Y}}{\sqrt{\sigma_{X}^{2}/m+\sigma_{Y}^{2}/n}}$$
在 $H_0: \mu_X = \mu_Y$ 为真时，$U \sim N(0,1)$。
* 对于 $H_{0}:\mu_{X}\ge\mu_{Y}\leftrightarrow H_{1}:\mu_{X}<\mu_{Y}$，拒绝域为 $D=\{U<-u_{\alpha}\}$。
* 对于 $H_{0}:\mu_{X}\le\mu_{Y}\leftrightarrow H_{1}:\mu_{X}>\mu_{Y}$，拒绝域为 $D=\{U>u_{\alpha}\}$。
* 对于 $H_{0}:\mu_{X}=\mu_{Y}\leftrightarrow H_{1}:\mu_{X}\ne\mu_{Y}$，拒绝域为 $D=\{|U|>u_{\alpha/2}\}$。

##### (B) $\sigma_X^2$ 和 $\sigma_Y^2$ 均未知时
* **(a) $\sigma_{X}^{2}=\sigma_{Y}^{2}=\sigma^{2}$ 情形**: 检验统计量为
    $$t=\frac{\overline{X}-\overline{Y}}{S_{W}\sqrt{1/m+1/n}}$$
    其中 $S_{W}^{2}=\frac{(m-1)S_{X}^{2}+(n-1)S_{Y}^{2}}{m+n-2}$。在 $H_0: \mu_X = \mu_Y$ 为真时，$t \sim t(n+m-2)$。
* **(b) m=n 情形 (成对数据)**: 令 $Z_i = X_i - Y_i$，问题转化为单样本t检验。检验统计量为
    $$t=t(\tilde{Z})=\frac{\overline{Z}}{S_{Z}^{*}/\sqrt{n}}$$
    在 $H_0: \mu_X = \mu_Y$ (即 $\mu_Z=0$) 为真时，$t \sim t(n-1)$。
* **(c) 一般情形**:
    * **大样本**: 当m和n都很大时，检验统计量 $U=\frac{\overline{X}-\overline{Y}}{\sqrt{S_{X}^{2}/m+S_{Y}^{2}/n}}$ 近似服从 $N(0,1)$。
    * **小样本**: 检验统计量 $\frac{\overline{X}-\overline{Y}-(\mu_{X}-\mu_{Y})}{\sqrt{S_{X}^{2}/m+S_{Y}^{2}/n}}$ 近似服从 $t(l)$ 分布，其中
        $$l=\frac{\{{S_{X}^{*}}^{2}/m+{S_{Y}^{*}}^{2}/n\}^{2}}{\frac{{S_{X}^{*}}^{4}}{m^{2}(m-1)}+\frac{{S_{Y}^{*}}^{4}}{n^{2}(n-1)}}$$

#### 2、比较方差的假设检验
##### (A) $\mu_X$ 和 $\mu_Y$ 均已知时
检验统计量为：
$$F=\frac{\sum_{i=1}^{m}(X_{i}-\mu_{X})^{2}/m}{\sum_{j=1}^{n}(Y_{j}-\mu_{Y})^{2}/n}$$
在 $H_0: \sigma_X^2 = \sigma_Y^2$ 为真时，$F \sim F(m,n)$。
* 对于 $H_{0}:\sigma_{X}^{2}\ge\sigma_{Y}^{2}\leftrightarrow H_{1}:\sigma_{X}^{2}<\sigma_{Y}^{2}$，拒绝域为 $D=\{F<F_{1-\alpha}(m,n)\}$。
* 对于 $H_{0}:\sigma_{X}^{2}\le\sigma_{Y}^{2}\leftrightarrow H_{1}:\sigma_{X}^{2}>\sigma_{Y}^{2}$，拒绝域为 $D=\{F>F_{\alpha}(m,n)\}$。
* 对于 $H_{0}:\sigma_{X}^{2}=\sigma_{Y}^{2}\leftrightarrow H_{1}:\sigma_{X}^{2}\ne\sigma_{Y}^{2}$，拒绝域为 $D=\{F<F_{1-\alpha/2}(m,n) \text{ 或 } F>F_{\alpha/2}(m,n)\}$。

##### (B) $\mu_X$ 和 $\mu_Y$ 均未知时
检验统计量为：
$$F=\frac{S_{X}^{2}}{S_{Y}^{2}}$$
在 $H_0: \sigma_X^2 = \sigma_Y^2$ 为真时，$F \sim F(m-1,n-1)$。
* 对于 $H_{0}:\sigma_{X}^{2}\ge\sigma_{Y}^{2}\leftrightarrow H_{1}:\sigma_{X}^{2}<\sigma_{Y}^{2}$，拒绝域为 $D=\{F<F_{1-\alpha}(m-1,n-1)\}$。
* 对于 $H_{0}:\sigma_{X}^{2}\le\sigma_{Y}^{2}\leftrightarrow H_{1}:\sigma_{X}^{2}>\sigma_{Y}^{2}$，拒绝域为 $D=\{F>F_{\alpha}(m-1,n-1)\}$。
* 对于 $H_{0}:\sigma_{X}^{2}=\sigma_{Y}^{2}\leftrightarrow H_{1}:\sigma_{X}^{2}\ne\sigma_{Y}^{2}$，拒绝域为 $D=\{F<F_{1-\alpha/2}(m-1,n-1) \text{ 或 } F>F_{\alpha/2}(m-1,n-1)\}$。

---

## §5.3 假设检验与区间估计

### 假设检验与置信区间的关系

**定理**:
* 设对每个 $\theta_{0}\in\Theta$，$A(\theta_{0})$ 为假设检验问题 $H_{0}:\theta=\theta_{0} \leftrightarrow H_{1}:\theta\ne\theta_{0}$ 的显著性水平为$\alpha$的检验的接受域。对每个$\tilde{x}$，令
    $$C(\tilde{x})=\{\theta_{0}:\tilde{x}\in A(\theta_{0})\}$$
    则随机集 $C(\tilde{X})$ 是参数$\theta$的置信水平为$1-\alpha$的置信集。

* 反过来，设 $C(\tilde{X})$ 是参数$\theta$的置信水平为$1-\alpha$的置信集。对任意的 $\theta_{0}\in\Theta$，令
    $$A(\theta_{0})=\{\tilde{x}:\theta_{0}\in C(\tilde{x})\}$$
    则 $A(\theta_{0})$ 为假设检验问题 $H_{0}:\theta=\theta_{0} \leftrightarrow H_{1}:\theta\ne\theta_{0}$ 的显著性水平为$\alpha$的检验的接受域。

---
### 单参数指数型分布族的假设检验

#### 一、单参数指数型分布族的性质

**定理**: 设总体X服从单参数指数型分布，$\tilde{X}=(X_{1},X_{2},\cdots,X_{n})$ 是取自该总体的样本，其联合密度函数(或分布列)为:
$$p(\tilde{x};\theta)=C^{*}(\theta)\cdot exp\{Q(\theta)\cdot U(\tilde{x})\}\cdot h^{*}(\tilde{x})$$
其中，$Q(\theta)$ 是$\theta$的严格增函数。若 $\psi(U(\tilde{x}))$ 是 $U(\tilde{x})$ 的一个非降函数，则 $E_{\theta}\psi(U(\tilde{X}))$ 也是$\theta$的一个非降函数。

**推论**: 设总体X服从单参数指数型分布，其联合密度函数(或分布列)为:
$$p(\tilde{x};\theta)=C^{*}(\theta)\cdot exp\{Q(\theta)\cdot U(\tilde{x})\}\cdot h^{*}(\tilde{x})$$
其中，$Q(\theta)$ 是$\theta$的严格增函数，则对任意给定的常数C，
$P_{\theta}\{U(\tilde{X})>C\}$ 和 $P_{\theta}\{U(\tilde{X})<C\}$ 分别是$\theta$的一个非降和非增函数。

#### 二、单参数指数型分布族的假设检验问题

设总体X服从单参数指数型分布，其联合密度函数(或分布列)为:
$$p(\tilde{x};\theta)=C^{*}(\theta)\cdot exp\{Q(\theta)\cdot U(\tilde{x})\}\cdot h^{*}(\tilde{x})$$
其中，$Q(\theta)$ 是$\theta$的严格增函数。

对于假设检验问题：
1.  $H_{0}:\theta\ge\theta_{0}\leftrightarrow H_{1}:\theta<\theta_{0}$
2.  $H_{0}:\theta\le\theta_{0}\leftrightarrow H_{1}:\theta>\theta_{0}$
3.  $H_{0}:\theta=\theta_{0}\leftrightarrow H_{1}:\theta\ne\theta_{0}$

常取检验统计量为 $U(\tilde{X})$。当 $E_{\theta}U(\tilde{X})$ 是$\theta$的单调非降函数时，上述三个假设检验问题的拒绝域形式分别为
* $\{U(\tilde{x})<C\}$
* $\{U(\tilde{x})>C\}$
* $\{U(\tilde{x})<C_{1} \text{ or } U(\tilde{x})>C_{2}\}$

拒绝域的临界值可以由 $\theta=\theta_{0}$ 时，$U(\tilde{X})$ 的分布求得。

#### 三、指数分布下参数假设检验问题
总体X服从参数为$\theta^{-1}$的指数分布 $(\theta>0)$，其联合密度函数为：
$$p(\tilde{x};\theta)=\theta^{-n}exp\{-\frac{1}{\theta}\cdot\sum_{i=1}^{n}x_{i}\}$$
检验统计量为 $T=\sum_{i=1}^{n}X_{i}=n\overline{X}$。
* 对于 $H_{0}:\theta\ge\theta_{0}\leftrightarrow H_{1}:\theta<\theta_{0}$，显著性水平为$\alpha$的检验的拒绝域为:
    $$\{\tilde{x}:2n\overline{x}/\theta_{0}<\chi_{1-\alpha}^{2}(2n)\}$$
* 对于 $H_{0}:\theta\le\theta_{0}\leftrightarrow H_{1}:\theta>\theta_{0}$，显著性水平为$\alpha$的检验的拒绝域为:
    $$\{\tilde{x}:2n\overline{x}/\theta_{0}>\chi_{\alpha}^{2}(2n)\}$$
* 对于 $H_{0}:\theta=\theta_{0}\leftrightarrow H_{1}:\theta\ne\theta_{0}$，显著性水平为$\alpha$的检验的拒绝域为:
    $$\{\tilde{x}:2n\overline{x}/\theta_{0}<\chi_{1-\alpha/2}^{2}(2n) \text{ 或 } 2n\overline{x}/\theta_{0}>\chi_{\alpha/2}^{2}(2n)\}$$

#### 四、两点分布下参数假设检验问题
设总体X服从二点分布 $B(1,p), 0<p<1$。检验统计量为 $U(\tilde{X})=\sum_{i=1}^{n}X_{i}$。
* $H_{0}:p\ge p_{0}\leftrightarrow H_{1}:p<p_{0}$，拒绝域为 $\{\tilde{x}:\sum x_{i}\le C^{\prime}\}$，其中 $C^{\prime}=max\{C\in\mathcal{N}:\sum_{i=0}^{C}\binom{n}{i}p_{0}^{i}(1-p_{0})^{n-i}\le\alpha\}$。
* $H_{0}:p\le p_{0}\leftrightarrow H_{1}:p>p_{0}$，拒绝域为 $\{\tilde{x}:\sum x_{i}\ge C^{\prime}\}$，其中 $C^{\prime}=min\{C\in\mathcal{N}:\sum_{i=C}^{n}\binom{n}{i}p_{0}^{i}(1-p_{0})^{n-i}\le\alpha\}$。
* $H_{0}:p=p_{0}\leftrightarrow H_{1}:p\ne p_{0}$，拒绝域为 $\{\tilde{x}:\sum x_{i}\le C_{1}^{\prime} \text{ 或 } \sum x_{i}\ge C_{2}^{\prime}\}$，其中 $C_{1}^{\prime}, C_{2}^{\prime}$ 按显著性水平$\alpha/2$确定。

#### 五、泊松分布下参数λ假设检验问题
$\tilde{X}=(X_{1},\dots,X_{n})$ 取自总体 $X\sim P(\lambda)$。检验统计量为 $T=\sum_{i=1}^{n}X_{i}$。
* 对于 $H_{0}:\lambda\ge\lambda_{0}\leftrightarrow H_{1}:\lambda<\lambda_{0}$，拒绝域为 $\{T\le C^{\prime}\}$，其中$C^{\prime}$为满足 $2n\lambda_{0}\ge\chi_{\alpha}^{2}(2(C+1))$ 的最大正整数。
* 对于 $H_{0}:\lambda\le\lambda_{0}\leftrightarrow H_{1}:\lambda>\lambda_{0}$，拒绝域为 $\{T\ge C^{\prime}\}$，其中$C^{\prime}$为满足 $2n\lambda_{0}\le\chi_{1-\alpha}^{2}(2C)$ 的最小正整数。
* 对于 $H_{0}:\lambda=\lambda_{0}\leftrightarrow H_{1}:\lambda\ne\lambda_{0}$，拒绝域为 $\{T\le C_{2}^{\prime} \text{ 或 } T\ge C_{1}^{\prime}\}$，其中$C_1', C_2'$ 按显著性水平$\alpha/2$确定。

当n充分大时，可取 $\frac{T-n\lambda_{0}}{\sqrt{n\lambda_{0}}}$ 为近似服从$N(0,1)$的检验统计量，并使用U-检验法。

## §5.4 (广义)似然比检验((Generalized) Likelihood Ratio Test)

#### 一、广义似然比检验

设 $\tilde{X}=(X_{1},\cdots,X_{n})$ 为来自总体 $X\sim p(x;\theta)$ 的样本。似然函数为
$$L(\theta;\tilde{x})=\prod_{i=1}^{n}p(x_{i};\theta), \quad \theta\in\Theta$$
考察假设检验问题
$$H_{0}:\theta\in\Theta_{0}\leftrightarrow H_{1}:\theta\in\Theta_{1}=\Theta_{0}^{c}$$

**定义 5.3.1**: 令
$$\lambda(\tilde{x})=\frac{sup_{\theta\in\Theta}L(\theta;\tilde{x})}{sup_{\theta\in\Theta_{0}}L(\theta;\tilde{x})}=\frac{L(\hat{\theta};\tilde{x})}{L(\hat{\theta_{0}};\tilde{x})}$$
其中
$$\hat{\theta_{0}}=arg~sup_{\theta\in\Theta_{0}}L(\theta;\tilde{x}), \quad \hat{\theta}=arg~sup_{\theta\in\Theta}L(\theta;\tilde{x})$$
称$\lambda(\tilde{x})$为 **(广义)似然比**((generalized) likelihood ratio)，而称 $\lambda(\tilde{X})$ 为检验问题的**(广义)似然比(检验)统计量**。拒绝域为 $\{\tilde{x}:\lambda(\tilde{x})>C\}$ 的检验称为**(广义)似然比检验**((generalized) likelihood ratio test)。

临界值C根据Neyman-Pearson原则确定，要求满足
$$sup_{\theta\in\Theta_{0}}P_{\theta}\{\lambda(\tilde{X})>C\}\le\alpha \quad \text{且尽可能接近}$$
有时，如果存在另一个统计量 $G(\tilde{X})$，使得 $\lambda(\tilde{X})$ 是 $G(\tilde{X})$ 的严格单调函数，并且 $G(\tilde{X})$ 的分布已知，那么可以据 $G(\tilde{X})$ 定出等价的拒绝域。

#### 二、分布的似然比检验

假设关于总体X的密度函数(或分布列)可提出如下两个假设:
$$H_{0}:p(x;\theta)=p_{0}(x;\theta)\leftrightarrow H_{1}:p(x;\theta)=p_{1}(x;\theta)$$
取样本 $\tilde{X}=(X_{1},\cdots,X_{n})$，似然比统计量为：
$$\lambda(\tilde{x})=\frac{L_{1}(\hat{\theta}_{1};\tilde{x})}{L_{0}(\hat{\theta}_{0};\tilde{x})}=\frac{sup_{\theta}\prod_{i=1}^{n}p_{1}(x_{i};\theta)}{sup_{\theta}\prod_{i=1}^{n}p_{0}(x_{i};\theta)}$$
检验规则为：
* 当 $\lambda(\tilde{x})>C$ 时, 拒绝 $H_0$；
* 当 $\lambda(\tilde{x})\le C$ 时, 接受 $H_0$。
其中，C满足 $P(\lambda(\tilde{X})>C|H_{0} \text{为真}) \le \alpha$ 且尽可能接近。

---
## §5.5 检验的p值 (p-value)

假设检验的结论通常是：在给定的显著性水平下，不是拒绝原假设就是保留原假设。然而，有时在不同的显著性水平下可能得到不同的结论。这种情况在应用中会带来一些麻烦。

### p值的定义
**定义 (Definition)**: 在统计假设检验中，p值是在原假设为真的前提下，获得一个至少与实际观测到的结果一样极端的结果的概率。(In statistical hypothesis testing, the p-value is the probability of obtaining a result at least as extreme as the one that was actually observed, assuming that the null hypothesis is true.)

通俗地说，p值是当原假设成立时，检验统计量取比观察到的结果更为极端的数值的概率。

p值是能用当前这组样本（的检验统计量的值）做出“拒绝$H_0$”这一结论的最小的显著性水平。

### p值的优点
* 它比较客观，避免了事先确定显著性水平。
* 由检验的p值与人们心目中的显著性水平$\alpha$进行比较，可以很容易作出检验的结论：
    * 如果 $\alpha \ge p$，则在显著性水平$\alpha$下拒绝$H_0$。
    * 如果 $\alpha < p$，则在显著性水平$\alpha$下保留$H_0$。

如今的统计软件如R/Python/SPSS/SAS等对检验问题一般都是对于给定的样本给出检验的p值。

### 处理假设检验问题的一般步骤
#### 方法一（传统方法）
1.  根据问题提出原假设$H_0$和备择假设$H_1$。
2.  确定检验统计量$T(X_1, \dots, X_n)$并根据备择假设确定拒绝域D的形式。
3.  选取适当的显著性水平$\alpha$，并求出临界值，使得 $sup P(\tilde{X} \in D | H_0 \text{为真}) \le \alpha$ 且尽可能地接近$\alpha$。
4.  由样本$\tilde{X}$的观察值算出检验统计量的值 $t=T(\tilde{x})$，并与临界点进行比较，若观察值落入拒绝域D，则拒绝原假设$H_0$，否则接受原假设。

#### 方法二（p值法）
1.  根据问题提出原假设$H_0$和备择假设$H_1$。
2.  确定检验统计量$T(X_1, \dots, X_n)$，并根据原假设和备择假设确定拒绝域D的形式。
3.  计算检验统计量的值，并由此得到p值。
4.  将p值与显著性水平$\alpha$进行比较，得出相应的结论。
