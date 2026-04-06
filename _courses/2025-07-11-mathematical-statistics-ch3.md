---
title: "Ch3 Point Estimation"
date: 2025-07-11
category: "Mathematical Statistics"
excerpt: "Unbiasedness, efficiency, consistency, method of moments, MLE, UMVUE, and Cramer-Rao lower bound."
---

## §3.1 引言

### 点估计的概念

**定义 3.1.1**: 用来估计未知参数$$\theta$$的统计量
$$\hat{\theta}=\hat{\theta}(X_{1},X_{2},\cdots,X_{n})$$
称为参数$$\theta$$的**点估计量** (point estimator)。而样本观察值的函数 $$\hat{\theta}=\hat{\theta}(x_{1},x_{2},\cdots,x_{n})$$ 称为参数的**估计值**。

### 评价估计量的基本准则

#### 1. 无偏性 (unbiasedness)
在平均意义下，$$\hat{\theta}$$的值不偏离$$\theta$$。

**定义 3.1.2**: 设 $$\hat{\theta}=\hat{\theta}(X_{1},X_{2},\cdots,X_{n})$$ 是参数$$\theta$$的一个估计量。假如
$$E(\hat{\theta})=\theta, \quad \forall\theta\in\Theta$$
则称$$\hat{\theta}$$为$$\theta$$的**无偏估计量** (unbiased estimator)。

* 通常 $$Bias_{\theta}\hat{\theta}=E_{\theta}(\hat{\theta})-\theta$$ 称为参数$$\theta$$的估计的**偏差**。无偏估计量没有系统偏差。

#### 渐近无偏性 (asymptotic unbiasedness)
**定义**: 设对每个自然数n，$$\hat{\theta_{n}}=\hat{\theta}(X_{1},X_{2},\cdots,X_{n})$$ 是$$\theta$$的估计量，假如
$$Bias_{\theta}(\hat{\theta}_{n})=E(\hat{\theta}_{n})-\theta\rightarrow0, \quad \text{as } n\rightarrow\infty, \quad \forall\theta\in\Theta$$
则称 $$\hat{\theta_{n}}$$ 是$$\theta$$的渐近无偏估计。

#### 2. 有效性 (efficiency)
**定义 3.1.3**: 设 $$\hat{\theta_{1}}=\hat{\theta_{1}}(X_{1},X_{2},\cdots,X_{n})$$ 和 $$\hat{\theta_{2}}=\hat{\theta_{2}}(X_{1},X_{2},\cdots,X_{n})$$ 是参数$$\theta$$的两个无偏估计，假如
$$Var_{\theta}(\hat{\theta_{1}})\le Var_{\theta}(\hat{\theta_{2}}), \quad \forall\theta\in\Theta$$
且至少存在一个 $$\theta\in\Theta$$，使得严格不等号成立，则称 $$\hat{\theta_{1}}$$ 比 $$\hat{\theta_{2}}$$ 有效。

#### 3. 均方误差 (Mean Squared Error, MSE)
**定义**: 设$$\hat{\theta}$$是参数$$\theta$$的一个估计量，其均方误差定义为
$$MSE_{\theta}(\hat{\theta})=E(\hat{\theta}-\theta)^{2}$$
它与方差和偏倚的关系是：
$$MSE_{\theta}(\hat{\theta})=Var(\hat{\theta})+Bias_{\theta}^{2}(\hat{\theta})$$

#### 4. 相合性 (Consistency)
**定义 3.1.4**: 设对每个自然数n，$$\hat{\theta_{n}}=\hat{\theta}(X_{1},X_{2},\cdots,X_{n})$$ 是$$\theta$$的估计量，假如
$$\hat{\theta}_{n}\rightarrow\theta, \quad n\rightarrow\infty, \quad \forall\theta\in\Theta$$
即对任意的 $$\epsilon>0$$，
$$P_{\theta}(|\hat{\theta}_{n}-\theta|\ge\epsilon)\rightarrow0, \quad n\rightarrow\infty, \quad \forall\theta\in\Theta$$
则称 $$\hat{\theta}_{n}$$ 是$$\theta$$的(弱)相合估计，亦称一致估计。

* 若 $$P_{\theta}(lim_{n\rightarrow\infty}\hat{\theta}_{n}=\theta)=1, \forall\theta\in\Theta$$，则称 $$\hat{\theta}_{n}$$ 是$$\theta$$的**强相合估计**。
* 若 $$lim_{n\rightarrow\infty}E_{\theta}|\hat{\theta}_{n}-\theta|^{r}=0, \forall\theta\in\Theta$$，则称$$\hat{\theta_{n}}$$是$$\theta$$的**r阶矩相合估计**，当 $$r=2$$ 时，称为**均方相合估计**。

#### 5. 渐近正态性 (asymptotic normality)
**定义**: 设 $$\hat{\theta}_{n}$$ 是$$\theta$$的估计量，假如
$$\sqrt{k_{n}}(\hat{\theta_{n}}-\theta)\Rightarrow N(0,\sigma^{2}), \quad \text{as } n\rightarrow\infty, \quad \forall\theta\in\Theta$$
则称 $$\hat{\theta_{n}}$$ 具有渐近正态性，$$\sigma^2$$称为渐近方差。

---
## §3.2 估计方法之一--矩法估计 (Method of Moments)

### 一、矩法的统计思想
用样本矩替换相应的总体矩。
* **总体k阶矩**: $$\mu_{k}=E(X^{k})=\int x^{k}dF(x)$$
* **样本k阶矩**: $$a_{n,k}=\int x^{k}dF_{n}(x)=\frac{1}{n}\sum_{i=1}^{n}X_{i}^{k}$$
* **总体k阶中心矩**: $$\nu_{k}=E(X-EX)^{k}$$
* **样本k阶中心矩**: $$m_{n,k}=\frac{1}{n}\sum_{i=1}^{n}(X_{i}-\overline{X})^{k}$$

**定义**: 设有总体分布族 $$\{F_{\theta};\theta\in\Theta\}$$，待估参数 $$g(\theta)=G(\mu_{1},\cdots,\mu_{k};\nu_{2},...,\nu_{s})$$。将式中的 $$\mu_i$$ 与 $$\nu_j$$ 分别用它们的矩估计量 $$a_{n,i}$$ 与 $$m_{n,j}$$ 代替，得
$$\hat{g(\theta)}=G(a_{n,1},\cdots,a_{n,k};m_{n,2},\cdots,m_{n,s})$$
则 $$\hat{g(\theta)}$$ 称为 $$g(\theta)$$ 的**矩法估计量** (moment estimator)。

### 二、分布中未知参数的矩法估计
#### 基本步骤:
1.  写出总体的前k阶矩，一般是这个待估参数的函数：$$\mu_{i}=g_{i}(\theta_{1},\theta_{2},\cdots,\theta_{k})$$
2.  写出待估参数关于总体矩的函数表达式：$$\theta_{j}=h_{j}(\mu_{1},\mu_{2},\cdots,\mu_{k})$$
3.  写出第(2)步中出现的总体矩相应的样本矩。
4.  以相应的样本矩替换第(2)步中出现的总体矩，从而得到参数的矩法估计量：$$\hat{\theta_{j}}=h_{j}(a_{n,1},a_{n,2},\cdots,a_{n,k})$$

#### 矩法估计的优缺点
优点：
- 简便且不要求事先知道总体分布
- 一般来说，矩法估计还是相应参数的相合估计量

缺点：
- 没有充分利用已知参数分布族提供的信息
- 在小样本场合无突出性质
- 不够稳健
- 在很多场合下，不具有唯一性

---

## §3.3 估计方法之二——极大似然估计法 (Maximum Likelihood Estimation)

### 一、似然函数
**定义**: 设 $$p(\tilde{x};\theta)$$ 为样本 $$\tilde{X}$$ 的联合密度函数(pdf)或分布列(pmf)，$$\theta\in\Theta$$ 是未知参数。对给定的$$\tilde{x}$$，参数$$\theta$$的函数
$$L(\theta;\tilde{x})=p(\tilde{x};\theta), \quad \theta\in\Theta$$
称为**似然函数** (likelihood function)。$$l(\theta;\tilde{x})=\log L(\theta;\tilde{x})$$ 称为**对数似然函数** (log likelihood function)。

### 二、极大似然估计 (MLE)

**定义**: 设总体分布族为 $$\mathcal{F}=\{p(x;\theta), \theta\in\Theta\}$$。似然函数为 $$L(\theta;\tilde{x})=\prod_{i=1}^{n}p(x_{i};\theta)$$。假如存在统计量 $$\hat{\theta}(\tilde{X})$$，使得
$$L(\hat{\theta}(\tilde{x});\tilde{x})=sup_{\theta\in\Theta}L(\theta;\tilde{x})$$
或等价地
$$l(\hat{\theta}(\tilde{x});\tilde{x})=sup_{\theta\in\Theta}l(\theta;\tilde{x})$$
则称 $$\hat{\theta}(\tilde{X})$$ 为$$\theta$$的**极大似然估计量** (Maximum likelihood estimator)，简记为MLE。

### 求MLE的基本步骤:
1.  写出似然函数 $$L(\theta;\tilde{x})$$。
2.  求 $$\hat{\theta}$$ 使得 $$L(\hat{\theta}(\tilde{x});\tilde{x})=sup_{\theta\in\Theta}L(\theta;\tilde{x})$$。
当似然函数关于参数可微时，求MLE常常转变为求下列方程的解：
* **似然方程**: $$\frac{\partial}{\partial\theta_{i}}L(\theta_{1},\theta_{2},\cdots,\theta_{k};\tilde{x})=0$$
* **对数似然方程**: $$\frac{\partial}{\partial\theta_{i}}l(\theta_{1},\theta_{2},\cdots,\theta_{k};\tilde{x})=0$$

### 极大似然估计的性质
* **不变原则**: 设 $$\hat{\theta}_{MLE}$$ 是$$\theta$$的MLE，则对任意一可测函数 $$g(\theta)$$，$$g(\hat{\theta}_{MLE})$$ 是 $$g(\theta)$$ 的MLE。
* **与充分统计量的关系**: 设 $$T(\tilde{X})$$ 是$$\theta$$的充分统计量，则$$\theta$$的MLE $$\hat{\theta}_{MLE}$$ 必是 $$T(\tilde{X})$$ 的函数。
* **渐近正态性**: 在适当的正则条件下，似然方程有一个解 $$\hat{\theta}_{n}$$ 满足：$$\hat{\theta}_{n}$$ 依概率收敛于真值 $$\theta_{0}$$，且
    $$\sqrt{n}(\hat{\theta_{n}}-\theta_{0})\longrightarrow N(0,\frac{1}{I(\theta_{0})}) \quad \text{当} \ n\rightarrow\infty$$
    其中 $$I(\theta)=E_{\theta}[(\frac{\partial \log p(X;\theta)}{\partial\theta})^{2}]$$ 为Fisher信息量。

## §3.4 一致最小方差无偏估计——优化准则

### 均方误差准则与最优无偏估计

**定义 3.4.1**: 设 $$\hat{g}=\hat{g}(\tilde{X})$$ 是 $$g(\theta)$$ 的一个估计量，称
$$E_{\theta}(\hat{g}-g(\theta))^{2}$$
为$$\hat{g}$$的**均方误差** (mean squared error)。

在所有估计量中，不存在一个在均方误差意义下的最优估计。因此，我们将寻找范围缩小到由无偏估计组成的类：
$$U_{g}=\{\hat{g}=\hat{g}(\tilde{X}):E_{\theta}\hat{g}=g(\theta),\forall\theta\in\Theta\}$$
当 $$U_g$$ 非空时，我们称参数 $$g(\theta)$$ 是**可估参数**。在 $$U_g$$ 中，统计量的均方误差就是它的方差。

**定义**: 设 $$g(\theta)$$ 是一个可估参数，$$U_g$$ 为 $$g(\theta)$$ 的无偏估计类。假如 $$\hat{g}^* = \hat{g}^*(\tilde{X})$$ 是这样一个无偏估计，对任何 $$\hat{g}(\tilde{X}) \in U_g$$，有
$$Var_{\theta}\{\hat{g}^{*}\}\le Var_{\theta}\{\hat{g}\}, \quad \forall\theta\in\Theta$$
则称 $$\hat{g}^*$$ 是 $$g(\theta)$$ 的最优无偏估计，也称**一致最小方差无偏估计** (uniform minimum variance unbiased estimator)，记为**UMVUE**。

### 寻找UMVUE的方法

#### Rao-Blackwell 定理

**定理 3.4.1 (Rao-Blackwell)**: 设 $$T=T(\tilde{X})$$ 是参数 $$\theta \in \Theta$$ 的充分统计量，$$\varphi(\tilde{X})$$ 是参数 $$g(\theta)$$ 的一个无偏估计，则
$$\hat{g}(T)=E\{\varphi(\tilde{X})|T\}$$
也是 $$g(\theta)$$ 的无偏估计，且
$$Var_{\theta}\{\hat{g}(T)\}\le Var_{\theta}\{\varphi(\tilde{X})\}, \quad \forall\theta\in\Theta$$
其中等号成立的充分必要条件是 $$P_{\theta}\{\varphi(\tilde{X})=\hat{g}(T)\}=1, \forall\theta\in\Theta$$。

**定理推论**:
* UMVUE一定是充分统计量的函数。
* 改善无偏估计方差的方法：对现有的无偏估计关于充分统计量求条件期望。

#### 零无偏估计法

**定理 3.4.2**: 设 $$\hat{g}=\hat{g}(\tilde{X})$$ 为 $$g(\theta)$$ 的一个无偏估计，$$Var_{\theta}(\hat{g})<\infty (\forall\theta\in\Theta)$$。记零无偏估计集
$$\mathcal{L}=\{l=l(\tilde{X}):E_{\theta}(l(\tilde{X}))=0,\forall\theta\in\Theta\}$$
若对 $$\forall l\in\mathcal{L}$$ 有
$$Cov_{\theta}(\hat{g},l)=E_{\theta}(\hat{g}\cdot l)=0, \quad \forall\theta\in\Theta$$
则 $$\hat{g}$$ 必为 $$g(\theta)$$ 的UMVUE。

**推论 3.4.1**: 设 $$T=T(\tilde{X})$$ 为$$\theta$$的充分统计量，$$h(T)$$是$$g(\theta)$$的一个无偏估计量。记
$$\mathcal{L}_{T}=\{\delta=\delta(T):E_{\theta}(\delta(T))=0,\forall\theta\in\Theta\}$$
若对 $$\forall \delta \in\mathcal{L}_{T}$$ 有
$$Cov_{\theta}(h(T),\delta(T))=E_{\theta}(h(T)\cdot\delta(T))=0, \quad \forall\theta\in\Theta$$
则$$h(T)$$必为$$g(\theta)$$的UMVUE。

#### 充分完备统计量法

##### 完备性 (Completeness)

**定义**: 设X的分布所在的分布族为 $$\mathcal{F}=\{p(x;\theta):\theta\in\Theta\}$$。假如对任一个实函数 $$\varphi(x)$$，由
$$E_{\theta}\varphi(X)=0, \quad \forall\theta\in\Theta$$
可推出
$$P_{\theta}\{\varphi(X)=0\}=1, \quad \forall\theta\in\Theta$$
则称此分布族 $$\mathcal{F}$$ 是**完备的**。

**定义**: 设有参数分布族 $$\mathcal{F}=\{p(x;\theta):\theta\in\Theta\}$$。T是一个统计量，假如对任一实函数 $$\varphi(t)$$，由
$$E_{\theta}\varphi(T)=0, \quad \forall\theta\in\Theta$$
可推出
$$P_{\theta}\{\varphi(T)=0\}=1, \quad \forall\theta\in\Theta$$
则称T是**完备统计量**。

##### Basu 定理

**定理 2.7.2 (Basu定理)**: 设总体X的分布族为 $$\mathcal{F}=\{p(x;\theta):\theta\in\Theta\}$$，$$T=T(\tilde{X})$$ 是一个完备且充分的统计量，若随机变量 $$V=V(\tilde{X})$$ 的分布与参数$$\theta$$无关，则对于任何$$\theta \in \Theta$$，T与V独立。

##### Lehmann-Scheffé 定理

**定理 3.4.3 (Lehmann-Scheffé)**: 设 $$S=S(\tilde{X})$$ 是参数的**充分完备统计量** (complete sufficient statistic)，则可估参数$$g(\theta)$$的UMVUE存在且唯一。若 $$\varphi(\tilde{X})$$ 是$$g(\theta)$$的一个无偏估计，则$$g(\theta)$$唯一的UMVUE就是
$$E[\varphi(\tilde{X})|S]$$

**寻找UMVUE的方法**:
1.  先找一个充分完备统计量S，再找一个可估参数的无偏估计$$\varphi(\tilde{X})$$，计算条件期望 $$E[\varphi(\tilde{X})|S]$$ 就得到UMVUE。
2.  先找一个充分完备统计量S，再找一个S的函数h(S)使得它是可估参数的无偏估计，h(S)即为所求的UMVUE。

**定理**: 设总体来自指数型分布族，样本的联合pdf或pmf为
$$p(\tilde{x};\theta)=c^{*}(\theta)exp\{\sum_{j=1}^{k}Q_{j}(\theta)T_{j}^{*}(\tilde{x})\}h^{*}(\tilde{x})$$
若 $$Q=(Q_{1}(\theta),\cdots,Q_{k}(\theta))$$ 的值域有非空的内部,则 $$T=(T_{1}^{*},...,T_{k}^{*})$$ 为**充分完备统计量**。

## §3.5 Cramér-Rao不等式——无偏估计方差的下界

### Fisher 信息量

设随机变量(或向量)X来自分布族 $$\mathcal{F}=\{p(x;\theta):\theta\in\Theta\}$$，则 $$\theta$$ 的 **Fisher Information** 定义为：
$$I(\theta):=Var_{\theta}\left\{\frac{\partial \log p(X;\theta)}{\partial\theta}\right\}=E_{\theta}\left[\frac{\partial \log p(X;\theta)}{\partial\theta}\right]^{2}$$
在正则条件下，Fisher信息量也可以表示为：
$$I(\theta)=-E_{\theta}\left[\frac{\partial^{2}\log p(X;\theta)}{\partial\theta^{2}}\right]$$

#### 样本的Fisher信息量
如果总体X的Fisher信息量存在，记为 $$I(\theta)$$。$$\tilde{X}=(X_{1},X_{2},...,X_{n})$$ 是来自X的简单随机样本，其联合pdf为 $$p(\tilde{x};\theta)$$，那么 $$\tilde{X}$$ 的Fisher Information为：
$$I_{n}(\theta)=Var_{\theta}\left\{\frac{\partial \log p(\tilde{X};\theta)}{\partial\theta}\right\}=nI(\theta)$$

### Cramér-Rao 不等式

**定理 3.5.1 (Cramér-Rao 不等式)**: 设总体X来自分布族 $$\{p(x;\theta),\theta\in\Theta\}$$，样本 $$\tilde{X}$$ 的联合密度函数为 $$p(\tilde{x};\theta)$$。设 $$\hat{g}=\hat{g}(\tilde{X})$$ 是 $$g(\theta)$$ 的无偏估计 (其中 $$g(\cdot)$$ 可微)，且 $$Var_{\theta}\{\hat{g}\}<\infty$$。在正则条件下，则：
$$Var_{\theta}\{\hat{g}\}\ge\frac{(g^{\prime}(\theta))^{2}}{E_{\theta}[\frac{\partial}{\partial\theta}\log p(\tilde{X};\theta)]^{2}}$$
上式右端称为参数$$g(\theta)$$的无偏估计方差的**Cramér-Rao下界** (简称C-R下界)。

**推论**: 设 $$\tilde{X}=(X_{1},X_{2},...,X_{n})$$ 为取自总体 $$X\sim p(x;\theta)$$ 的i.i.d.样本，在正则条件下，对 $$g(\theta)$$ 的任一无偏估计 $$\hat{g}(\tilde{X})$$，有：
$$Var_{\theta}\{\hat{g}\}\ge\frac{(g^{\prime}(\theta))^{2}}{nE_{\theta}[\frac{\partial}{\partial\theta}\log p(X;\theta)]^{2}}$$

### 有效估计和估计的效率

**结论**: 如果C-R定理中的条件满足，则达到C-R下界的无偏估计一定是UMVUE。

**定义 3.5.2 (效率)**: 设 $$\hat{g}(\tilde{X})$$ 是 $$g(\theta)$$ 的无偏估计，称
$$e_{n}=\frac{[g^{\prime}(\theta)]^{2}/[nI(\theta)]}{Var_{\theta}\{\hat{g}(\tilde{X})\}}$$
为无偏估计 $$\hat{g}(\tilde{X})$$ 的**效率** (efficiency)。
* 如果 $$e_n=1$$，则称 $$\hat{g}(\tilde{X})$$ 为 $$g(\theta)$$ 的**有效估计**。
* 如果 $$e_{n}\rightarrow1$$，当 $$n\rightarrow\infty$$，则称 $$\hat{g}(\tilde{X})$$ 为 $$g(\theta)$$ 的**渐近有效估计**。

如果 $$\hat{g}(\tilde{X})$$ 是 $$g(\theta)$$ 的估计，满足
$$\sqrt{n}\{\hat{g}(\tilde{X})-g(\theta)\}\rightarrow N(0,V(\theta)), \quad \text{当} \ n\rightarrow\infty$$
也把
$$e=\frac{[g^{\prime}(\theta)]^{2}/I(\theta)}{V(\theta)}$$
称为**渐近效率**。