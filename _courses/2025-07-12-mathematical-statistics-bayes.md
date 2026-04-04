---
title: "Bayesian Inference"
date: 2025-07-12
category: "Mathematical Statistics"
excerpt: "Prior and posterior distributions, Bayesian point estimation, credible intervals, Bayesian testing, and Inverse-Gamma distribution."
---

# Bayes 点估计法

### 一、统计推断中可用的三种信息
1.  **总体信息 (population information)**: 总体分布或总体所属分布族所提供的信息。
2.  **样本信息 (sample information)**: 样本提供给我们的信息。
3.  **先验信息 (prior information)**: 在抽样之前就有的有关统计推断问题的信息，常常是过去同类统计推断问题提供的信息。

* **经典统计学派**: 在统计推断中只利用总体信息和样本信息。
* **Bayes统计学派**: 建议在统计推断中在利用总体信息和样本信息的同时，还利用先验信息。

### 二、先验分布和后验分布
Bayes统计学派认为：未知参数$\theta$是个随机变量，它有一个分布 (pdf或pmf) $\pi(\theta)$，称为**先验分布**。

#### Bayes统计的三个基本假设:
* **假设I**: 总体X有一个概率分布 $p(x|\theta)$，这是给定$\theta$后的一个条件概率分布。
* **假设II**: 当给定$\theta$后，从总体中随机抽取一个样本 $\tilde{X}=(X_{1},\dots,X_{n})$，其联合条件概率分布为 $p(\tilde{x}|\theta)=\prod_{i=1}^{n}p(x_{i}|\theta)$。
* **假设III**: 未知参数$\theta$是一个随机变量，它有一个已知的先验分布 $\pi(\theta)$。

样本$\tilde{X}$和参数$\theta$的联合概率分布是
$$p(\tilde{x},\theta)=p(\tilde{x}|\theta)\pi(\theta)$$
给定样本 $\tilde{X}=\tilde{x}$，参数$\theta$的条件概率分布
$$\pi(\theta|\tilde{x})=\frac{p(\tilde{x},\theta)}{p(\tilde{x})} = \frac{p(\tilde{x}|\theta)\pi(\theta)}{\int p(\tilde{x}|\theta)\pi(\theta)d\theta}$$
称为**后验分布** (posterior distribution)。

**结论**: 用充分统计量 $T=T(\tilde{X})$ 代替样本$\tilde{X}$所得的后验分布是一样的，即 $\pi(\theta|\tilde{x}) = \pi(\theta|t)$，其中 $t=T(\tilde{x})$ 
从而
$$\pi(\theta|\tilde{x})=\pi(\theta|t) = \frac{p(t|\theta)\pi(\theta)}{\int p(t|\theta)\pi(\theta)d\theta}$$

### 三、Bayes点估计
Bayes方法对参数的统计推断是建立在后验分布 $\pi(\theta|\tilde{x})$ 的基础上的。
* **众数型Bayes估计**: 后验分布 $\pi(\theta|\tilde{x})$ 的众数。
* **中位数型Bayes估计**: 后验分布 $\pi(\theta|\tilde{x})$ 的中位数。
* **期望型Bayes估计**: 后验分布 $\pi(\theta|\tilde{x})$ 的期望。
    $$\hat{\theta}_{B}(\tilde{x})=E(\theta|\tilde{x})=\int\theta\cdot\pi(\theta|\tilde{x})d\theta$$

#### 先验分布的选取
1. “同等无知”
2. **“共轭分布法”**
**定义**: 设$\theta$是某分布的一个参数，$\pi(\theta)$是其先验分布。假如由抽样信息算得的后验分布 $\pi(\theta|\tilde{x})$与$\pi(\theta)$属于同一类分布族，则称$\pi(\theta)$是的**共轭先验分布**。

| 总体分布 | 参数 | 共轭先验分布 |
| :--- | :--- | :--- |
| 二项分布 | 成功概率 | beta分布 |
| 泊松分布 | 均值 | gamma分布 |
| 指数分布 | 均值 | 倒gamma分布 |
| 指数分布 | 均值倒数 | gamma分布 |
| 正态分布(方差已知) | 均值 | 正态分布 |
| 正态分布(均值已知) | 方差 | 倒gamma分布 |

### Bayes 预测
设 $X\sim p(x|\theta)$，在获取数据 $T=t$ 后，对具有概率密度函数 $g(z|\theta)$ 的随机变量Z的未来观察值作出预测。
在Bayes框架下，给定$T=t$，Z的**后验预测密度函数**为
$$\hat{p}_{Z|T}(z|t)=\int_{\Theta}g(z|\theta)\pi(\theta|t)d\theta$$
### 无信息先验
Bayes分析的一个重要特点就是在统计推断时要利用先验信息，但是实际中常常会出现这样的情况：没有先验信息或者只有极少的先验信息可利用，但仍然要用Bayes方法。

此时所需要的是一种**无信息先验** (noninformative prior)，即对参数空间中的任何一点没有偏爱的先验信息分布(“同等无知”)。

#### 均匀分布与广义先验分布
* 若 $\Theta=\{\theta_{1},...,\theta_{l}\}$ 为有限集，无信息先验给每个元素以相同的概率，即 $P(\theta=\theta_{i})=1/l, i=1,...,l$。
* 若 $\Theta=[a,b]$ 为有限区间，则取无信息先验为区间 $[a,b]$ 上的均匀分布 $U(a,b)$。
* 问题是若参数空间无界，无信息先验如何选取？例如，总体分布为 $N(\theta,1)$，此时 $\theta=(-\infty,\infty)$。若无信息先验，可取先验分布为 $\pi(\theta)\equiv1$，但它不是通常的密度函数，这就要引入广义先验分布的概念。

**定义 (广义先验)**: 设样本 $\tilde{X}\sim p(\tilde{x}|\theta)$，若 $\pi(\theta)$ 满足下列条件:
* $\pi(\theta)\ge0$ 且 $\int_{\Theta}\pi(\theta)d\theta=\infty$。
* 后验密度函数(或后验概率分布列)
    $$\pi(\theta|\tilde{x})=\frac{p(\tilde{x}|\theta)\pi(\theta)}{\int_{\Theta}p(\tilde{x}|\theta)\pi(\theta)d\theta}$$
    是正常的密度函数(或概率分布列)。
    则称$\pi(\theta)$为$\theta$的**广义先验密度** (improper prior density)。

#### 位置参数的无信息先验
设总体 $X\sim f(x-\theta)$, $-\infty<\theta<\infty$，$\theta$为**位置参数** (location parameter)。对X做平移变换 $Y=X+c$，同时对$\theta$也作同样的平移变换 $\eta=\theta+c$，显然$Y \sim f(y-\eta)$，$\eta$为位置参数。$(X,\theta)$ 与 $(Y,\eta)$ 的统计结构相同，因此主张它们有相同的无信息先验是合理的。这样$\theta$与$\eta=\theta+c$有相同的分布，从而$\theta$的先验分布不依赖于原点的选择，它在等长度区间内的先验概率应当一样，所以可取
$$\pi(\theta)\equiv1$$
它是一个广义先验密度。

#### 尺度参数的无信息先验
设总体 $X\sim\frac{g(x/\sigma)}{\sigma}$，$\sigma>0$ 为**尺度参数** (scale parameter)。对X作尺度变换 $Y=cX$ ($c>0$)，同时对$\sigma$也作同样的尺度变换 $\eta=c\sigma$，易知道 $Y\sim\frac{g(y/\eta)}{\eta}$，$\eta$为尺度参数。$(X,\sigma)$与$(Y,\eta)$的统计结构相同，因此主张它们有相同的无信息先验是合理的。这样$\sigma$与$\eta=c\sigma$有相同的分布，从而
$$\pi(\sigma)d\sigma = \pi(\eta)d\eta = \pi(c\sigma)cd\sigma$$
由此得 $\pi(c)=\pi(1)/c$，故可取$\sigma$的无信息先验为
$$\pi(\sigma)=\frac{1}{\sigma}, \quad \sigma>0$$
它是一个广义先验密度。

### 加权平均均方误差与Bayes 估计
考察使得均方误差 $R(\delta,\theta)=E_{\theta}\{(\theta-\delta(\tilde{X}))^{2}\}$ 关于$\theta$的加权平均最小的问题:
$$min_{\delta}\int R(\delta,\theta)\pi(\theta)d\theta$$
其中 $\pi(\theta)$ 为权函数，可看作先验分布。

**结论**: 对给定的权函数 $\pi(\theta)$，在加权均方误差(或称Bayes均方风险)意义下，**期望型Bayes估计**是最优估计。
$$E[\theta|\tilde{X}]=arg~min_{\delta(\tilde{X})}\int E_{\theta}\{\theta-\delta(\tilde{X})\}^{2}\pi(\theta)d\theta$$
若损失函数为绝对值损失 $L(\delta, \theta) = |\theta-\delta(\tilde{x})|$，则最优估计为后验分布的中位数。

#### Loss Function Optimality
* **损失函数 (Loss function)**: $L(\delta, \theta)$，表示决策$\delta$在参数为$\theta$时造成的损失。
* **风险函数 (Risk function)**: $R(\delta,\theta)=E_{\tilde{X}|\theta}[L(\delta(\tilde{X}),\theta)]=E[L(\delta(\tilde{X}),\theta)|\theta]$
* **Bayes风险 (Bayesian risk)**: $R_{\pi}(\delta)=\int_{\Theta}R(\delta,\theta)\pi(\theta)d\theta=E_{\tilde{X}}[R(\delta|\tilde{x})]$
* **后验风险 (Posterior risk)**: $R(\delta|\tilde{x})=E_{\theta|\tilde{x}}[L(\delta,\theta)]=\int_{\Theta}L(\delta,\theta)\pi(\theta|\tilde{x})d\theta$

**定理 (后验风险最小原则)**: 如果 $\delta_\pi$ 是在后验风险意义下的最优决策，即 $R(\delta_{\pi}|\tilde{x})=inf_{\delta}R(\delta|\tilde{x})$，则是Bayes风险意义下的最优决策。
$$arg~inf_{\delta}R(\delta)=arg~inf_{\delta}R(\delta|\tilde{x})$$

**一般损失下的Bayes估计**

- **加权平方损失下的Bayes估计**:  $L(\delta,\theta)=w(\theta)(\delta-\theta)^{2}$ ，$\theta$的Bayes估计值为
$$\hat{\theta_{B}}=\frac{E[\theta w(\theta)|\tilde{x}]}{E[w(\theta)|\tilde{x}]}$$
   其中 $w(\theta)$ 在参数空间上恒为正。

- **平方损失下的Bayes估计**：$L(\delta,\theta)=(\delta-\theta)^{2}$ , $\theta$的Bayes估计值为
$$\hat{\theta}_{B}(\tilde{x})=E(\theta|\tilde{x})=\int\theta\cdot\pi(\theta|\tilde{x})d\theta$$
- **绝对值损失下的Bayes估计**： $L(\delta,\theta)=|\delta-\theta|$ , $\theta$的Bayes估计值为 **后验中位数**

---
# 贝叶斯区间估计与假设检验

## Bayesian Intervals (贝叶斯区间估计)

在Bayes统计理论中，参数$\theta$是随机变量，我们可以讨论$\theta$落入某个区间的概率。
设 $\pi(\theta|\tilde{x})$ 是给定 $\tilde{X}=\tilde{x}$ 时，$\theta$的后验分布，称
$$P(\theta\in A|\tilde{x})=\int_{A}\pi(\theta|\tilde{x})d\theta$$
为A关于$\theta$的**可信概率** (credible probability)。

### 可信区间

**定义**: 设参数$\theta$的后验分布为 $\pi(\theta|\tilde{x})$。对给定的$1-\alpha$，若存在 $\hat{\theta}_{L}=\hat{\theta}_{L}(\tilde{x})$ 和 $\hat{\theta}_{U}=\hat{\theta}_{U}(\tilde{x})$，使得
$$P(\hat{\theta_{L}}\le\theta\le\hat{\theta_{U}}|\tilde{x})\ge1-\alpha$$
则称区间 $[\hat{\theta}_{L},\hat{\theta}_{U}]$ 为参数$\theta$的可信水平为$1-\alpha$的**可信区间** (Credible Interval)。

若
$$P(\theta\ge\hat{\theta}_{L}|\tilde{x})\ge1-\alpha$$
则称$\hat{\theta}_L$为参数$\theta$的可信水平为$1-\alpha$的**单侧可信下限**。同样可以定义单侧可信上限。

### 置信区间与贝叶斯可信区间的含义区别

* **置信区间**: 置信水平为90%的置信区间 $[\hat{\theta}_{L}(\tilde{x}),\hat{\theta}_{U}(\tilde{x})]$ 意味着，在一系列独立的重复试验中，平均有90%的这样构造出的区间会包含未知的参数真值。因此，这一个具体的区间客观上有90%的机会包含真值。
* **可信区间**: 可信水平为90%的可信区间 $[\hat{\theta}_{L}(\tilde{x}),\hat{\theta}_{U}(\tilde{x})]$ 意味着，试验者结合了先验信息和样本信息后，主观上有90%的把握说，未知参数$\theta$落入这一个具体的区间中。

---

## Bayesian Tests (贝叶斯检验)

考察假设检验问题：
$$H_{0}:\theta\in\Theta_{0}\leftrightarrow H_{1}:\theta\in\Theta_{1}$$
在求得$\theta$的后验分布 $\pi(\theta|\tilde{x})$ 后，计算$H_0$和$H_1$发生的概率：
$$\alpha_{0}=P(H_{0} \text{ is true}|\tilde{x})=P(\theta\in\Theta_{0}|\tilde{x})$$
$$\alpha_{1}=P(H_{1} \text{ is true}|\tilde{x})=P(\theta\in\Theta_{1}|\tilde{x})$$

### 检验法则

当 $\frac{\alpha_{0}}{\alpha_{1}}\ge1$ 时接受$H_0$，否则拒绝$H_0$。

若 $H_1 = \Theta_0^c$，则拒绝域为
$$\{\tilde{x}:P(\theta\in\Theta_{0}^{c}|\tilde{x})>\frac{1}{2}\}$$
如果希望保护原假设，拒绝域有时也可定义为
$$\{\tilde{x}:P(\theta\in\Theta_{0}^{c}|\tilde{x})>1-\alpha\}$$



# 倒Gamma分布 (Inverse-Gamma Distribution)

#### 概率密度函数 (PDF)
一个服从倒Gamma分布的随机变量 $X$，其概率密度函数由两个参数决定：形状参数 $\alpha > 0$ 和尺度参数 $\beta > 0$。其PDF为：
$$f(x; \alpha, \beta) = \frac{\beta^\alpha}{\Gamma(\alpha)} x^{-\alpha-1} \exp\left(-\frac{\beta}{x}\right), \quad x > 0$$
其中 $\Gamma(\alpha)$ 是标准的Gamma函数。

#### 数字特征 (Numerical Characteristics)
* **均值 (Mean)**:
    $$E[X] = \frac{\beta}{\alpha-1}, \quad \text{要求 } \alpha > 1$$
* **方差 (Variance)**:
    $$Var(X) = \frac{\beta^2}{(\alpha-1)^2(\alpha-2)}, \quad \text{要求 } \alpha > 2$$
* **众数 (Mode)**:
    $$\text{Mode}[X] = \frac{\beta}{\alpha+1}$$

#### 性质与关系
* **与Gamma分布的关系**:
    倒Gamma分布与Gamma分布关系密切。如果一个随机变量 $Y$ 服从 **Gamma分布** $Y \sim \Gamma(\alpha, \beta)$（其中 $\beta$ 是速率参数 rate parameter），那么它的倒数 $X = 1/Y$ 就服从**倒Gamma分布** $X \sim \text{Inverse-Gamma}(\alpha, \beta)$（其中 $\beta$ 是尺度参数 scale parameter）。

* **与卡方分布的关系**:
    **尺度倒卡方分布** (Scaled-inverse-chi-squared distribution) 是倒Gamma分布的一个特例。具体来说，自由度为 $\nu$、尺度参数为 $\sigma^2$ 的尺度倒卡方分布等价于 $\text{Inverse-Gamma}(\frac{\nu}{2}, \frac{\nu\sigma^2}{2})$ 分布。

* **作为共轭先验**:
    在贝叶斯推断中，倒Gamma分布扮演着至关重要的角色。当一个正态分布总体的均值已知而方差未知时，倒Gamma分布是其**方差**的**共轭先验分布**。这意味着，如果选择倒Gamma分布作为方差的先验分布，那么在观测到数据后，其后验分布仍然是倒Gamma分布，这极大地简化了计算。