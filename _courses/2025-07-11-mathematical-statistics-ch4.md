---
title: "Ch4 Interval Estimation"
date: 2025-07-11
category: "Mathematical Statistics"
excerpt: "Confidence intervals, pivotal quantity method, normal population parameters, and general construction methods."
---

## §4.1 区间估计的基本概念

### 一、区间估计
**定义 4.1.1**: $$\theta\in\Theta\subset\mathbb{R}$$ 是待估参数，设 $$\hat{\theta}_{L}=\hat{\theta}_{L}(\tilde{X})$$ 和 $$\hat{\theta}_{U}=\hat{\theta}_{U}(\tilde{X})$$ 是参数空间上取值的样本$$\tilde{X}$$的两个统计量，且满足 $$\hat{\theta}_{L}(\tilde{X})\le\hat{\theta_{U}}(\tilde{X})$$。则称随机区间 $$[\hat{\theta}_{L}(\tilde{X}), \hat{\theta}_{U}(\tilde{X})]$$ 为参数$$\theta$$的**区间估计(量)** (interval estimation)。

### 二、区间估计的评价标准
#### 1. 置信度
**定义**: 参数$$\theta$$的区间估计 $$[\hat{\theta}_{L}(\tilde{X}),\hat{\theta_{U}}(\tilde{X})]$$ 的**涵盖概率** (coverage probability) (也称置信度)，是指它涵盖待估参数$$\theta$$的概率：
$$P_{\theta}\{\theta\in[\hat{\theta}_{L}(\tilde{X}),\hat{\theta_{U}}(\tilde{X})]\}$$
**定义 4.1.2**: 参数$$\theta$$的区间估计 $$[\hat{\theta}_{L}(\tilde{X}),\hat{\theta_{U}}(\tilde{X})]$$ 的置信度在参数空间上的下确界：
$$inf_{\theta\in\Theta}P_{\theta}\{\theta\in[\hat{\theta}_{L}(\tilde{X}),\hat{\theta_{U}}(\tilde{X})]\}$$
称为该区间估计的**置信系数** (confidence coefficient)。

#### 2. 精确度
最常用的标准是随机区间 $$[\hat{\theta}_{L}(\tilde{X}), \hat{\theta}_{U}(\tilde{X})]$$ 的**平均长度**:
$$E\{\hat{\theta}_{U}(\tilde{X})-\hat{\theta_{L}}(\tilde{X})\}$$
希望随机区间的平均长度越小越好。

**Neyman原则**: 在保证置信系数达到指定要求的前提下，尽可能提高精确度。

### 三、置信区间
**定义 4.1.3**: 对给定的 $$\alpha (0<\alpha<1)$$，参数$$\theta$$的置信系数不低于$$1-\alpha$$的区间估计 $$[\hat{\theta}_{L}(\tilde{X}),\hat{\theta_{U}}(\tilde{X})]$$ 称为$$\theta$$的置信水平为$$1-\alpha$$的(双侧)**置信区间** (confidence interval)，即
$$P_{\theta}\{\theta\in[\hat{\theta_{L}}(\tilde{X}),\hat{\theta_{U}}(\tilde{X})]\}\ge1-\alpha, \quad \forall\theta\in\Theta$$
此时 $$\hat{\theta}_{L}(\tilde{X})$$ 和 $$\hat{\theta}_{U}(\tilde{X})$$ 分别称为$$\theta$$的置信水平为$$1-\alpha$$的(双侧)**置信下限**和**置信上限**。

### 四、同等置信区间
**定义**: 设 $$[\hat{\theta}_{L}(\tilde{X}), \hat{\theta}_{U}(\tilde{X})]$$ 是参数$$\theta$$的一个区间估计，假定对给定的$$\alpha (0<\alpha<1)$$，有
$$inf_{\theta\in\Theta}P_{\theta}\{\theta\in[\hat{\theta}_{L}(\tilde{X}),\hat{\theta_{U}}(\tilde{X})]\}=1-\alpha$$
则称 $$[\hat{\theta}_{L}(\tilde{X}), \hat{\theta}_{U}(\tilde{X})]$$ 是$$\theta$$的置信水平为$$1-\alpha$$的**同等(双侧)置信区间**。

### 五、置信限 (单侧区间估计)
**定义 4.1.4**: 对给定的 $$\alpha (0<\alpha<1)$$，
* 若参数的区间估计 $$(-\infty, \hat{\theta_{U}}(\tilde{X})]$$ 的置信系数不低于$$1-\alpha$$，即
    $$P_{\theta}\{\theta\le\hat{\theta}_{U}\}\ge1-\alpha, \quad \forall\theta\in\Theta$$
    则称 $$\hat{\theta}_{U}(\tilde{X})$$ 为$$\theta$$的置信水平为$$1-\alpha$$的**单侧置信上限**。
* 若参数的区间估计 $$[\hat{\theta}_{L}(\tilde{X}),+\infty)$$ 的置信系数不低于$$1-\alpha$$，即
    $$P_{\theta}\{\theta\ge\hat{\theta}_{L}\}\ge1-\alpha, \quad \forall\theta\in\Theta$$
    则称$$\hat{\theta}_{L}(\tilde{X})$$为$$\theta$$的置信水平为$$1-\alpha$$的**单侧置信下限**。

**Lemma**: 设 $$\hat{\theta}_{L}(\tilde{X})$$ 和 $$\hat{\theta}_{U}(\tilde{X})$$ 分别是参数$$\theta$$的置信水平为 $$1-\alpha_{1}$$ 和 $$1-\alpha_{2}$$ 的单侧置信下、上限，且对任何样本 $$\tilde{X}$$，都有 $$\hat{\theta}_{L}(\tilde{X})\le\hat{\theta}_{U}(\tilde{X})$$。那么 $$[\hat{\theta_{L}}(\tilde{X}),\hat{\theta_{U}}(\tilde{X})]$$ 是$$\theta$$的置信水平为$$1-(\alpha_{1}+\alpha_{2})$$ 的双侧置信区间 ($$0<\alpha_{1}+\alpha_{2}<1$$)。

### 六、置信域
**定义**: 设待估参数 $$\theta=(\theta_{1},\cdots,\theta_{k})\in\Theta\subset\mathbb{R}^{k}$$。假设样本 $$\tilde{X}=(X_{1},\cdots,X_{n})$$ 的集合函数 $$S(\tilde{X})$$ 满足：
* 对任一样本观察值，S(x)是$$\Theta$$的一个子集。
* 对给定的$$\alpha (0<\alpha<1)$$，
    $$P_{\theta}\{\theta\in S(\tilde{X})\}\ge1-\alpha, \quad \forall\theta\in\Theta$$
    则称S(X)是$$\theta$$的置信水平为$$1-\alpha$$的**置信域** (或置信集)。

---
## §4.2 枢轴(变)量法 (Pivotal Quantity)

### 一、枢轴(变)量法构造同等置信区间的步骤:
1.  构造一个样本$$\tilde{X}$$和待估参数$$\theta$$的函数 $$G(\tilde{X},\theta)$$，要求G的分布不依赖于任何未知参数。一般称具有这种性质的函数 $$G(\tilde{X},\theta)$$ 为**枢轴(变)量**。
2.  对给定的$$\alpha (0<\alpha<1)$$，适当地选取两个常数c和d,使得
    $$P_{\theta}\{c\le G(\tilde{X},\theta)\le d\}=1-\alpha$$
3.  若能将不等式 $$c\le G(\tilde{X},\theta)\le d$$ 等价地转化为 $$\hat{\theta}_{L}(\tilde{X})\le\theta\le\hat{\theta}_{U}(\tilde{X})$$，那么 $$[\hat{\theta}_{L}(\tilde{X}),\hat{\theta}_{U}(\tilde{X})]$$ 是$$\theta$$的置信水平为$$1-\alpha$$的同等置信区间。

### 二、正态总体参数的置信区间
设 $$X_{1}, \dots, X_n$$ 是取自正态总体 $$N(\mu, \sigma^2)$$ 的样本。

#### 1. 正态总体均值 $$\mu$$ 的置信区间
* **$$\sigma^2$$已知时**: 枢轴量为 
  $$G(\tilde{X},\mu)=\frac{\overline{X}-\mu}{\sigma/\sqrt{n}}\sim N(0,1)$$
    $$\mu$$的置信水平为$$1-\alpha$$的最优置信区间为：
    $$[\overline{X}-u_{\alpha/2}\frac{\sigma}{\sqrt{n}},\overline{X}+u_{\alpha/2}\frac{\sigma}{\sqrt{n}}]$$
* **$$\sigma^2$$未知时**: 枢轴量为 $$T=\frac{\overline{X}-\mu}{S/\sqrt{n}}\sim t(n-1)$$
    $$\mu$$的置信水平为$$1-\alpha$$的最优置信区间为：
    $$[\overline{X}-t_{\alpha/2}(n-1)\cdot S/\sqrt{n},\overline{X}+t_{\alpha/2}(n-1)\cdot S/\sqrt{n}]$$

#### 2. 正态总体方差 $$\sigma^2$$ 的置信区间
* **$$\mu$$已知时**: 枢轴量为 $$\frac{\sum_{i=1}^{n}(X_{i}-\mu)^{2}}{\sigma^{2}}\sim\chi^{2}(n)$$
    $$\sigma^2$$的置信水平为$$1-\alpha$$的同等置信区间为：
    $$[\sum_{i=1}^{n}(X_{i}-\mu)^{2}/\chi_{\alpha/2}^{2}(n),\sum_{i=1}^{n}(X_{i}-\mu)^{2}/\chi_{1-\alpha/2}^{2}(n)]$$
* **$$\mu$$未知时**: 枢轴量为 $$\frac{(n-1)S^{2}}{\sigma^{2}}\sim\chi^{2}(n-1)$$
    $$\sigma^2$$的置信水平为$$1-\alpha$$的同等置信区间为：
    $$[(n-1)S^{2}/\chi_{\alpha/2}^{2}(n-1),(n-1)S^{2}/\chi_{1-\alpha/2}^{2}(n-1)]$$

#### 3. 两个正态总体
设 $$X \sim N(\mu_{1},\sigma_{1}^{2})$$ 与 $$Y\sim N(\mu_{2},\sigma_{2}^{2})$$ 是两个相互独立的正态总体，从两个总体中分别抽取样本 $$\tilde{X}=(X_{1},\dots,X_{m})$$ 和 $$\tilde{Y}=(Y_{1},\dots,Y_{n})$$，且两样本相互独立。

##### (一) 均值差 $$\delta=\mu_{1}-\mu_{2}$$ 的置信区间

**(A) $$\sigma_{1}^{2}$$ 和 $$\sigma_{2}^{2}$$ 都已知的情形**
* **枢轴量**:
    $$G=\frac{(\overline{X}-\overline{Y})-(\mu_1-\mu_2)}{\sqrt{\sigma_{1}^{2}/m+\sigma_{2}^{2}/n}}\sim N(0,1)$$
* **置信区间**: $$\delta$$ 的置信水平为 $$1-\alpha$$ 的最优同等置信区间为
    $$[(\overline{X}-\overline{Y})-u_{\alpha/2}\sqrt{\sigma_{1}^{2}/m+\sigma_{2}^{2}/n}, (\overline{X}-\overline{Y})+u_{\alpha/2}\sqrt{\sigma_{1}^{2}/m+\sigma_{2}^{2}/n}]$$

**(B) $$\sigma_{1}^{2}$$ 和 $$\sigma_{2}^{2}$$ 都未知，但 $$\sigma_{1}^{2}=\sigma_{2}^{2}=\sigma^{2}$$**
* **枢轴量**:
    $$t=\frac{(\overline{X}-\overline{Y})-(\mu_1-\mu_2)}{S_{W}\sqrt{1/m+1/n}}\sim t(m+n-2)$$
    其中 $$S_{W}^{2}=\frac{(m-1)S_{X}^{2}+(n-1)S_{Y}^{2}}{m+n-2}$$。
* **置信区间**: $$\delta$$ 的置信水平为 $$1-\alpha$$ 的同等置信区间为
    $$[(\overline{X}-\overline{Y})-t_{\alpha/2}(m+n-2)S_{W}\sqrt{1/m+1/n}, (\overline{X}-\overline{Y})+t_{\alpha/2}(m+n-2)S_{W}\sqrt{1/m+1/n}]$$

**(C) $$\sigma_{1}^{2}$$ 和 $$\sigma_{2}^{2}$$ 都未知，但方差比 $$\sigma_{2}^{2}/\sigma_{1}^{2}=\lambda$$ 已知**
* **枢轴量**:
    $$t=\frac{(\overline{X}-\overline{Y})-(\mu_1-\mu_2)}{\sqrt{\frac{(m-1)S_X^2+(n-1)S_Y^2/\lambda}{m+n-2}}\sqrt{1/m+\lambda/n}} \sim t(m+n-2)$$
* **置信区间**: $$\delta$$ 的置信水平为 $$1-\alpha$$ 的同等置信区间的上下限为
    $$(\overline{X}-\overline{Y}) \pm t_{\alpha/2}(m+n-2)\sqrt{\frac{(m-1)S_X^2+(n-1)S_Y^2/\lambda}{m+n-2}(1/m+\lambda/n)}$$

**(D) $$\sigma_{1}^{2}$$ 和 $$\sigma_{2}^{2}$$ 都未知，但 $$m=n$$ (成对数据)**
* **枢轴量**: 令 $$Z_i = X_i - Y_i$$，则 $$t=\frac{\overline{Z}-(\mu_1-\mu_2)}{S_Z/\sqrt{n}} \sim t(n-1)$$。
* **置信区间**: $$\delta$$ 的置信水平为 $$1-\alpha$$ 的同等置信区间为
    $$[\overline{Z}-t_{\alpha/2}(n-1)S_Z/\sqrt{n}, \overline{Z}+t_{\alpha/2}(n-1)S_Z/\sqrt{n}]$$

**(E) 一般情形 (Behrens-Fisher 问题)**
* **大样本**: 当m和n都充分大时，近似枢轴量为
    $$U=\frac{(\overline{X}-\overline{Y})-(\mu_1-\mu_2)}{\sqrt{S_{X}^{2}/m+S_{Y}^{2}/n}} \rightarrow N(0,1)$$
    $$\delta$$ 的置信系数近似为 $$1-\alpha$$ 的置信区间为
    $$[(\overline{X}-\overline{Y})-u_{\alpha/2}\sqrt{S_{X}^{2}/m+S_{Y}^{2}/n}, (\overline{X}-\overline{Y})+u_{\alpha/2}\sqrt{S_{X}^{2}/m+S_{Y}^{2}/n}]$$
* **小样本**: 近似枢轴量为
    $$T=\frac{(\overline{X}-\overline{Y})-(\mu_1-\mu_2)}{\sqrt{S_{X}^{2}/m+S_{Y}^{2}/n}}$$
    其近似服从 $$t(l)$$ 分布，其中自由度 $$l$$ 为
    $$l=\frac{({S_{X}^{2}/m+S_{Y}^{2}/n})^{2}}{\frac{({S_{X}^{2}/m})^2}{m-1}+\frac{({S_{Y}^{2}/n})^2}{n-1}}$$
    $$\delta$$ 的置信系数近似为 $$1-\alpha$$ 的置信区间为
    $$[(\overline{X}-\overline{Y})-t_{\alpha/2}(l)\sqrt{S_{X}^{2}/m+S_{Y}^{2}/n}, (\overline{X}-\overline{Y})+t_{\alpha/2}(l)\sqrt{S_{X}^{2}/m+S_{Y}^{2}/n}]$$

##### (二) 方差比 $$\sigma_{1}^{2}/\sigma_{2}^{2}$$ 的置信区间

**(A) $$\mu_{1}$$ 和 $$\mu_{2}$$ 都已知时**
* **枢轴量**:
    $$F=\frac{\sum_{i=1}^{m}(X_{i}-\mu_{1})^{2}/m}{\sum_{j=1}^{n}(Y_{j}-\mu_{2})^{2}/n} \cdot \frac{\sigma_2^2}{\sigma_1^2} \sim F(m,n)$$
* **置信区间**: $$\sigma_{1}^{2}/\sigma_{2}^{2}$$ 的置信水平为 $$1-\alpha$$ 的同等置信区间为
    $$\left[\frac{\sum(X_i-\mu_1)^2/m}{\sum(Y_j-\mu_2)^2/n} \frac{1}{F_{\alpha/2}(m,n)}, \frac{\sum(X_i-\mu_1)^2/m}{\sum(Y_j-\mu_2)^2/n} \frac{1}{F_{1-\alpha/2}(m,n)}\right]$$

**(B) $$\mu_{1}$$ 和 $$\mu_{2}$$ 都未知时**
* **枢轴量**:
    $$F=\frac{S_{X}^{2}/S_{Y}^{2}}{\sigma_{1}^{2}/\sigma_{2}^{2}}\sim F(m-1,n-1)$$
* **置信区间**: $$\sigma_{1}^{2}/\sigma_{2}^{2}$$ 的置信水平为 $$1-\alpha$$ 的同等置信区间为
    $$\left[\frac{S_{X}^{2}}{S_{Y}^{2}}\frac{1}{F_{\alpha/2}(m-1,n-1)}, \frac{S_{X}^{2}}{S_{Y}^{2}}\frac{1}{F_{1-\alpha/2}(m-1,n-1)}\right]$$
    利用 $$F_{1-\alpha/2}(m-1,n-1)=\frac{1}{F_{\alpha/2}(n-1,m-1)}$$ 可简化。


### 三、非正态总体的区间估计 (大样本方法)
当样本容量充分大时，可用渐近分布来构造近似的置信区间。
例如，由中心极限定理，对未知均值$$\mu$$，$$\frac{\overline{X}-\mu}{\sigma/\sqrt{n}}$$ 近似服从 $$N(0,1)$$。若$$\sigma$$未知，可用样本标准差$$S$$代替，构造近似枢轴量。


---

## §4.3 寻找枢轴量的一般方法——将统计量枢轴量化

如果T的连续的分布函数为 $$F(t;\theta)$$, 则 $$F(T;\theta)\sim U(0,1)$$。
$$F(T;\theta)$$ 可作为枢轴量，并且
$$P(\alpha/2\le F(T;\theta)\le1-\alpha/2)=1-\alpha/2-\alpha/2=1-\alpha$$
这时，若要求置信水平为 $$1-\alpha$$ 的置信区间，只要等价变形
$$\alpha/2\le F(T;\theta)\le1-\alpha/2$$
即可。

### 寻找置信区间的一般方法

#### 引理

**引理**: 设 $$F(x)$$ 是随机变量X的分布函数，则
$$P\{F(X)\le y\}\le y\le P\{F(X-0)<y\}, \quad \forall0<y<1$$
特别地,若 $$F(x)$$ 是连续函数,则 $$F(X)\sim U(0,1)$$。

#### 定理

设 $$\tilde{X}=(X_{1},X_{2},\cdots,X_{n})$$ 是来自总体 $$F_{\theta}(x)$$ 的一个样本，其中 $$\theta\in\Theta\subset R$$。$$T=T(\tilde{X})$$ 是一个统计量，其分布函数记为 $$G(t;\theta)$$。

**定理 1 (连续型分布函数情形)**: 设统计量 $$T=T(\tilde{X})$$ 的分布函数 $$G(t;\theta)$$ 是t的连续函数，也是$$\theta$$的连续函数。对 $$\alpha_{1}+\alpha_{2}=\alpha (0<\alpha<1)$$ 定义函数 $$\theta_{1}(t)$$ 和 $$\theta_{2}(t)$$ 如下:
$$1-G(t-0;\theta_{1}(t))=P_{\theta_{1}(t)}\{T\ge t\}=\alpha_{1}$$
$$G(t;\theta_{2}(t))=P_{\theta_{2}(t)}\{T\le t\}=\alpha_{2}$$
1. 若 $$G(t;\theta)$$ 是$$\theta$$的严格减函数，则 $$[\theta_{1}(T),\theta_{2}(T)]$$ 是$$\theta$$的置信水平为 $$1-\alpha$$ 的同等置信区间。
2. 若 $$G(t;\theta)$$ 是$$\theta$$的严格增函数，则 $$[\theta_{2}(T),\theta_{1}(T)]$$ 是$$\theta$$的置信水平为 $$1-\alpha$$ 的同等置信区间。

**定理 2 (分布函数关于参数连续的场合)**: 设统计量 $$T=T(\tilde{X})$$ 的分布函数 $$G(t;\theta)$$ 是$$\theta$$的连续函数。对 $$\alpha_{1}+\alpha_{2}=\alpha (0<\alpha<1)$$ 定义函数 $$\theta_1(t)$$ 和 $$\theta_2(t)$$ 如下:
$$1-G(t-0;\theta_{1}(t))=P_{\theta_{1}(t)}\{T\ge t\}=\alpha_{1}$$
$$G(t;\theta_{2}(t))=P_{\theta_{2}(t)}\{T\le t\}=\alpha_{2}$$
* 若 $$G(t;\theta)$$ 是$$\theta$$的严格减函数，则 $$[\theta_{1}(T),\theta_{2}(T)]$$ 是$$\theta$$的置信水平为$$1-\alpha$$的置信区间。
* 若 $$G(t;\theta)$$ 是$$\theta$$的严格增函数，则 $$[\theta_{2}(T),\theta_{1}(T)]$$ 是$$\theta$$的置信水平为$$1-\alpha$$的置信区间。
(注: 不一定是同等置信区间。)

#### 一般情形的结果

设 $$\tilde{X}=(X_{1},X_{2},\cdots,X_{n})$$ 是来自总体 $$F_{\theta}(x)$$ 的一组样本，其中 $$\theta\in\Theta\subset R$$。$$T=T(\tilde{X})$$ 是一个统计量，其分布函数 $$G(t;\theta)$$ 是参数的严格单调函数，这时置信区间可按下述方法构造:
* 取 $$\theta_{1}(t)$$ s.t. $$P_{\theta_{1}(t)}(T\ge t)\le\alpha_{1}$$ 但概率尽可能接近 $$\alpha_{1}$$。
* 取 $$\theta_{2}(t)$$ s.t. $$P_{\theta_{2}(t)}(T\le t)\le\alpha_{2}$$ 但概率尽可能接近 $$\alpha_{2}$$。

则 $$[\theta_{1}(T),\theta_{2}(T)]$$ 或者 $$[\theta_{2}(T), \theta_{1}(T)]$$ 是$$\theta$$的置信水平为 $$1-\alpha_{1}-\alpha_{2}$$ 的置信区间。