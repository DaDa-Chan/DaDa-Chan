---
title: "Appendix: Common Sufficient & Complete Statistics"
date: 2025-07-19
category: "Mathematical Statistics"
excerpt: "A reference table of complete sufficient statistics and MLEs for common exponential family distributions."
---

完备统计量通常存在于**指数型分布族 (Exponential Family)** 中。指数型分布族是一类具有良好性质的分布族，在统计推断中扮演着核心角色。

对于属于指数型分布族的分布，通常存在一个**完备充分统计量**。

以下是一些常见的指数型分布族及其对应的完备统计量：

1.  **正态分布 (Normal Distribution) $$N(\mu, \sigma^2)$$**
    * 如果 $$\sigma^2$$ 已知，$$\mu$$ 未知，则 $$\sum_{i=1}^n X_i$$ 或样本均值 $$\bar{X}$$ 是完备充分统计量。
    * 如果 $$\mu$$ 已知，$$\sigma^2$$ 未知，则 $$\sum_{i=1}^n (X_i - \mu)^2$$ 是完备充分统计量。
    * 如果 $$\mu$$ 和 $$\sigma^2$$ 都未知，则 $$(\sum_{i=1}^n X_i, \sum_{i=1}^n X_i^2)$$ 是完备充分统计量。等价地，$$(\bar{X}, S^2)$$ 也是完备充分统计量（其中 $$S^2$$ 是样本方差）。

2.  **伯努利分布 (Bernoulli Distribution) $$B(1, p)$$**
    * 对于参数 $$p$$ (成功概率)，样本和 $$\sum_{i=1}^n X_i$$ 是完备充分统计量。

3.  **二项分布 (Binomial Distribution) $$B(n, p)$$**
    * 对于参数 $$p$$，样本和 $$\sum_{i=1}^k X_i$$ (假设有 $$k$$ 个独立的二项试验) 是完备充分统计量。

4.  **泊松分布 (Poisson Distribution) $$P(\lambda)$$**
    * 对于参数 $$\lambda$$，样本和 $$\sum_{i=1}^n X_i$$ 是完备充分统计量。

5.  **指数分布 (Exponential Distribution) $$Exp(\lambda)$$**
    * 对于参数 $$\lambda$$，样本和 $$\sum_{i=1}^n X_i$$ 是完备充分统计量。

6.  **伽马分布 (Gamma Distribution) $$Gamma(\alpha, \beta)$$**
    * 如果 $$\alpha$$ 已知，$$\beta$$ 未知，则 $$\sum_{i=1}^n X_i$$ 是完备充分统计量。
    * 如果 $$\beta$$ 已知，$$\alpha$$ 未知，则 $$\sum_{i=1}^n \ln(X_i)$$ 是完备充分统计量。
    * 如果 $$\alpha$$ 和 $$\beta$$ 都未知，则 $$(\sum_{i=1}^n X_i, \sum_{i=1}^n \ln(X_i))$$ 是完备充分统计量。

7.  **贝塔分布 (Beta Distribution) $$Beta(\alpha, \beta)$$**
    * 如果 $$\alpha$$ 已知，$$\beta$$ 未知，则 $$\sum_{i=1}^n \ln(1-X_i)$$ 是完备充分统计量。
    * 如果 $$\beta$$ 已知，$$\alpha$$ 未知，则 $$\sum_{i=1}^n \ln(X_i)$$ 是完备充分统计量。
    * 如果 $$\alpha$$ 和 $$\beta$$ 都未知，则 $$(\sum_{i=1}^n \ln(X_i), \sum_{i=1}^n \ln(1-X_i))$$ 是完备充分统计量。

8.  **均匀分布 (Uniform Distribution) $$U(a, b)$$**
    * 对于 $$U(0, \theta)$$，最大值统计量 $$\max(X_1, \dots, X_n)$$ 是完备充分统计量。
    * 对于 $$U(\theta, \theta+1)$$ (长度已知为1)，$$(\min(X_i), \max(X_i))$$ 是完备充分统计量，但单一统计量通常不完备。
    * 对于 $$U(a, b)$$，最小值和最大值统计量 $$(\min(X_1, \dots, X_n), \max(X_1, \dots, X_n))$$ 是完备充分统计量。
---

### 离散分布 (Discrete Distributions)

| 分布 (Distribution) | 概率质量函数 (PMF) | 参数 (Parameter) | 最大似然估计 (MLE) |
| :--- | :--- | :--- | :--- |
| **伯努利分布 (Bernoulli)** | $$P(X=x) = p^x(1-p)^{1-x}$$ for $$x \in \{0, 1\}$$ | $$p$$ (成功概率) | $$\hat{p} = \frac{1}{n}\sum_{i=1}^{n}X_i = \bar{X}$$ |
| **二项分布 (Binomial)** | $$P(X=k) = \binom{n}{k} p^k(1-p)^{n-k}$$ | $$p$$ (成功概率) | $$\hat{p} = \frac{\sum_{i=1}^{N}X_i}{nN} = \frac{\bar{X}}{n}$$ (其中N是试验次数) |
| **泊松分布 (Poisson)** | $$P(X=k) = \frac{\lambda^k e^{-\lambda}}{k!}$$ | $$\lambda$$ (率参数) | $$\hat{\lambda} = \frac{1}{n}\sum_{i=1}^{n}X_i = \bar{X}$$ |
| **几何分布 (Geometric)** | $$P(X=k) = (1-p)^{k-1}p$$ | $$p$$ (成功概率) | $$\hat{p} = \frac{n}{\sum_{i=1}^{n}X_i} = \frac{1}{\bar{X}}$$ |
| **负二项分布 (Negative Binomial)** | $$P(X=k) = \binom{k-1}{r-1}p^r(1-p)^{k-r}$$ | $$p$$ (成功概率) | $$\hat{p} = \frac{r}{\bar{X}}$$ |

---

### 连续分布 (Continuous Distributions)

| 分布 (Distribution)      | 概率密度函数 (PDF)                                                           | 参数 (Parameter)                   | 最大似然估计 (MLE)                                                                                                          |
| :--------------------- | :--------------------------------------------------------------------- | :------------------------------- | :-------------------------------------------------------------------------------------------------------------------- |
| **正态分布 (Normal)**      | $$f(x) = \frac{1}{\sqrt{2\pi\sigma^2}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}$$ | $$\mu$$ (均值)<br>$$\sigma^2$$ (方差)    | $$\hat{\mu} = \frac{1}{n}\sum_{i=1}^{n}X_i = \bar{X}$$<br>$$\hat{\sigma}^2 = \frac{1}{n}\sum_{i=1}^{n}(X_i - \bar{X})^2$$ |
| **均匀分布 (Uniform)**     | $$f(x) = \frac{1}{b-a}$$ for $$x \in [a, b]$$                              | $$a, b$$ (区间端点)                    | $$\hat{a} = \min(X_1, ..., X_n)$$<br>$$\hat{b} = \max(X_1, ..., X_n)$$                                                    |
| **指数分布 (Exponential)** | $$f(x) = \lambda e^{-\lambda x}$$ for $$x \ge 0$$                          | $$\lambda$$ (率参数)                  | $$\hat{\lambda} = \frac{n}{\sum_{i=1}^{n}X_i} = \frac{1}{\bar{X}}$$                                                     |
| **伽马分布 (Gamma)**       | $$f(x) = \frac{\beta^\alpha}{\Gamma(\alpha)}x^{\alpha-1}e^{-\beta x}$$   | $$\alpha$$ (形状参数)<br>$$\beta$$ (率参数) | 无显式解 (No closed-form solution for $$\alpha$$);<br>$$\hat{\beta} = \frac{\hat{\alpha}}{\bar{X}}$$                          |
| **拉普拉斯分布 (Laplace)**   | $$f(x) = \frac{1}{2b}e^{-\frac{x-\mu}{b}}$$                              | $$\mu$$ (位置参数)<br>$$b$$ (尺度参数)       | $$\hat{\mu} = \text{median}(X_1, ..., X_n)$$<br>$$\hat{b} = \frac{1}{n}\sum_{i=1}^{n}X_i - \hat{\mu}$$                    |
| **柯西分布 (Cauchy)**      | $$f(x) = \frac{1}{\pi\gamma[1+(\frac{x-x_0}{\gamma})^2]}$$               | $$x_0$$ (位置参数)<br>$$\gamma$$ (尺度参数)  | 无显式解 (No closed-form solution)                                                                                        |

**注意:**

* $$X_1, ..., X_n$$ 是来自相应分布的独立同分布 (i.i.d.) 样本。
* $$\bar{X} = \frac{1}{n}\sum_{i=1}^{n}X_i$$ 是样本均值。
* 对于某些分布（如伽马分布和柯西分布），其参数的最大似然估计没有简单的封闭解析解 (closed-form solution)，需要通过数值优化方法（如牛顿法）来求解。
* 正态分布方差的 MLE $$\hat{\sigma}^2 = \frac{1}{n}\sum_{i=1}^{n}(X_i - \bar{X})^2$$ 是有偏估计 (biased estimator)。其无偏估计 (unbiased estimator) 是 $$S^2 = \frac{1}{n-1}\sum_{i=1}^{n}(X_i - \bar{X})^2$$。