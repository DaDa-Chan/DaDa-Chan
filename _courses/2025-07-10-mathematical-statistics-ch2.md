---
title: "Ch2 Sampling Distributions"
date: 2025-07-10
category: "Mathematical Statistics"
excerpt: "Normal sampling distributions, order statistics, common distribution families, exponential family, and sufficient statistics."
---

## §2.2 正态总体样本均值和样本方差的分布

### **已知结论**
总体 $X\sim F$，$X=(X_{1},X_{2},\cdots,X_{n})$ 为来自该总体的简单随机样本，$\overline{X}$和$S^2$为其样本均值与样本方差，若总体的二阶矩存在，并记 $EX=\mu, VarX=\sigma^{2}$，则有：
$$E\overline{X}=\mu, Var\overline{X}=\sigma^{2}/n$$
$$E(S^{2})=\sigma^{2}$$

### **Lemma**
设在两个随机变量 $\tilde{X}=(X_{1},X_{2},\cdots,X_{n})^{\prime}$ 与 $\tilde{Y}=(Y_{1},Y_{2},\cdots,Y_{n})^{\prime}$ 间有一个线性变换 $\tilde{Y}=A\tilde{X}$，其中 $A=(a_{ij})$ 为 $n \times n$ 方阵, 则：
$$E\tilde{Y}=A(E\tilde{X})$$
$$Var\tilde{Y}=A(Var\tilde{X})A^{\prime}$$

### **核心定理**
设 $X_{1},X_{2},\cdots,X_{n}$ 是取自正态总体 $N(\mu,\sigma^{2})$ 的一组简单随机样本。$\overline{X}$和$S^2$ 为其样本均值与样本方差, 则：
1.  $\overline{X}\sim N(\mu,\sigma^{2}/n)$
2.  $$\frac{(n-1)S^{2}}{\sigma^{2}}=\frac{nS_{n}^{2}}{\sigma^{2}}=\frac{1}{\sigma^{2}}\sum_{i=1}^{n}(X_{i}-\overline{X})^{2}\sim\chi^{2}(n-1)$$
3.  $\overline{X}$与$S^{2}$相互独立.



### 渐近分布 (Asymptotic distribution)
在大多数场合，精确抽样分布不易求出，或精确抽样分布过于复杂而难于应用，这时常常求统计量 $T=T(X_{1},X_{2},\cdots,X_{n})$ 当 $n \rightarrow \infty$ 时的极限分布，这种极限分布称为渐近分布。当n较大时，用渐近分布当作抽样分布的近似。

**例1**: 设 $X_{1},X_{2},\cdots,X_{n}$ 取自正态总体 $X\sim N(0,\sigma^{2})$ 的简单随机样本,则 $\sqrt{n}\overline{X}/\sigma\rightarrow N(0,1)$ 和 $S_{n}^{2}\rightarrow\sigma^{2}$。所以
$$T=\sqrt{n}\overline{X}/S_{n}=\frac{\sqrt{n}\overline{X}/\sigma}{S_{n}/\sigma}\rightarrow N(0,1), \quad 当 n\rightarrow\infty$$
即 $T\sim AN(0,1)$。

**例2**: 设 $X_{1},X_{2},\cdots,X_{n}$ 取自总体 $X\sim F(x)$ 的简单随机样本，且 $\mu_{2k}=EX^{2k}$ 存在。则
$$a_{n,k}\sim AN(\mu_{k},\frac{\mu_{2k}-(\mu_{k})^{2}}{n})$$

---

## §2.3 次序统计量及其分布

### 一、次序统计量

**定义**: 样本 $X_{1},X_{2},\cdots,X_{n}$ 按从小到大的顺序重排为 $X_{(1)}\le X_{(2)}\le...\le X_{(n)}$，则 $(X_{(1)},X_{(2)},...,X_{(n)})$ 称为次序统计量(order statistics)，特别地，$X_{(1)}$ 称为最小次序统计量，$X_{(n)}$ 称为最大次序统计量。


### 二、次序统计量的分布
设总体为连续分布,分布函数为$F(x)$,概率密度函数为$p(x)$则：

* $(X_{(1)},X_{(2)},...,X_{(n)})$ 的联合密度函数为:
    $$g(y_{1},y_{2},\cdots,y_{n})=n!p(y_{1})\cdots p(y_{n}) \quad y_{1}\le...\le y_{n}$$
* 最大次序统计量 $X_{(n)}$ 的分布函数为 $F^{n}(y)$，密度为 $np(y)F^{n-1}(y)$。
* 最小次序统计量 $X_{(1)}$ 的分布函数为 $1-(1-F(y))^{n}$，密度为 $np(y)(1-F(y))^{n-1}$。
* $X_{(k)}$ 的密度函数为:
    $$\frac{n!}{(k-1)!(n-k)!}p(y)F^{k-1}(y)(1-F(y))^{n-k}$$
* $(X_{(i)},X_{(j)})(i<j)$ 的联合密度函数为:
    $$g_{ij}(y_{i},y_{j})=\frac{n!}{(i-1)!(j-i-1)!(n-j)!}p(y_{i})p(y_{j}) \times F^{i-1}(y_{i})(F(y_{j})-F(y_{i}))^{j-i-1}(1-F(y_{j}))^{n-j}, \quad y_i \le y_j$$
* 特别地，$(X_{(1)},X_{(n)})$ 的联合密度函数为:
    $$g_{1n}(y_{1},y_{n})=n(n-1)(F(y_{n})-F(y_{1}))^{n-2}p(y_{1})p(y_{n}), \quad y_1 \le y_n$$

**定理**: $X_{(j_{1})},X_{(j_{2})}\cdots,X_{(j_{r})}$(其中 $j_{1}<j_{2}<\cdots<j_{r}$) 的联合密度函数为
$$g(y_{j_{1}},y_{j_{2}},\cdots,y_{j_{r}}) = \frac{n!}{(j_{1}-1)!(j_{2}-j_{1}-1)!\cdots(j_{r}-j_{r-1}-1)!(n-j_{r})!} \times p_{1}^{j_{1}-1}p_{2}^{j_{2}-j_{1}-1}\cdots p_{r}^{j_{r}-j_{r-1}-1}p_{r+1}^{n-j_{r}}\times p(y_{j_{1}})p(y_{j_{2}})\cdots p(y_{j_{r}})$$
其中 $a=y_{j_0} \le \cdots \le y_{j_{r+1}}=b$，$p_k = F(y_{j_k}) - F(y_{j_{k-1}})$。

### 三、极差、样本中位数与分位数

* **样本极差 (sample range)**: $R_{n}=X_{(n)}-X_{(1)}$。其密度函数为：
    $$p_{R_{n}}(r)=\int_{-\infty}^{\infty}n(n-1)(F(r+z)-F(z))^{n-2}p(z)p(r+z)dz, \quad r>0$$
* **样本中位数 (sample median)**:
    $$m_{e}=\begin{cases}X_{(\frac{n+1}{2})} & \text{当n为奇数} \\ (X_{(\frac{n}{2})}+X_{(\frac{n}{2}+1)})/2 & \text{当n为偶数} \end{cases}$$
* **总体分位数**: 对于$0<p<1$，若 $F(\xi_{p}-0)\le p$ 但 $F(\xi_{p})\ge p$，则称 $\xi_p$ 为总体X的p分位数。
* **样本分位数**: 
    **定义**: 设 $X_{1},X_{2},\cdots,X_{n}$ 是取自总体 $F(x)$ 的一个样本，$0<p<1$, 称
    $$m_{p}=X_{([np])}+(n+1)(p-\frac{[np]}{n+1})(X_{([np]+1)}-X_{([np])})$$
    为该样本的(下侧)p分位数。
* **样本四分位数 (quartile)**:
    * $m_{0.25}$：下四分位数
    * $m_{0.75}$：上四分位数
    * $IQR=m_{0.75}-m_{0.25}$：四分位距(interquartile range)

**定理**: 设总体的p分位数为 $\xi_p$，又设总体的密度函数 $f(x)$ 在 $\xi_p$ 处连续且不为零, 那么当 $n\rightarrow\infty$ 时, 有:
$$\sqrt{n}(m_{p}-\xi_{p})\Rightarrow N(0,\frac{p(1-p)}{f^{2}(\xi_{p})})$$
即
$$m_{p}\sim AN(\xi_{p},\frac{p(1-p)}{nf^{2}(\xi_{p})})$$
特别地，样本中位数 $m_e$:
$$m_{e}\sim AN(\xi_{0.5},\frac{1}{4nf^{2}(\xi_{0.5})})$$

---


## §2.4 常用分布与分布族

### Γ 分布 (Gamma Distribution)

* **定义**
    具有下列密度函数的分布称为Γ分布：
    $$p(x;\alpha,\lambda)=\frac{\lambda^{\alpha}}{\Gamma(\alpha)}x^{\alpha-1}e^{-\lambda x}, \quad x>0$$
    其中 $α > 0$ 称为“形状参数”，$λ > 0$ 称为“速率参数” 。

* **性质**
    * **分布函数**：当 $x > 0$ 时，
        $$F(x;\alpha,\lambda)=\int_{0}^{x}\frac{\lambda^{\alpha}y^{\alpha-1}}{\Gamma(\alpha)}e^{-\lambda y}dy=\int_{0}^{\lambda x}\frac{t^{\alpha-1}}{\Gamma(\alpha)}e^{-t}dt = F(\lambda x;\alpha,1)$$ 
    * **矩母函数**：$Ee^{tX}=(1-\frac{t}{\lambda})^{-\alpha}$, $t<\lambda$ 。
    * **特征函数**：$\varphi(t)=Ee^{itX}=(1-\frac{it}{\lambda})^{-\alpha}$ 。

* **数字特征**
    * **k阶矩**：$EX^{k}=\frac{\Gamma(k+\alpha)}{\lambda^{k}\Gamma(\alpha)}$ 。
    * **负指数阶矩**：$EX^{-\beta}=\frac{\lambda^{\beta}\Gamma(\alpha-\beta)}{\Gamma(\alpha)},$ $\beta<\alpha$ 。
    * **均值**：$E(\Gamma(\alpha,\lambda))=\frac{\alpha}{\lambda}$ 。
    * **方差**：$Var(\Gamma(\alpha,\lambda))=\frac{\alpha}{\lambda^{2}}$ 。

* **定理**
    1.  设随机变量 $X_1, X_2$ 相互独立, $X_{i}\sim\Gamma(\alpha_{i},\lambda)$, $i=1,2,$ 则 $X_{1}+X_{2}\sim\Gamma(\alpha_{1}+\alpha_{2},\lambda)$ 。
    2.  设 $X\sim\Gamma(\alpha,\lambda)$ 则 $Y=X/k\sim\Gamma(\alpha,k\lambda)$,其中 $k>0$ 。

* **与其他分布的关系**
    * $\Gamma(1,\lambda)$ 为指数分布 $E(\lambda)$ 。
    * $\Gamma(n/2,1/2)$ 就是 $\chi^{2}(n)$ 分布 。



### χ² 分布 (Chi-squared Distribution)

* **构造性定义**
    设 $X_{1},X_{2},...,X_{n}$ i.i.d, $\sim N(0,1)$,则 $\xi=\sum_{i=1}^{n}X_{i}^{2}$ 的分布定义为具有自由度n的χ²分布, 记作 $\xi\sim\chi^{2}(n)$ 。

* **性质**
    * **密度函数**：$p(x;n)=\frac{1}{2^{n/2}\Gamma(n/2)}x^{n/2-1}e^{-x/2},$ $x>0$ 。
    * **可加性**: 若 $\xi_{1}\sim\chi^{2}(n_{1}),$ $\xi_{2}\sim\chi^{2}(n_{2}),$ 且相互独立, 则 $\xi_{1}+\xi_{2}\sim\chi^{2}(n_{1}+n_{2})$ 。
    * **特征函数**：$\varphi(t;n)=(1-2it)^{-n/2}$ 。

* **数字特征**
    * **均值**：$E\xi=n$ 。
    * **方差**：$Var\xi=2n$ 。

* **与其他分布的关系**
    * χ²分布是Γ分布的一个特例：$\chi^{2}(n)$ 分布是 $\Gamma(n/2,1/2)$ 。
    * 若 $X_1 \sim N(0,1)$，则 $Y=X_1^2 \sim \Gamma(1/2, 1/2)$，即 $\chi^2(1)$ 。



### 非中心 χ² 分布 (Non-central Chi-squared Distribution)

* **构造性定义**
    设 $X_{1},...,X_{n}$ 相互独立, $X_{i}\sim N(a_{i},1)$，$a_i$ 不全为0。则 $\xi=\sum_{i=1}^{n}X_{i}^{2}$ 的分布定义为具有n自由度、非中心参数为 $\lambda=a_{1}^{2}+\cdots+a_{n}^{2}$ 的非中心χ²分布, 记为 $\sim\chi^{2}(n,\lambda)$ 。

* **性质**
    * **特征函数**：$\varphi(t;n)=(1-2it)^{-n/2}exp\{\frac{it\lambda}{1-2it}\}$ 。
    * **可加性**：若 $\xi_{j}\sim\chi_{n_{j},\lambda_{j}}^{2},j=1,...,k,$ 且相互独立, 则 $\sum_{j=1}^{k}\xi_{j}\sim\chi_{n,\lambda}^{2}$，其中 $n=\sum_{j=1}^{k}n_{j},\lambda=\sum_{j=1}^{k}\lambda_{j}$ 。

* **数字特征**
    * **均值**：$E\xi=n+\lambda$ 。
    * **方差**：$Var(\xi)=2(n+2\lambda)$ 。

---

### t 分布 (Student's t Distribution)

* **构造性定义**
    设 $X\sim N(0,1)$, $K\sim\chi^{2}(n)$, 且X与K相互独立,则 $T=\frac{X}{\sqrt{K/n}}$ 的分布定义为自由度为n的t分布, 记作 $T\sim t(n)$ 。

* **性质**
    * **密度函数**：$p(t;n)=\frac{\Gamma((n+1)/2)}{\sqrt{n\pi}\Gamma(n/2)}(1+\frac{t^{2}}{n})^{-(n+1)/2},-\infty<t<\infty$ 。
    * **尾部特征**：t分布的尾部比标准正态分布的“尾重” 。
    * **极限分布**：当 $n\rightarrow\infty$ 时, t(n)依分布收敛到 $N(0,1)$ 分布 。

* **数字特征**
    * **矩**：只有当 $r<n$ ($n>1$) 时, r阶矩才存在 。奇数阶矩为0，偶数阶矩为：
        $$ET^{r}=EN(0,1)^{r}E(\chi^{2}(n)/n)^{-r/2} = \frac{n^{\frac{r}{2}}}{\sqrt{\pi}}\frac{\Gamma(\frac{r+1}{2})\Gamma(\frac{n-r}{2})}{\Gamma(\frac{n}{2})}, \quad r<n \text{ 为偶数}$$ 
    * **均值**：$E(t(n))=0$，$n\ge2$ 。
    * **方差**：$Var(t(n))=\frac{n}{n-2},n\ge3$ 。

* **定理**
    1.  **单样本T统计量**：设 $X_{1},...,X_{n}$ 为来自总体 $X\sim N(\mu,\sigma^{2})$ 的样本, $\overline{X}$ 为样本均值, S² 为样本方差,则 
    $$T=\frac{\sqrt{n}(\overline{X}-\mu)}{S}\sim t(n-1)$$
    2.  **双样本T统计量**：设 $X_{1},...,X_{m}$ i.i.d $\sim N(\mu_{X},\sigma^{2})$，$Y_{1},...,Y_{n}$ i.i.d. $\sim N(\mu_{Y},\sigma^{2})$ 且两组样本相互独立。记 $S_{W}^{2}=\frac{1}{m+n-2}\{(m-1)S_{X}^{2}+(n-1)S_{Y}^{2}\}$ ，则
        $$\frac{\overline{X}-\overline{Y}-(\mu_{X}-\mu_{Y})}{S_{W}\sqrt{\frac{1}{m}+\frac{1}{n}}}\sim t(m+n-2)$$



### 非中心 t 分布 (Non-central t Distribution)

* **构造性定义**
    设 $X\sim N(\delta,1),$ $K\sim\chi^{2}(n)$ 且X与K相互独立,则 $T=\frac{X}{\sqrt{K/n}}$ 的分布定义为自由度为n、非中心参数为δ的非中心t分布, 记为 $T\sim t(n,\delta)$。

---

### F 分布 (F Distribution)

* **构造性定义**
    设 $K_{1}\sim\chi^{2}(m)$, $K_{2}\sim\chi^{2}(n)$, 且K₁与K₂相互独立,则 $F=\frac{K_{1}/m}{K_{2}/n}$ 的分布定义为具有自由度(m,n)的F分布, 记为 $F\sim F(m,n)$ 。

* **性质**
    - **密度函数**：$p(x;m,n)=\frac{\Gamma(\frac{m+n}{2})}{\Gamma(\frac{m}{2})\Gamma(\frac{n}{2})}m^{m/2}n^{n/2}\frac{x^{m/2-1}}{(n+mx)^{(m+n)/2}},x>0$ 。
    * **倒数性质**：若 $F\sim F(n,m)$, 则 $1/F\sim F(m,n)$ 。

* **数字特征**
    * **均值**：$EF(n,m)=\frac{m}{m-2}$ 。

* **定理**
    * 设 $X_1, ..., X_m$ i.i.d. $\sim N(\mu_X, \sigma_X^2)$ 和 $Y_1, ..., Y_n$ i.i.d. $\sim N(\mu_Y, \sigma_Y^2)$ 为两组独立的样本，其样本方差分别为 $S_X^2$ 和 $S_Y^2$，则 
    $$F=\frac{{S_{X}}^{2}/\sigma_{X}^{2}}{{S_{Y}}^{2}/\sigma_{Y}^{2}}\sim F(m-1,n-1)$$ 



### 非中心 F 分布 (Non-central F Distribution)

* **构造性定义**
    设 $K_{1}\sim\chi^{2}(m,\lambda),$ $K_{2}\sim\chi^{2}(n),$ 且K₁与K₂相互独立,则 $F=\frac{K_{1}/m}{K_{2}/n}$ 的分布定义为具有自由度(m,n)和非中心参数为λ的非中心F分布, 记为 $F\sim F(m,n,\lambda)$ 。

---

### β 分布 (Beta Distribution)

* **定义**
    具有下列密度函数的分布称为beta分布：
    $$p(x;a,b)=\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}x^{a-1}(1-x)^{b-1}, \quad 0<x<1$$ 
    其中a>0, b>0是形状参数 。

* **数字特征**
    * **k阶矩**：$EX^{k}=\frac{B(a+k,b)}{B(a,b)}=\frac{\Gamma(a+b)\Gamma(a+k)}{\Gamma(a)\Gamma(a+b+k)}$ 。
    * **均值**：$E(beta(a,b))=\frac{a}{a+b}$ 。
    * **方差**：$Var(beta(a,b))=\frac{ab}{(a+b)^{2}(a+b+1)}$ 。

* **与其他分布的关系**
    * 若a=b=1, 则为(0,1)区间上的均匀分布U(0,1) 。
    * 若a=b=1/2, 则为反正弦分布 。
    * 若 $X_1, ..., X_n$ i.i.d. $\sim U(0,1)$，则其次序统计量 $X_{(k)}\sim\beta(k,n-k+1)$ 
    * 若 $Y\sim\beta(a,b)$ 则 $X=\frac{Y}{1-Y}\sim Z(a,b)$ 
    * 若$X\sim\beta(a,b)$ 则 $Y=\frac{b}{a} \frac{X}{1-X} \sim F(2a,2b)$ （上条 + Z分布第二条）
    * 若 $X\sim Z(a,b)$ 则 $Y=\frac{X}{1+X}\sim\beta(a,b)$ 
    * 若 $X_1 \sim \Gamma(\alpha_1, \lambda)$ 和 $X_2 \sim \Gamma(\alpha_2, \lambda)$ 相互独立，则 $\frac{X_{1}}{X_{1}+X_{2}}\sim\beta(\alpha_{1},\alpha_{2})$ 
    * 若 $F\sim F(n,m)$，则
    $$\frac{\frac{n}{m}F}{1+\frac{n}{m}F}\sim\beta(n/2,m/2)$$

---

### Fisher Z 分布

* **定义**
    具有下列密度函数的分布称为Z分布：
    $$p(x;a,b)=\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}\frac{x^{a-1}}{(1+x)^{a+b}}, \quad x>0$$ 
    其中 a>0, b>0 。记为 $Z(a,b)$ 。

* **数字特征**
    * **k阶矩**：$EX^{k}=\frac{B(a+k,b-k)}{B(a,b)}=\frac{\Gamma(a+k)}{\Gamma(a)}\cdot\frac{\Gamma(b-k)}{\Gamma(b)}$ 。
    * **均值**：$EX=\frac{a}{b-1},b>1$ 。
    * **二阶矩**：$EX^{2}=\frac{(a+1)a}{(b-1)(b-2)}, b>2$ 。

* **与其他分布的关系**
    * **与Γ分布**: 设 $X_{1}\sim\Gamma(\alpha_{1},\lambda),X_{2}\sim\Gamma(\alpha_{2},\lambda)$ 且相互独立,则 $Y_{2}=X_{1}/X_{2}\sim Z(\alpha_{1},\alpha_{2})$ 。
    - **与F分布**: 若 $F\sim F(n,m)$, 则 $\frac{n}{m}F\sim Z(\frac{n}{2},\frac{m}{2})$  , $Z\sim Z(a,b)$, 则$\frac{b}{a}Z\sim F(2a,2b)$
    * **与β分布**: 见β分布部分 。

---

### 上α分位数

* **定义**
    设X的概率密度函数为f(x),对于给定的正数α, $0<\alpha<1,$ 若存在实数 $x_\alpha$ 满足 $P(X\ge x_{\alpha})=\int_{x_{\alpha}}^{+\infty}f(x)dx=\alpha$，则称点 $x_\alpha$ 为X的上α分位点 。

* **常用分布的上α分位数**
    * **标准正态分布**：$u_\alpha$ 满足 $P\{X\ge u_{\alpha}\} = 1 - \Phi(u_\alpha) = \alpha$ 。
    * **χ²(n)分布**：$\chi_\alpha^2(n)$ 满足 $\int_{\chi_{\alpha}^{2}(n)}^{\infty}f(x;n)dx=\alpha$ 。
    * **t(n)分布**：$t_\alpha(n)$ 满足 $\int_{t_{\alpha}(n)}^{\infty}f(x;n)dx=\alpha$ 。其具有对称性 $t_{1-\alpha}(n)=-t_{\alpha}(n)$ 。
    * **F(m,n)分布**：$F_\alpha(m,n)$ 满足 $\int_{F_{\alpha}(m,n)}^{\infty}f(x;m,n)dx=\alpha$ 。其具有倒数关系 $F_{1-\alpha}(m,n)=\frac{1}{F_{\alpha}(n,m)}$ 。

## §2.6 指数型分布族 (Exponential family)

**定义**: 设有参数分布族 $\mathcal{F}=\{p(x;\theta):\theta\in\Theta\}$，$p(x;\theta)$ 为分布的密度函数(pdf)或分布列(pmf)。若 $p(x;\theta)$ 可表示成如下形式
$$p(x;\theta)=c(\theta)exp\{\sum_{j=1}^{k}Q_{j}(\theta)T_{j}(x)\}h(x)$$
则称此分布族称为指数型分布族，或简称为指数族。其中 $c(\theta)>0, Q_{j}(\theta)$ $(j=1,...,k)$ 为定义在参数空间上的函数，与x无关，$h(x)\ge0$, $T_{j}(x)$ $(j=1,...,k)$ 为与θ无关的函数。

如果令 $\lambda_{j}=Q_{j}(\theta)$, 若$c(\theta)$可表示成 $\tilde{\lambda}=(\lambda_{1},...,\lambda_{k})$ 的函数 $c^{*}(\tilde{\lambda})$，那么 $p(x;\theta)$ 可表示成
$$p(x;\tilde{\lambda})=c^{*}(\tilde{\lambda})exp\{\sum_{j=1}^{k}\lambda_{j}T_{j}(x)\}h(x)$$
这种形式称为指数型分布族的**自然形式** (natural form)，此时
$$\Lambda=\{\tilde{\lambda}:\int exp\{\sum_{j=1}^{k}\lambda_{j}T_{j}(x)\}h(x)dx<\infty\}$$
称为**自然参数空间** (natural parametric space).

**定理**: 如果总体分布族是指数型分布族，那么从中抽取的简单随机样本的分布族也是指数型分布族。因为如果总体X的pdf或pmf为 $c(\theta)exp\{\sum_{j=1}^{k}Q_{j}(\theta)T_{j}(x)\}h(x)$，那么样本 $(X_{1},X_{2},...,X_{n})$ 的联合pdf或pmf为
$$c^{n}(\theta)exp\{\sum_{j=1}^{k}Q_{j}(\theta)\{T_{j}(x_{1})+...+T_{j}(x_{n})\}\}h(x_{1})\cdots h(x_{n})$$

**结论**: 指数族分布的支撑集与参数θ无关。由定义不难看出，对于指数族，支撑集 $\{x:p(x;\theta)>0\}=\{x:h(x)>0\}$，显然与θ无关。

---

## §2.7 充分统计量 (Sufficient Statistics)

**定义**: 设有一个分布族 $\mathcal{F}=\{F_\theta\}$，$(X_{1},X_{2},\cdots,X_{n})$ 是从某总体$F \in \mathcal{F}$中抽取的一个样本。$T=T(X_{1},X_{2},...,X_{n})$ 为一个(一维或多维的)统计量，如果当给定 $T=t$ 后, 样本 $(X_{1},X_{2},...,X_{n})$ 的条件分布与总体分布F(或参数θ)无关,则称T为此分布族的**充分统计量** (sufficient statistics)。

**定理 (因子分解定理)**: 设样本X的联合pdf或pmf为 $p(\tilde{x};\theta)$，其中θ为未知参数, 则$T=T(\tilde{X})$为(关于θ的)充分统计量，当且仅当
$$p(\tilde{x};\theta)=g(T(\tilde{x});\theta)h(\tilde{x})$$
其中 $g(t;\theta)$ 是定义在统计量 $T(\tilde{X})$ 取值空间上的函数，$h(\tilde{x})$ 与θ无关。

**结论**: 假定样本的pdf或pmf所在的分布族是指数型分布族，即
$$p(\tilde{x};\theta)=c^{*}(\theta)exp\{\sum_{j=1}^{k}Q_{j}(\theta)T_{j}^{*}(\tilde{x})\}h^{*}(\tilde{x})$$
则 $\tilde{T}=(T_{1}^{*},...,T_{k}^{*})$ 为充分统计量。

**定义 (极小充分统计量)**: 设S是分布族 $\mathcal{F}$ 的充分统计量，假如对 $\mathcal{F}$ 的任意一个充分统计量T,存在一个函数 $f(\cdot)$，使得 $S=f(T)$，则称S是此分布族$\mathcal{F}$的**极小充分统计量**。

**定理**: 设总体是指数型分布族，则样本的联合pdf或pmf可以写为
$$p(\tilde{x};\theta)=c^{*}(\theta)exp\{\sum_{j=1}^{k}Q_{j}(\theta)T_{j}^{*}(\tilde{x})\}h^{*}(\tilde{x})$$
若 $\tilde{Q}(\theta)=(Q_{1}(\theta),\cdots,Q_{k}(\theta))$ 的值域有非空的内部, 则 $\tilde{T}=(T_{1}^{*},...,T_{k}^{*})$ 为极小充分统计量。

---

### 常见样本分布的充分统计量

* **二项分布** $B(1,p)$: 对来自该总体的样本 $X_1, ..., X_n$，统计量 $T = \sum_{i=1}^n X_i$ 是参数 $p$ 的充分统计量。
* **几何分布** $P(X=x)=\theta(1-\theta)^{x}$: 对来自该总体的样本 $X_1, ..., X_n$，统计量 $T=\sum_{i=1}^{n}X_{i}$ 是参数 $\theta$ 的充分统计量。
* **指数分布** $E(\lambda)$: 对来自该总体的样本 $X_1, ..., X_n$，样本均值 $\overline{X}$ 是参数 $\lambda$ 的充分统计量。
* **均匀分布** $U(0,\theta)$: 对来自该总体的样本 $X_1, ..., X_n$，最大次序统计量 $X_{(n)}$ 是参数 $\theta$ 的充分统计量。
* **正态分布** $N(\mu,\sigma^2)$: 对来自该总体的样本 $X_1, ..., X_n$，
    * 若 $\mu, \sigma^2$ 均未知, $(\sum_{i=1}^{n}X_{i}, \sum_{i=1}^{n}X_{i}^{2})$ 是 $(\mu, \sigma^2)$ 的充分统计量，这等价于 $(\overline{X}, S^2)$ 是充分统计量。
* **任意分布**: 对来自任意总体的样本 $X_1, ..., X_n$，次序统计量 $(X_{(1)}, ..., X_{(n)})$ 是充分统计量。