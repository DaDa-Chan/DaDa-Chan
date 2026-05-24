---
title: "Regression, Classification, and ERM"
date: 2026-01-20
category: "STAT 154"
excerpt: "Statistical learning framework, ERM, linear regression, binary classification, logistic regression, LDA/QDA, and GLMs. Based on Hugo Chardon's STAT 154/254 lecture notes, Spring 2026 UC Berkeley."
---

<div style="display:none">
$$
\newcommand{\R}{\mathbb{R}}
\newcommand{\N}{\mathbb{N}}
\newcommand{\E}{\mathbb{E}}
\newcommand{\Var}{\mathrm{Var}}
\newcommand{\Tr}{\mathrm{Tr}}
\newcommand{\rank}{\mathrm{rank}}
\newcommand{\sign}{\mathrm{sign}}
\newcommand{\KL}{\mathrm{KL}}
\newcommand{\iid}{\overset{\mathrm{i.i.d.}}{\sim}}
\newcommand{\asto}{\xrightarrow{\;\mathrm{a.s.}\;}}
\newcommand{\ip}[1]{\langle #1 \rangle}
\newcommand{\norm}[1]{\| #1 \|}
\newcommand{\abs}[1]{| #1 |}
\newcommand{\what}[1]{\widehat{#1}}
\newcommand{\Flin}{\mathcal{F}_{\mathrm{lin}}}
\newcommand{\Fclass}{\mathcal{F}}
\newcommand{\argmin}{\operatorname*{arg\,min}}
\newcommand{\argmax}{\operatorname*{arg\,max}}
$$
</div>

*Based on Hugo Chardon's handwritten lecture notes, STAT 154/254 Spring 2026, UC Berkeley. Reorganised with errata corrections.*

---

## Notation

<div class="note-rem" style="padding:12px 16px; font-size:0.93em;">
<span class="note-label" style="color:#7b34a6;">Convention on $\E$</span>
Throughout, $\E$ without subscript denotes expectation w.r.t.\ the joint distribution $P$ of $(X,Y)$, unless stated otherwise. When only a subset of variables is random (e.g.\ fixed-design regression), we write $\E_Y[\cdot]$.
</div>

| Symbol | Meaning |
|--------|---------|
| $(X,Y)\sim P$ | Generic test pair; $P$ unknown on $\mathcal{X}\times\mathcal{Y}$ |
| $(X_1,Y_1),\dots,(X_n,Y_n)\iid P$ | Training sample |
| $\ell(\hat y, y)\geq 0$ | Loss function |
| $L(f) = \E_{(X,Y)\sim P}[\ell(f(X),Y)]$ | Population risk of predictor $f$ |
| $\what L_n(f)=\frac{1}{n}\sum_{i=1}^n\ell(f(X_i),Y_i)$ | Empirical risk |
| $f_{\mathrm B}$, $L_{\mathrm B}$ | Bayes predictor and Bayes risk |
| $f^{\ast}\in\argmin_{f\in\Fclass}L(f)$ | Best predictor in model $\Fclass$ |
| $\what f_n\in\argmin_{f\in\Fclass}\what L_n(f)$ | ERM predictor |
| $\Sigma=\E[XX^\top]$ | Second-moment matrix of $X$ |
| $\what\Sigma_n=\frac{1}{n}\sum_i X_iX_i^\top$ | Empirical second-moment matrix |
| $\mathbb{X}\in\R^{n\times d}$ | Design matrix (rows $=X_i^\top$) |
| $\eta(x)=\E[Y\mid X=x]$ | Regression function |
| $\sigma(t)=1/(1+e^{-t})$ | Sigmoid / logistic function |
| $\mathbf 1(\cdot)$ | Indicator function |

---

## 1. Statistical Learning Framework and Linear Regression

### 1.1 The Supervised Learning Problem

<div class="note-box note-def">
<span class="note-label">Definition 1.1 — Supervised learning problem</span>

Given: (i) an input space $\mathcal X$ (typically $\R^d$) and outcome space $\mathcal Y$; (ii) a loss function $\ell:\mathcal Y\times\mathcal Y\to[0,\infty)$; (iii) a distribution $P$ on $\mathcal X\times\mathcal Y$ generating i.i.d.\ training data $(X_1,Y_1),\dots,(X_n,Y_n)$; (iv) a model class $\Fclass\subset\{f:\mathcal X\to\mathcal Y\}$;

design an algorithm that uses the training data to construct $\what f_n\in\Fclass$ with small **population risk**

$$L(\what f_n) = \E_{(X,Y)\sim P}\bigl[\ell(\what f_n(X),Y)\bigr].$$

Here $(X,Y)$ is an *independent test pair* drawn from the same $P$.
</div>

**Common losses.**

- *Square loss* (regression): $\ell(\hat y,y)=(\hat y - y)^2$.
- *0/1 loss* (classification): $\ell(\hat y,y)=\mathbf 1(\hat y\ne y)$.
- *Logistic loss*: $\ell(y,\hat y)=\log(1+e^{-y\hat y})$, with $y\in\{-1,1\}$.
- *Hinge loss*: $\ell(y,\hat y)=(1-y\hat y)_+$.

### 1.2 Excess Risk Decomposition

<div class="note-box note-def">
<span class="note-label">Definition 1.2 — Bayes risk and Bayes predictor</span>

The **Bayes risk** is $L_{\mathrm B}=\inf_{\text{all }f}L(f)$. Any minimiser $f_{\mathrm B}$ is called a **Bayes predictor**.
</div>

For any $f\in\Fclass$, the excess risk decomposes as

$$\begin{equation}\label{eq:bias-var}
\underbrace{L(f)-L_{\mathrm B}}_{\text{excess risk}}
\;=\;
\underbrace{L(f)-\inf_{g\in\Fclass}L(g)}_{\text{estimation error}}
\;+\;
\underbrace{\inf_{g\in\Fclass}L(g)-L_{\mathrm B}}_{\text{approximation error}}.
\end{equation}$$

As $\Fclass$ grows, estimation error increases while approximation error decreases: this is the **bias–variance trade-off**.

<div class="note-box note-def">
<span class="note-label">Definition 1.3 — Well-specified model</span>

$\Fclass$ is **well-specified** if $f_{\mathrm B}\in\Fclass$; otherwise it is **misspecified**.
</div>

### 1.3 Empirical Risk Minimisation (ERM)

The population risk $L(f)=\E[\ell(f(X),Y)]$ cannot be computed (since $P$ is unknown). The **ERM paradigm** replaces it by

$$\begin{equation}\label{eq:emp-risk}
\what L_n(f) = \frac{1}{n}\sum_{i=1}^n \ell(f(X_i),Y_i),
\end{equation}$$

and minimises over $\Fclass$: $\what f_n\in\argmin_{f\in\Fclass}\what L_n(f)$.

<div class="note-box note-ex">
<span class="note-label">Example — MLE as ERM</span>

Parametric model $\{p_\theta:\theta\in\Theta\}$. The MLE $\what\theta_n = \argmin_\theta\frac{1}{n}\sum_i(-\log p_\theta(Z_i))$ is ERM with loss $\ell(\theta,z)=-\log p_\theta(z)$. The excess risk is $\KL(P_{\theta^{\ast}}\|P_\theta)$.
</div>

<div class="note-box note-ex">
<span class="note-label">Example — Mean estimation as ERM</span>

For $Z\sim P$ on $\R^d$, $\E[Z]=\argmin_{m\in\R^d}\E[\norm{Z-m}^2]$. The ERM is the sample mean $\bar Z_n=\frac{1}{n}\sum_i Z_i$.
</div>

### 1.4 Linear Algebra Essentials

**Inner product and norm.** $\ip{u,v}=\sum_{j=1}^d u_jv_j$, $\norm{u}=\sqrt{\ip{u,u}}$. **Pythagoras:** $\norm{u+v}^2=\norm{u}^2+\norm{v}^2+2\ip{u,v}$.

**Transpose identity.** For $A\in\R^{p\times q}$, $u\in\R^q$, $v\in\R^p$: $\ip{Au,v}=\ip{u,A^\top v}$.

**Frobenius.** $\ip{A,B}_F=\Tr(AB^\top)$, $\norm{A}_F=\sqrt{\Tr(AA^\top)}$.

**Spectral theorem.** Any real symmetric $M$ admits $M=V\Lambda V^\top$ with $V$ orthogonal and $\Lambda=\mathrm{diag}(\lambda_1,\dots,\lambda_d)$.

<div class="note-box note-def">
<span class="note-label">Definition 1.4 — Elliptic norm</span>

For $M\succ 0$, $\norm{u}_M=\sqrt{u^\top M u}=\norm{M^{1/2}u}$. In the eigenbasis of $M$: $\norm{u}_M^2=\sum_j\lambda_j\tilde u_j^2$, where $\tilde u=V^\top u$.
</div>

### 1.5 Linear Regression

**Setup.** $\mathcal X=\R^d$, $\mathcal Y=\R$, square loss, linear class $\Flin=\{f_\theta:x\mapsto\ip{\theta,x},\;\theta\in\R^d\}$. Write $L(\theta)\coloneqq L(f_\theta)=\E_{(X,Y)\sim P}[(Y-\ip{\theta,X})^2]$.

#### The Optimal Linear Predictor

<div class="note-box note-prop">
<span class="note-label">Proposition 1.1 — Optimal linear predictor</span>

If $\Sigma=\E[XX^\top]\succ 0$, then $L(\cdot)$ is strictly convex and the unique minimiser is

$$\theta^{\ast} = \Sigma^{-1}\E[YX].$$
</div>

<div class="note-proof">
<em>Proof.</em> The gradient is $\nabla L(\theta) = -2\,\E\bigl[(Y-\ip{\theta,X})\,X\bigr]$. The Hessian is $\nabla^2 L(\theta)=2\Sigma\succ 0$, so $L$ is strictly convex. Setting $\nabla L(\theta^{\ast})=0$ gives $\E[YX]=\Sigma\,\theta^{\ast}$. $\square$
</div>

<div class="note-box note-err">
<span class="note-label">Errata — Gradient of a quadratic form (original notes, p. 19)</span>

The original notes state: for $F(\theta)=\tfrac{1}{2}\ip{A\theta,\theta}$ with $A$ symmetric, "$\nabla F(\theta)=\tfrac{1}{2} A\theta$".

**Correction.** $\nabla F(\theta) = A\theta$ (not $\frac{1}{2}A\theta$). Verification: $F(\theta+h)=F(\theta)+\ip{A\theta,h}+\frac{1}{2}\ip{Ah,h}$, so $\nabla F(\theta)=A\theta$.
</div>

<div class="note-box note-rem">
<span class="note-label">Remark — Gradient rules used throughout</span>

1. **Linear:** $F(\theta)=\ip{u,\theta}\Rightarrow\nabla F=u$.
2. **Quadratic:** $F(\theta)=\frac{1}{2}\theta^\top A\theta\Rightarrow\nabla F(\theta)=A\theta$, $\nabla^2 F=A$ ($A$ symmetric).
3. **Chain rule:** $F=g\circ\varphi\Rightarrow\nabla F(\theta)=g'(\varphi(\theta))\,\nabla\varphi(\theta)$.
</div>

#### The OLS Estimator

Since $\theta^{\ast}$ involves unknown $\Sigma$ and $\E[YX]$, we replace expectations by sample averages:

$$\begin{equation}\label{eq:OLS}
\what\theta_n^{\,\mathrm{LS}}
= \what\Sigma_n^{-1}\cdot\frac{1}{n}\sum_{i=1}^n Y_i X_i
= (\mathbb X^\top\mathbb X)^{-1}\mathbb X^\top Y.
\end{equation}$$

**Geometric view.** The empirical risk is $\what L_n(\theta)=\frac{1}{n}\norm{Y-\mathbb X\theta}^2$. Minimising over $\theta$ is the orthogonal projection of $Y\in\R^n$ onto $\mathrm{Im}(\mathbb X)$.

<div class="note-box note-err">
<span class="note-label">Errata — "positive definite" vs "positive semidefinite" (original notes, p. 17)</span>

The notes state: "a symmetric matrix $A$ is *symmetric* if $u^\top Au\geq 0$ for all $u$."

**Correction.** Should read: "$A$ is **positive semidefinite** if $u^\top Au\geq 0$ for all $u$."
</div>

#### Conditional Expectation and the Bayes Predictor

<div class="note-box note-prop">
<span class="note-label">Proposition 1.2 — Bayes predictor under square loss</span>

If $\E[Y^2]<\infty$, the Bayes predictor under square loss is the regression function $\eta(x)=\E[Y\mid X=x]$. For any $f$:

$$\E\bigl[(Y-f(X))^2\bigr] = \E\bigl[(Y-\eta(X))^2\bigr] + \E\bigl[(\eta(X)-f(X))^2\bigr].$$
</div>

<div class="note-proof">
<em>Proof.</em> Write $Y-f(X) = (Y-\eta(X))+(\eta(X)-f(X))$. Expanding and taking $\E$, the cross term vanishes by the tower property: $\E\bigl[(Y-\eta(X))(\eta(X)-f(X))\bigr] = \E\bigl[(\eta(X)-f(X))\,\underbrace{\E[Y-\eta(X)\mid X]}_{=0}\bigr]=0$. $\square$
</div>

**Well-specification for linear regression.** $\Flin$ is well-specified iff $\eta(x)=\ip{\theta^{\ast},x}$ for some $\theta^{\ast}$, equivalently $Y=\ip{\theta^{\ast},X}+\varepsilon$ with $\E[\varepsilon\mid X]=0$.

#### Relevance of Linearity: Embeddings

If $Y=\sum_{j=1}^d\theta_j^{\ast}\varphi_j(X)$ for some basis $\varphi_1,\dots,\varphi_d$, define $\Phi(X)=(\varphi_1(X),\dots,\varphi_d(X))$; then $Y=\ip{\theta^{\ast},\Phi(X)}$ is linear *in the features* $\Phi(X)$. Examples: polynomial regression, Fourier bases, expert aggregation.

### 1.6 Risk Bounds for OLS

#### Excess Risk Has a Closed Form

<div class="note-box note-prop">
<span class="note-label">Proposition 1.3 — Excess risk = elliptic distance</span>

For any $\theta\in\R^d$:

$$\begin{equation}\label{eq:excess-elliptic}
L(\theta)-L(\theta^{\ast}) = \norm{\theta-\theta^{\ast}}_\Sigma^2.
\end{equation}$$
</div>

<div class="note-proof">
<em>Proof.</em> Expand $L(\theta)-L(\theta^{\ast}) = \E[\ip{\theta^{\ast}-\theta,X}^2] + 2\E[(Y-\ip{\theta^{\ast},X})\ip{\theta^{\ast}-\theta,X}]$. The first term is $\norm{\theta-\theta^{\ast}}_\Sigma^2$. The second equals $2\ip{\theta^{\ast}-\theta,\;\E[(Y-\ip{\theta^{\ast},X})X]}=0$ by the optimality condition. $\square$
</div>

**Consequence.** $L(\what\theta_n)-L(\theta^{\ast})=\norm{\what\theta_n-\theta^{\ast}}_\Sigma^2\asto 0$ by LLN, proving consistency.

#### Fixed-Design, Well-Specified Model

**Model.** $Y_i=\ip{\theta^{\ast},x_i}+\xi_i$, $x_i\in\R^d$ fixed, $\xi_i$ i.i.d.\ with $\E[\xi_i]=0$, $\E[\xi_i^2]=\sigma^2$. Let $\Sigma_n=\frac{1}{n}\mathbb X^\top\mathbb X$.

<div class="note-box note-thm">
<span class="note-label">Theorem 1.1 — Exact excess risk, fixed design</span>

If $\Sigma_n\succ 0$, the OLS satisfies

$$\begin{equation}\label{eq:fixed-design}
\E_Y\bigl[L(\what\theta_n)\bigr]-L(\theta^{\ast}) = \sigma^2\,\frac{d}{n}.
\end{equation}$$
</div>

<div class="note-proof">
<em>Proof.</em> Since $Y=\mathbb X\theta^{\ast}+\xi$, we have $\what\theta_n-\theta^{\ast}=(\mathbb X^\top\mathbb X)^{-1}\mathbb X^\top\xi$. By Proposition 1.3 (with $\Sigma_n$):

$$\E_Y[L(\what\theta_n)]-L(\theta^{\ast}) = \E_\xi\bigl[\norm{\what\theta_n-\theta^{\ast}}_{\Sigma_n}^2\bigr] = \frac{1}{n}\,\E_\xi\bigl[\xi^\top\Pi\,\xi\bigr],$$

where $\Pi=\mathbb X(\mathbb X^\top\mathbb X)^{-1}\mathbb X^\top$ is the projector onto $\mathrm{Im}(\mathbb X)$. Using $\E[\xi\xi^\top]=\sigma^2 I_n$: $\E[\xi^\top\Pi\xi]=\sigma^2\Tr(\Pi)=\sigma^2 d$. $\square$
</div>

#### Random-Design, Well-Specified Model

<div class="note-box note-prop">
<span class="note-label">Proposition 1.4 — Random design excess risk</span>

Under $Y=\ip{\theta^{\ast},X}+\xi$, $\E[\xi\mid X]=0$, $\E[\xi^2\mid X]=\sigma^2$, if $n\geq d$:

$$\E\bigl[L(\what\theta_n)\bigr]-L(\theta^{\ast}) = \frac{\sigma^2}{n}\;\E\bigl[\Tr(\what\Sigma_n^{-1}\Sigma)\bigr].$$

Under mild distributional assumptions this is of order $\sigma^2 d/n$.
</div>

#### Model Selection via Truncation

Consider $Y=\sum_{k=1}^\infty\theta_k^{\ast}\varphi_k(X)+\xi$ and truncated models $\Fclass_K$. Let $R_K=\sum_{k>K}(\theta_k^{\ast})^2$ (energy of truncated tail).

<div class="note-box note-err">
<span class="note-label">Errata — Approximation error (original notes, eq. (1.4.16))</span>

The original notes claim $\inf_{f\in\Fclass_K}L(f)-L(f^{\ast})=R_K+\sigma^2$.

**Correction.** The Bayes risk is $L(f^{\ast})=\sigma^2$ (irreducible noise), and the best risk in $\Fclass_K$ is $\sigma^2+R_K$. Hence the approximation error is

$$\inf_{f\in\Fclass_K}L(f)-L(f^{\ast}) = R_K,$$

*without* the $+\sigma^2$ term. The **total excess risk** is $R_K + (R_K+\sigma^2)\frac{K}{n}$, not $(R_K+\sigma^2)(1+\frac{K}{n})$.
</div>

### 1.7 Uniform Convergence and ERM

<div class="note-box note-prop">
<span class="note-label">Proposition 1.5 — ERM excess risk via uniform convergence</span>

$$\begin{equation}\label{eq:unif-conv}
L(\what f_n) - L(f^{\ast}) \;\leq\; 2\sup_{f\in\Fclass}\abs{\what L_n(f)-L(f)}.
\end{equation}$$
</div>

<div class="note-proof">
<em>Proof.</em> Write $L(\what f_n)-L(f^{\ast}) = [L(\what f_n)-\what L_n(\what f_n)] + \underbrace{[\what L_n(\what f_n)-\what L_n(f^{\ast})]}_{\leq\,0} + [\what L_n(f^{\ast})-L(f^{\ast})]$. The middle term is $\leq 0$ by definition of ERM; bound each remaining term by $\sup_{f\in\Fclass}|\what L_n(f)-L(f)|$. $\square$
</div>

<div class="note-box note-rem">
<span class="note-label">Remark — Why well-specified analysis is useful even under misspecification</span>

Since $\{\text{well-specified distributions}\}\subset\{\text{all distributions}\}$, a worst-case lower bound over well-specified distributions is also a lower bound in general.
</div>

---

## 2. Binary Classification

### 2.1 Setup and Bayes Classifier

Convention: $Y\in\{0,1\}$, loss $\ell_{0/1}(\hat y,y)=\mathbf 1(\hat y\ne y)$, risk $L(f)=P(f(X)\ne Y)$. The regression function is $\eta(x)=P(Y=1\mid X=x)$.

<div class="note-box note-thm">
<span class="note-label">Theorem 2.1 — Bayes classifier and Bayes risk</span>

1. $f_{\mathrm B}(x)=\mathbf 1(\eta(x)>1/2)$ (predict the most likely label).
2. $L_{\mathrm B}=\E[\min\{\eta(X),\,1-\eta(X)\}]$.
3. For any classifier $f$: $L(f)-L_{\mathrm B}=\E\bigl[|2\eta(X)-1|\cdot\mathbf 1(f(X)\ne f_{\mathrm B}(X))\bigr]$.
</div>

<div class="note-proof">
<em>Proof of (1) and (2).</em> For any $f:\mathcal X\to\{0,1\}$,

$$P(f(X)\ne Y\mid X) = \mathbf 1(f(X)=0)\,\eta(X) + \mathbf 1(f(X)=1)\,(1-\eta(X)) \geq \min\{\eta(X),\,1-\eta(X)\},$$

with equality iff $f=f_{\mathrm B}$. Taking $\E_X$ gives (2). $\square$
</div>

<div class="note-box note-err">
<span class="note-label">Errata — Incomplete sentence (original notes, p. 38)</span>

The notes state "if $X$ is such that it is concentrated around values $x$ where $\eta(x)$, then…" — should read "…where $\eta(x)\approx 1/2$, then even the Bayes risk is close to $1/2$."
</div>

### 2.2 Plug-in Classification via Regression

<div class="note-box note-def">
<span class="note-label">Definition 2.1 — Plug-in classifier</span>

Given an estimator $\what\eta_n$ of $\eta$, the **plug-in classifier** is $\what f_{\what\eta_n}(x) = \mathbf 1(\what\eta_n(x) > 1/2)$.
</div>

<div class="note-box note-prop">
<span class="note-label">Proposition 2.1 — Plug-in excess risk bound</span>

$$L(\what f_{\what\eta_n}) - L_{\mathrm B} \;\leq\; 2\,\E\bigl[|\what\eta_n(X)-\eta(X)|\bigr] \;\leq\; 2\sqrt{\E\bigl[(\what\eta_n(X)-\eta(X))^2\bigr]}.$$
</div>

<div class="note-proof">
<em>Proof.</em> By Theorem 2.1(3), the excess risk is $\E[|2\eta(X)-1|\cdot\mathbf 1(\what f\ne f_{\mathrm B})]$. When $\what f(x)\ne f_{\mathrm B}(x)$, the estimates are on opposite sides of $1/2$, so $|2\eta(x)-1|\leq 2|\what\eta_n(x)-\eta(x)|$. The second inequality is Jensen. $\square$
</div>

**Local averaging estimators** for $\eta$:

1. **Regressogram**: partition $\mathcal X$ into bins, average labels per bin.
2. **$k$-NN**: $\what\eta_{k\text{-NN}}(x)=\frac{1}{k}\sum_{i\in N_k(x)}Y_i$.
3. **Kernel smoothing (Nadaraya–Watson)**: $\what\eta_h(x) = \frac{\sum_i Y_i K(\frac{X_i-x}{h})}{\sum_i K(\frac{X_i-x}{h})}$.

All share a complexity parameter controlling the bias–variance trade-off.

### 2.3 PAC Bounds for Finite Models

<div class="note-box note-thm">
<span class="note-label">Theorem 2.2 — PAC bound for finite classes</span>

Let $\Fclass=\{f_1,\dots,f_M\}$. For any $\delta\in(0,1)$, with probability $\geq 1-\delta$:

$$\begin{equation}\label{eq:PAC-finite}
L(\what f_n)-L(f^{\ast}) \;\leq\; 2\sqrt{\frac{\log M+\log(2/\delta)}{2n}}.
\end{equation}$$
</div>

<div class="note-proof">
<em>Proof.</em> By Proposition 1.5 it suffices to bound $\sup_j|\what L_n(f_j)-L(f_j)|$. For each fixed $f_j$, Hoeffding's inequality gives $P(|\what L_n(f_j)-L(f_j)|>t)\leq 2e^{-2nt^2}$. Union bound over $M$ classifiers: $P(\sup_j|\what L_n(f_j)-L(f_j)|>t)\leq 2Me^{-2nt^2}$. Set $2Me^{-2nt^2}=\delta$ and solve for $t$. $\square$
</div>

**Key messages:** (i) Model complexity enters only as $\log M$; (ii) sample complexity to achieve excess risk $\leq\varepsilon$ is $n\gtrsim (\log M+\log(1/\delta))/\varepsilon^2$.

### 2.4 Logistic Regression

**Motivation.** ERM with 0/1 loss over linear classifiers is NP-hard. Replace by a **convex surrogate**.

<div class="note-box note-def">
<span class="note-label">Definition 2.2 — Logistic loss and logistic regression</span>

With $Y\in\{-1,+1\}$ and linear score $s_\theta(x)=\ip{\theta,x}$:

$$\ell_{\log}(y,s) = \log(1+e^{-ys}).$$

The logistic regression estimator is $\what\theta_n = \argmin_\theta \frac{1}{n}\sum_{i=1}^n\log(1+e^{-Y_i\ip{\theta,X_i}})$.
</div>

**Probabilistic interpretation.** The logistic model assumes $P(Y=1\mid X=x)=\sigma(\ip{\theta^{\ast},x})$. Under this model, $\what\theta_n$ is the MLE.

**Classification calibration.** The logistic loss is *classification-calibrated*: small excess logistic risk implies small excess classification risk (Bartlett–Jordan–McAuliffe, 2006).

### 2.5 LDA and QDA

**Setting.** $K$-class classification, $Y\in\{1,\dots,K\}$. **Generative model**: $X\mid Y=k\;\sim\;\mathcal N(\mu_k,\Sigma_k)$; prior $p_k=P(Y=k)$.

<div class="note-box note-def">
<span class="note-label">Definition 2.3 — LDA and QDA</span>

**LDA** (homoskedastic): $\Sigma_k=\Sigma$ for all $k$. **QDA** (heteroskedastic): each class has its own $\Sigma_k$.
</div>

<div class="note-box note-prop">
<span class="note-label">Proposition 2.2 — Discriminant functions</span>

The Bayes classifier is $f_{\mathrm B}(x)=\argmax_k\;\delta_k(x)$, where:

- **LDA:** $\delta_k(x) = \log p_k + \ip{\Sigma^{-1}\mu_k,\,x} - \frac{1}{2}\norm{\mu_k}_{\Sigma^{-1}}^2$ (linear in $x$).
- **QDA:** $\delta_k(x) = \log p_k - \frac{1}{2}\log\det\Sigma_k - \frac{1}{2}\norm{x-\mu_k}_{\Sigma_k^{-1}}^2$ (quadratic in $x$).
</div>

<div class="note-proof">
<em>Proof.</em> Apply Bayes' rule: $P(Y=k\mid X=x)\propto p_k\,\gamma_k(x)$ where $\gamma_k$ is the $\mathcal N(\mu_k,\Sigma_k)$ density. Take $\log$ and drop terms not depending on $k$. For LDA, the $-\frac{1}{2} x^\top\Sigma^{-1}x$ term is common to all classes and drops out. $\square$
</div>

**Estimation.** Replace $p_k$, $\mu_k$, $\Sigma_k$ by empirical versions: $\what p_k=|I_k|/n$, $\what\mu_k=\frac{1}{|I_k|}\sum_{i\in I_k}X_i$, $\what\Sigma_k=\frac{1}{|I_k|}\sum_{i\in I_k}(X_i-\what\mu_k)(X_i-\what\mu_k)^\top$.

**Complexity trade-off.** LDA: $O(d^2)$ parameters (shared $\Sigma$). QDA: $O(Kd^2)$ parameters — lower bias but higher variance; fragile in high dimensions.

**Naive Bayes.** Assumes $X^{(1)},\dots,X^{(d)}$ conditionally independent given $Y=k$: $f_{X|Y=k}(x)=\prod_{j=1}^d f_{k,j}(x^{(j)})$. Reduces to $d$ univariate problems (rate $n^{-2/3}$ instead of $n^{-2/d}$).

### 2.6 Poisson Regression and GLMs

<div class="note-box note-def">
<span class="note-label">Definition 2.4 — Poisson regression</span>

$Y\mid X\;\sim\;\mathrm{Poisson}(\lambda(X))$ with $\log\lambda(X)=\ip{\theta^{\ast},X}$. The MLE minimises

$$\what L_n(\theta) = \frac{1}{n}\sum_{i=1}^n\bigl[e^{\ip{\theta,X_i}}-Y_i\ip{\theta,X_i}\bigr] + \text{const}.$$
</div>

<div class="note-box note-def">
<span class="note-label">Definition 2.5 — Generalised linear model</span>

$Y\mid X$ belongs to an exponential family $p(y,\eta)=\exp(y\eta-A(\eta)+c(y))$, and the natural parameter is linear: $\eta(X)=\ip{\theta^{\ast},X}$. Equivalently, $g(\E[Y\mid X])=\ip{\theta^{\ast},X}$ where $g=(A')^{-1}$ is the **canonical link**.
</div>

| GLM | $\mathcal Y$ | Family | Link $g$ | $A(\eta)$ |
|-----|------|--------|----------|-----------|
| Linear regression | $\R$ | Gaussian | identity | $\eta^2/2$ |
| Logistic regression | $\{0,1\}$ | Bernoulli | logit | $\log(1+e^\eta)$ |
| Poisson regression | $\mathbb{N}$ | Poisson | $\log$ | $e^\eta$ |

---

## References

- F. Bach. *Learning Theory from First Principles*. MIT Press, 2024.
- P. L. Bartlett, M. I. Jordan, J. D. McAuliffe. Convexity, classification, and risk bounds. *JASA*, 101(473), 2006.
- L. Devroye, L. Györfi, G. Lugosi. *A Probabilistic Theory of Pattern Recognition*. Springer, 1996.
- J. Friedman, T. Hastie, R. Tibshirani. *The Elements of Statistical Learning*, 2nd ed. Springer, 2009.
- J. Mourtada. Exact minimax risk for linear least squares. *Ann. Statist.*, 50(4), 2022.
- V. N. Vapnik, A. Ya. Chervonenkis. On uniform convergence. *Theory Probab. Appl.*, 16(2), 1971.
