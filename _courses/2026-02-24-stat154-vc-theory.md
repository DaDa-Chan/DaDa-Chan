---
title: "Vapnik–Chervonenkis Theory"
date: 2026-02-24
category: "STAT 154"
excerpt: "VC dimension and shattering numbers as a combinatorial measure of model complexity. Symmetrisation, Rademacher complexity, Sauer–Shelah lemma, and the final ERM bound. Based on Hugo Chardon's STAT 154/254 lecture notes."
---

<div style="display:none">
$$
\newcommand{\FF}{\mathcal{F}}
\newcommand{\XX}{\mathcal{X}}
\newcommand{\Lhat}{\widehat{L}}
\newcommand{\fhat}{\widehat{f}}
\newcommand{\1}{\mathbf{1}}
\newcommand{\E}{\mathbb{E}}
\newcommand{\PP}{\mathbb{P}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\VC}{\mathrm{VC}}
$$
</div>

---

## 1. Starting Point: Bound for ERM in Binary Classification

**Setting.** Let $\XX$ be an arbitrary input space. A *model* $\FF$ is a collection of classifiers $f:\XX\to\{0,1\}$.

**ERM predictor.** Given i.i.d.\ data $(X_1,Y_1),\dots,(X_n,Y_n)$, the ERM predictor is

$$\fhat_n \;\in\; \operatorname*{arg\,min}_{f\in\FF}\; \frac{1}{n}\sum_{i=1}^{n}\1\{f(X_i)\neq Y_i\} \;=:\; \Lhat_n(f).$$

The population risk is $L(f) = \PP\bigl(f(X)\neq Y\bigr)$.

**Recall (1).** If $f^*\in\operatorname*{arg\,min}_{f\in\FF} L(f)$ is the best classifier in $\FF$, then

$$\begin{equation}\label{eq:excess-risk}
L(\fhat_n) - L(f^*) \;\leq\; 2\,\sup_{f\in\FF}\bigl|\Lhat_n(f)-L(f)\bigr|.
\end{equation}$$

Our goal is to bound the supremum, i.e.\ to prove *uniform convergence* of $\Lhat_n$ to $L$ over $\FF$.

**Recall (2).** When $\FF$ is *finite*, a simple union bound gives

$$\sup_{f\in\FF}\bigl|\Lhat_n(f)-L(f)\bigr| \;\lesssim\; \sqrt{\frac{\log|\FF|}{n}}.$$

---

## 2. VC Theory: Dealing with Infinite $\FF$

VC theory handles the supremum when $\FF$ is *not* finite (e.g.\ $\XX=\R^d$ and $\FF=\{\text{linear functions}\}$).

The key idea: measure the complexity of $\FF$ by the number of *distinct classification patterns* it can produce.

### 2.1 Shattering Numbers and VC Dimension

**Classification pattern.** Given a sample $x_1,\dots,x_n\in\XX$ and classifier $f\in\FF$, we obtain a *classification pattern* $(f(x_1),\dots,f(x_n))\in\{0,1\}^n$.

<div class="note-box note-def">
<span class="note-label">Definition — Shattering number</span>

For every $k\geq 1$, the $k$-th **shattering number** of $\FF$ is

$$S_k(\FF) \;=\; \sup_{x_1,\dots,x_k}\;\bigl|\bigl\{\bigl(f(x_1),\dots,f(x_k)\bigr) :\; f\in\FF\bigr\}\bigr|.$$

This is the maximum number of classification patterns achievable by *some* $f\in\FF$, for *some* sample $x_1,\dots,x_k$. It is a **worst-case** quantity.
</div>

Since all $f\in\FF$ are $\{0,1\}$-valued, $S_k(\FF)\leq 2^k$.

**Intuition.** If $S_k(\FF)=2^k$, then $\FF$ can *memorise* any training sample of size $k$ (every possible labelling is realised), so there is high risk of overfitting. When $S_k(\FF)=2^k$, we say **$\FF$ shatters** $x_1,\dots,x_k$.

<div class="note-box note-def">
<span class="note-label">Definition — VC dimension</span>

The **VC dimension** of $\FF$ is

$$\VC(\FF) \;:=\; \sup\bigl\{k\geq 1 :\; S_k(\FF)=2^k\bigr\} \;=\; \text{max \# of points that can be shattered by } \FF.$$
</div>

### 2.2 Interpretation / Common Mistakes

**$\VC(\FF)\geq d$.**
By definition, $S_d(\FF)=2^d$. Therefore, **to show that $\VC(\FF)\geq d$, it is enough to exhibit a single combination $x_1,\dots,x_d$ that $\FF$ can shatter.**

**$\VC(\FF)< d$.**
It suffices to prove that $S_d(\FF) < 2^d$, i.e.\ **no combination of $d$ points can be shattered by $\FF$.**

---

## 3. Bounding the Supremum via Empirical Processes

More precisely, we bound

$$\E\Bigl[\sup_{f\in\FF}\bigl|\Lhat_n(f)-L(f)\bigr|\Bigr] \;=\; \E\Bigl[\sup_{f\in\FF}\;\Bigl|\frac{1}{n}\sum_{i=1}^{n} \bigl\{\1\{f(X_i)\neq Y_i\}-\PP(f(X)\neq Y)\bigr\}\Bigr|\Bigr].$$

This is the *uniform deviation* of the empirical frequency from its true probability.

### 3.1 Symmetrisation for Empirical Processes

Write $Z_{i,f} = \1\{f(X_i)\neq Y_i\}$ and $Z = \sup_{f\in\FF}\bigl|\frac{1}{n}\sum_{i=1}^{n} Z_{i,f} - \E[Z_{1,f}]\bigr|$.

<div class="note-box note-thm">
<span class="note-label">Lemma — Symmetrisation</span>

$$\E[Z] \;\leq\; 2\,\E\Bigl[\sup_{f\in\FF}\;\Bigl| \frac{1}{n}\sum_{i=1}^{n}\varepsilon_i\, Z_{i,f}\Bigr|\Bigr],$$

where $\varepsilon_1,\dots,\varepsilon_n$ are i.i.d.\ **Rademacher** r.v.'s ($\PP(\varepsilon_i=1)=\PP(\varepsilon_i=-1)=\tfrac{1}{2}$), independent of all $Z_{i,f}$'s.
</div>

<div class="note-box note-rem">
<span class="note-label">Remark</span>

This is called "symmetrisation" because $\varepsilon_i Z_{i,f}$ is a *symmetric* random variable (same distribution as its opposite), so $\sum_{i=1}^{n}\varepsilon_i Z_{i,f}$ is also symmetric.
</div>

<div class="note-proof">
<em>Proof.</em> Let $\{Z'_{1,f},\dots,Z'_{n,f}:\; f\in\FF\}$ be independent copies of $\{Z_{1,f},\dots,Z_{n,f}:\; f\in\FF\}$. By identical distribution and independence, $\E[Z_{i,f}] = \E[Z'_{i,f}]$, which equals $\E[Z'_{i,f}\mid Z_{1,f},\dots,Z_{n,f}]$ by independence. Thus for all $f\in\FF$,

$$\begin{align}
\Bigl|\sum_{i=1}^{n} Z_{i,f} - \E[Z_{i,f}]\Bigr| &= \Bigl|\sum_{i=1}^{n}\E\bigl[Z_{i,f}-Z'_{i,f}\mid Z_{i,f}\bigr]\Bigr| \\
&= \Bigl|\E\Bigl[\sum_{i=1}^{n}\varepsilon_i(Z_{i,f}-Z'_{i,f})\;\Big|\; Z_{1,f},\dots,Z_{n,f}\Bigr]\Bigr|.
\end{align}$$

By Jensen's inequality and the triangle inequality:

$$\E\Bigl[\sup_{f\in\FF}\;\Bigl|\cdots\Bigr|\Bigr] \leq \E\Bigl[\sup_{f\in\FF}\;\Bigl|\sum_i\varepsilon_i Z_{i,f}\Bigr| + \sup_{f\in\FF}\;\Bigl|\sum_i\varepsilon_i Z'_{i,f}\Bigr|\Bigr] = 2\,\E\Bigl[\sup_{f\in\FF}\;\Bigl|\sum_i\varepsilon_i Z_{i,f}\Bigr|\Bigr].$$

Dividing by $n$ concludes the proof. $\square$
</div>

---

## 4. Intermediate Summary and Loss Patterns

Combining the excess risk bound $\eqref{eq:excess-risk}$ with the symmetrisation lemma:

$$\begin{equation}\label{eq:star}
\E\bigl[L(\fhat_n)-L^*\bigr] \;\leq\; 4\,\E\Bigl[\sup_{f\in\FF}\;\Bigl| \frac{1}{n}\sum_{i=1}^{n}\varepsilon_i\,\1\{f(X_i)\neq Y_i\}\Bigr|\Bigr].
\end{equation}$$

**Loss patterns equal classification patterns.** For two classifiers $f,g\in\FF$ and any $i$:

$$\1\{f(X_i)\neq Y_i\} \neq \1\{g(X_i)\neq Y_i\} \quad\Longleftrightarrow\quad f(X_i)\neq g(X_i).$$

Therefore **the number of possible loss patterns equals the number of classification patterns**, regardless of the labels $Y_1,\dots,Y_n$. This makes it clear **why shattering numbers (hence VC dimension) control the generalisation ability of the ERM predictor**.

---

## 5. Rademacher Complexity

<div class="note-box note-def">
<span class="note-label">Definition — Rademacher complexity</span>

Given a subset $\mathcal{B}\subseteq\R^n$, the **Rademacher complexity** of $\mathcal{B}$ is

$$R_n(\mathcal{B}) \;:=\; \E_\varepsilon\Bigl[\sup_{\beta\in\mathcal{B}}\; \Bigl|\frac{1}{n}\sum_{i=1}^{n}\varepsilon_i\,\beta_i\Bigr|\Bigr], \qquad \varepsilon_1,\dots,\varepsilon_n \text{ i.i.d.\ Rademacher}.$$
</div>

**Geometric intuition.** Given vectors $\beta\in\mathcal{B}$, randomly flip the signs of their coordinates, sum them, take the supremum over $\mathcal{B}$, and average over the sign flips.

### 5.1 Connecting to the ERM Bound

Let $\ell\circ\FF(\mathcal{D}_n) = \bigl\{\bigl(\ell(f(X_1),Y_1),\dots,\ell(f(X_n),Y_n)\bigr):\;f\in\FF\bigr\}$ denote the loss pattern set. The symmetrisation inequality states:

$$\E\Bigl[\sup_{f\in\FF}\bigl|\Lhat_n(f)-L(f)\bigr|\Bigr] \;\leq\; 2\,\E\bigl[R_n\bigl(\ell\circ\FF(\mathcal{D}_n)\bigr)\bigr].$$

### 5.2 Sub-Gaussian Bound on Rademacher Complexity

<div class="note-box note-thm">
<span class="note-label">Lemma — Max of sub-Gaussians</span>

Let $\mathcal{Z}_1,\dots,\mathcal{Z}_M$ be a collection of $K$-sub-Gaussian random variables. Then

$$\E\Bigl[\max_{1\leq i\leq M}|\mathcal{Z}_i|\Bigr] \;\leq\; \sqrt{2K^2\log M}.$$
</div>

In our case $\mathcal{B}$ consists of vectors with entries in $\{0,1\}$, so all $\beta\in\mathcal{B}$ satisfy $\|\beta\|^2\leq n$. Therefore

$$R_n\bigl(\ell\circ\FF(\mathcal{D}_n)\bigr) \;\leq\; \sqrt{\frac{2\log\bigl(2\,|\ell\circ\FF(\mathcal{D}_n)|\bigr)}{n}}.$$

### 5.3 Counting Loss Patterns via Shattering Numbers

Recall that the number of loss patterns equals the number of classification patterns: $|\ell\circ\FF(\mathcal{D}_n)| = |\FF(\mathcal{D}_n)|$, which is at most $S_n(\FF)$. In the end:

$$\begin{equation}\label{eq:shattering-bound}
\E\Bigl[\sup_{f\in\FF}\bigl|\Lhat_n(f)-L(f)\bigr|\Bigr] \;\leq\; 2\sqrt{\frac{2\log\bigl(2\,S_n(\FF)\bigr)}{n}}.
\end{equation}$$

---

## 6. The Sauer–Shelah Lemma

Recall $S_n(\FF)\leq 2^n$. The VC dimension is the *last* order at which this is an equality. Assume $\VC(\FF)=d<+\infty$.

**Case $n\leq d$.** $S_n(\FF)=2^n$, giving an uninformative bound.

**Case $n>d$.** When $n\geq 2d$:

<div class="note-box note-thm">
<span class="note-label">Sauer–Shelah Lemma</span>

If $n\geq 2d$, then

$$S_n(\FF) \;\leq\; \Bigl(\frac{en}{d}\Bigr)^d,$$

so $S_n(\FF)$ becomes **polynomial** in $n$ (instead of exponential). In particular, $\log\bigl(S_n(\FF)\bigr) \leq d\,\log\!\bigl(\frac{en}{d}\bigr)$.
</div>

### Final Bound for ERM

Combining $\eqref{eq:shattering-bound}$ with the Sauer–Shelah lemma:

<div class="note-box note-thm">
<span class="note-label">Final Bound for ERM</span>

$$\E\bigl[L(\fhat_n) - L^*\bigr] \;\lesssim\; \sqrt{\frac{\VC(\FF)\,\log\!\bigl(\tfrac{n}{\VC(\FF)}\bigr)}{n}}.$$
</div>

**Summary.** By looking at all the classification patterns $\FF$ can produce on a dataset of size $n$, we obtain $S_n(\FF)$—a *combinatorial* measure of the richness of $\FF$. Once $n$ exceeds the VC dimension $d$, $S_n(\FF)$ grows only polynomially in $n$ (rather than exponentially), and the ERM bound becomes meaningful.

