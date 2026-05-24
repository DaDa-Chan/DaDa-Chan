---
title: "The Kernel Trick"
date: 2026-03-10
category: "STAT 154"
excerpt: "From dot products to kernels. Feature maps, positive semi-definite kernels, reproducing kernel Hilbert spaces, the reproducing property, canonical feature map, and the Representer Theorem. Based on Hugo Chardon's STAT 154/254 lecture notes §SVM 14–21."
---

<div style="display:none">
$$
\newcommand{\R}{\mathbb{R}}
\newcommand{\Hcal}{\mathcal{H}}
\newcommand{\Xcal}{\mathcal{X}}
\newcommand{\ip}[2]{\langle #1,\, #2 \rangle}
\newcommand{\ipH}[2]{\langle #1,\, #2 \rangle_{\mathcal{H}}}
\newcommand{\norm}[1]{\| #1 \|}
\newcommand{\normH}[1]{\| #1 \|_{\mathcal{H}}}
$$
</div>

*Typeset from handwritten lecture notes (SVM §14–21).*

---

## 1. Motivation: From Dot Products to Kernels

Recall that in the SVM dual formulation, the optimal Lagrange multiplier vector $\mu^{\ast}$ is determined by solving a quadratic program whose data-dependence enters *solely* through the matrix of pairwise dot products

$$\widetilde{K} = \bigl(\ip{y_i x_i}{\, y_j x_j}\bigr)_{1\leq i,j\leq n}.$$

More precisely, defining $\tilde{x}_i = y_i x_i$ and $\mathbb{1}=(1,\dots,1)^\top\in\R^n$, the dual reads

$$\mu^{\ast} = \operatorname*{argmax}_{\mu\in\R_+^n} \Bigl\{\ip{\mathbb{1}}{\mu} - \tfrac{1}{2}\,\norm{\mu}_{\widetilde{K}}^2\Bigr\}.$$

This is a canonical quadratic program with a simple convex constraint (the non-negative orthant $\R_+^n$), efficiently solvable by libraries such as `cvxpy`.

<div class="note-box note-rem">
<span class="note-label">Remark — Dimensionality advantage</span>

In a high-dimensional setting where $n \ll d$, solving a problem in $\R^n$ (via the dual) rather than in $\R^d$ (via the primal) is far more efficient.
</div>

Two key observations motivate the kernel trick:

1. **Data dependence.** The dual depends on the data *only* through the Gram matrix of dot products $\widetilde{K} = \bigl(\ip{x_i}{x_j}\bigr)_{1\leq i,j\leq n}$ (after absorbing the labels $y_i$).
2. **Non-separability.** If the data are not linearly separable in $\R^d$, we can lift them into a higher-dimensional space where a linear separator may exist.

---

## 2. Feature Maps and Kernels

### 2.1 Lifting to Higher Dimensions

When data in $\R^d$ are not linearly separable, map them into a richer feature space via a *feature map*

$$\varphi\colon \R^d \longrightarrow \Hcal,$$

where $\Hcal$ is a (possibly infinite-dimensional) Hilbert space equipped with an inner product $\ipH{\cdot}{\cdot}$ and induced norm $\normH{u} = \sqrt{\ipH{u}{u}}$. In this lifted space, the data $\varphi(x_1),\dots,\varphi(x_n)$ may become linearly separable.

### 2.2 The Kernel Function

Since the dual depends on the data only through pairwise inner products, after lifting we need only the quantities $\ipH{\varphi(x_i)}{\varphi(x_j)}$. We do *not* need to compute $\varphi(x)$ explicitly—we only need the function

$$\boxed{K(x_i, x_j) \;=\; \ipH{\varphi(x_i)}{\varphi(x_j)}.}$$

This is the **kernel function** associated with the feature map $\varphi$.

### 2.3 Properties of Kernels

Because $K$ arises from an inner product, it must inherit:

1. **Symmetry.** $K(x, x') = K(x', x)$ for all $x,x'\in\R^d$.
2. **Positive semi-definiteness (PSD).** For every $n\in\mathbb{N}$ and every $x_1,\dots,x_n\in\R^d$, the Gram matrix $\bigl(K(x_i, x_j)\bigr)_{1\leq i,j\leq n}$ is PSD: for all $\alpha\in\R^n$,

$$\sum_{i,j=1}^{n} \alpha_i\,\alpha_j\,K(x_i,x_j) \;\geq\; 0.$$

---

## 3. Reproducing Kernel Hilbert Spaces (RKHS)

The notion that rigorously bridges kernels and feature spaces is that of a **Reproducing Kernel Hilbert Space**.

<div class="note-box note-def">
<span class="note-label">Definition — Reproducing Kernel</span>

Let $\Hcal$ be a Hilbert space of functions $f\colon\Xcal\to\R$, and let $K\colon\Xcal\times\Xcal\to\R$ be a kernel on $\Xcal$. We say $K$ is a *reproducing kernel* for $\Hcal$ if $K$ is symmetric and:

1. For all $x\in\Xcal$, the function $K(x,\cdot)\in\Hcal$.
2. For all $f\in\Hcal$ and $x\in\Xcal$:

$$f(x) = \ipH{K(x,\cdot)}{f}.$$

Condition 2 is called the **reproducing property**. It makes function evaluation a linear operation in $\Hcal$.
</div>

<div class="note-box note-rem">
<span class="note-label">Remark — Intuition</span>

The reproducing property says that applying a function $f$ to a point $x$ is the same as taking the inner product of $f$ with the "feature" $K(x,\cdot)$. If $f\colon\R^d\to\R$ is linear, then there exists a unique $\theta\in\R^d$ such that $f(x)=\ip{\theta}{x}$ for all $x$. The RKHS generalises this: when we want to allow nonlinear functions, we lift to a space $\Hcal$ where the functions again admit a linear (inner-product) representation — that is precisely the combination of feature maps and the reproducing property.
</div>

### 3.1 Correspondence Between Kernels and Hilbert Spaces

The following two propositions establish a one-to-one correspondence between symmetric PSD kernels and RKHS.

<div class="note-box note-prop">
<span class="note-label">Proposition 3.1 — RKHS ⟹ PSD kernel</span>

Let $\Hcal$ be an RKHS over $\Xcal$ with reproducing kernel $K$. Then $K$ is positive semi-definite.
</div>

<div class="note-box note-prop">
<span class="note-label">Proposition 3.2 — PSD kernel ⟹ unique RKHS (Moore–Aronszajn)</span>

Let $K$ be a positive semi-definite kernel on $\Xcal$. Then there exists a unique RKHS $\Hcal$ of functions $\Xcal\to\R$ for which $K$ is the reproducing kernel.
</div>

### 3.2 The Canonical Feature Map

Given a Hilbert space $\Hcal$ with reproducing kernel $K$, there is a natural choice of feature map:

$$\varphi(x) \;=\; K(x,\cdot) \;\in\;\Hcal.$$

This construction is **canonical**: two space/kernel pairs $(\Hcal_0,K_0)$ and $(\Hcal_1,K_1)$ with $\varphi_i(x)=K_i(x,\cdot)$ satisfy

$$K_0 = K_1 \quad\Longleftrightarrow\quad \varphi_0(x) = \varphi_1(x)\;\;\text{for all }x.$$

---

## 4. The Representer Theorem

### 4.1 Kernel SVM: the Primal in $\Hcal$

After lifting via a feature map $\varphi\colon\Xcal\to\Hcal$, the soft-SVM primal becomes

$$\min_{f\in\Hcal}\;\Biggl\{ \sum_{i=1}^{n} h\bigl(y_i\,\ipH{f}{\varphi(x_i)}\bigr) \;+\; \lambda\,\normH{f}^2 \Biggr\}, \tag{$\star$}$$

where $h(u)=(1-u)_+$ is the hinge loss. The dual problem involves only the kernel matrix $K_n = \bigl(K(x_i,x_j)\bigr)_{1\leq i,j\leq n}\in\R^{n\times n}$.

Problem $(\star)$ requires optimising over $\Hcal$, which may be infinite-dimensional. The representer theorem guarantees that the solution lives in a *finite-dimensional* subspace.

<div class="note-box note-thm">
<span class="note-label">Theorem — Representer Theorem</span>

Let $\varphi\colon\Xcal\to\Hcal$ be a feature map and $x_1,\dots,x_n\in\Xcal$ input vectors. Then the solution to $(\star)$ lies in $\operatorname{Span}\bigl(\{\varphi(x_i)\}_{1\leq i\leq n}\bigr)$, i.e.\ it is of the form

$$f^{\ast} = \sum_{i=1}^{n} \alpha_i\,\varphi(x_i)$$

for some coefficients $\alpha_1,\dots,\alpha_n\in\R$.
</div>

<div class="note-box note-rem">
<span class="note-label">Remark — Practical consequence</span>

Although $\Hcal$ can be infinite-dimensional, the representer theorem reduces the optimisation to finding $n$ real coefficients $\alpha\in\R^n$. The only data-dependent quantity needed is the $n\times n$ kernel matrix $K_n$.
</div>

### 4.2 Summary

Starting from a classification problem in $\R^d$ that may not be linearly separable:

1. Choose a PSD kernel $K$ (equivalently, choose an RKHS $\Hcal$ and the canonical feature map $\varphi(x)=K(x,\cdot)$).
2. Replace every inner product $\ip{x_i}{x_j}$ in the algorithm by $K(x_i,x_j) = \ipH{\varphi(x_i)}{\varphi(x_j)}$.
3. By the representer theorem, the optimiser is $f^{\ast}=\sum_i \alpha_i\,\varphi(x_i)$, so only the $n\times n$ kernel matrix $K_n$ is needed — $\varphi$ is never computed explicitly.

