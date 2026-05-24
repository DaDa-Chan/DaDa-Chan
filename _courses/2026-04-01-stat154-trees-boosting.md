---
title: "Decision Trees, Bagging, Random Forests, and Boosting"
date: 2026-04-01
category: "STAT 154"
excerpt: "Decision trees (regression and classification), cost-complexity pruning, bagging and OOB error, random forests, gradient boosting, AdaBoost, loss functions for boosting, XGBoost, and a bias–variance perspective on ensembles. Based on Hugo Chardon's STAT 154/254 lecture notes."
---

<div style="display:none">
$$
\newcommand{\R}{\mathbb{R}}
\newcommand{\norm}[1]{\| #1 \|}
\newcommand{\ind}{\mathbf{1}}
\newcommand{\argmin}{\operatorname*{argmin}}
\newcommand{\argmax}{\operatorname*{argmax}}
\newcommand{\sign}{\operatorname{sign}}
$$
</div>

*Typeset from STAT 154/254 lecture slides (Instructor: Hugo Chardon). Main references: ISL Chap. 8; ESL Chap. 9–10, 15–16; Breiman (2001).*

**Notation.** Training data $(x_1,y_1),\dots,(x_n,y_n)$ with $x_i\in\R^d$ and $y_i\in\R$ (regression) or $y_i\in\{1,\dots,K\}$ (classification). $x_{i,j}$ denotes the $j$-th feature of observation $i$. A decision tree partitions the feature space into $J$ disjoint regions $R_1,\dots,R_J$; $\hat{y}_{R_j}$ is the prediction in region $R_j$. For classification, $\hat{p}_{jk}$ is the proportion of training observations in region $R_j$ from class $k$.

---

## 1. Decision Trees

### 1.1 The Idea

Tree-based methods recursively partition the feature space into axis-aligned rectangular regions and fit a constant within each region. At a given internal node, a split "$x_j < s$" sends observations left; those with $x_j\geq s$ go right. Terminal nodes (leaves) contain the predictions.

**Strengths:** simple, interpretable, handle qualitative predictors without dummy variables.  
**Weakness:** lower predictive accuracy and high variance — a small perturbation of the data can change the entire tree structure.

### 1.2 Regression Trees

#### Objective

We seek a partition $R_1,\dots,R_J$ of $\R^d$ minimising the residual sum of squares (RSS):

$$\text{RSS} \;=\; \sum_{j=1}^{J}\sum_{i:\,x_i\in R_j}(y_i - \hat{y}_{R_j})^2, \qquad \hat{y}_{R_j} = \frac{1}{|R_j|}\sum_{i:\,x_i\in R_j}y_i.$$

Finding the globally optimal partition is intractable; instead we use a greedy, top-down procedure.

#### Recursive Binary Splitting

For a candidate feature $j$ and threshold $s$, define $R_1(j,s) = \{x: x_j < s\}$, $R_2(j,s) = \{x: x_j \geq s\}$. At each step, choose $j$ and $s$ to minimise the combined RSS of the two children:

$$\min_{j,\,s}\;\Biggl[\sum_{i:\,x_i\in R_1(j,s)}\!\!(y_i - \hat{y}_{R_1})^2 \;+\; \sum_{i:\,x_i\in R_2(j,s)}\!\!(y_i - \hat{y}_{R_2})^2\Biggr].$$

For each fixed $j$, the minimisation over $s$ only requires scanning the $n$ observed values of $x_{i,j}$, so the cost per split is $O(nd)$.

<div class="note-algo">
<div class="algo-title">Algorithm 1: Recursive Binary Splitting (Regression)</div>

**Input:** Training data $(x_i,y_i)_{i=1}^n$; minimum leaf size $n_{\min}$  
**Output:** A regression tree $T$

1. Initialise the partition $\mathcal{P} \gets \{R_0\}$ where $R_0 = \R^d$ (all data in one region)
2. **While** there exists a leaf $R\in\mathcal{P}$ with $|R| > n_{\min}$:
   - For each such leaf $R$, compute the best split $(j_R^{\ast}, s_R^{\ast})$ by minimising RSS over all features $j$ and thresholds $s$
   - Select the leaf $R^{\ast}$ whose best split yields the largest decrease in RSS
   - Replace $R^{\ast}$ in $\mathcal{P}$ by $R^{\ast}_{\text{left}} = \{x\in R^{\ast}: x_{j^{\ast}} < s^{\ast}\}$ and $R^{\ast}_{\text{right}} = \{x\in R^{\ast}: x_{j^{\ast}}\geq s^{\ast}\}$
3. Assign prediction $\hat{y}_{R} = \frac{1}{|R|}\sum_{i:\,x_i\in R}y_i$ in each leaf $R\in\mathcal{P}$
</div>

#### Cost-Complexity Pruning

Growing a very deep tree $T_0$ leads to overfitting. For a tuning parameter $\alpha\geq 0$, **cost-complexity pruning** finds the subtree $T\subseteq T_0$ minimising

$$C_\alpha(T) \;=\; \sum_{j=1}^{|T|}\sum_{i:\,x_i\in R_j}(y_i - \hat{y}_{R_j})^2 \;+\; \alpha\,|T|,$$

where $|T|$ is the number of terminal nodes. As $\alpha$ increases, branches are pruned in a nested fashion.

<div class="note-algo">
<div class="algo-title">Algorithm 2: Building a Regression Tree with Pruning</div>

**Input:** Training data $(x_i,y_i)_{i=1}^n$; number of CV folds $K$  
**Output:** A pruned regression tree $T^{\ast}$

1. Grow a large tree $T_0$ by recursive binary splitting (Algorithm 1)
2. For a grid of $\alpha$ values, compute the best subtree $T_\alpha\subseteq T_0$ minimising $C_\alpha(T)$
3. Use $K$-fold cross-validation to select $\alpha^{\ast}$:
   - For each fold $k=1,\dots,K$: grow $T_0^{(-k)}$ on all folds except $k$; prune to $T_\alpha^{(-k)}$ for each $\alpha$; evaluate MSE on the held-out $k$-th fold
   - Set $\alpha^{\ast} = \argmin_\alpha\;\text{average CV error}(\alpha)$
4. **Return** $T^{\ast} = T_{\alpha^{\ast}}$
</div>

### 1.3 Classification Trees

For classification, each leaf predicts the majority class. RSS is not appropriate; instead we use measures of **node impurity**.

<div class="note-box note-def">
<span class="note-label">Definition — Split criteria for classification</span>

Let $\hat{p}_{jk}$ be the proportion of training observations in region $R_j$ from class $k$.

- **Gini index:** $G_j = \sum_{k=1}^{K}\hat{p}_{jk}(1-\hat{p}_{jk})$.
- **Entropy:** $D_j = -\sum_{k=1}^{K}\hat{p}_{jk}\log\hat{p}_{jk}$.

Both are small when a node is *pure* and maximal when all classes are equally represented. Numerically they are quite similar.
</div>

The splitting procedure is identical to the regression case, replacing RSS by the weighted Gini or entropy of the two children.

<div class="note-box note-rem">
<span class="note-label">Remark — Advantages and disadvantages of single trees</span>

**Advantages:** easy to explain and interpret; natural graphical representation; handle qualitative predictors without encoding.

**Disadvantages:** lower predictive accuracy than ensemble methods; high variance (non-robust to small data perturbations). These weaknesses motivate the ensemble methods below.
</div>

---

## 2. Bagging (Bootstrap Aggregation)

### 2.1 Motivation: Variance Reduction by Averaging

An **ensemble method** combines many simple *base learners* into a single stronger model. Decision trees are natural base learners: they have low bias (when grown deep) but high variance. Averaging many high-variance, low-bias estimators reduces variance without increasing bias.

### 2.2 The Bagging Procedure

<div class="note-algo">
<div class="algo-title">Algorithm 3: Bagging for Regression Trees</div>

**Input:** Training data $(x_i,y_i)_{i=1}^n$; number of bootstrap replicates $B$  
**Output:** Bagged predictor $\hat{f}_{\mathrm{bag}}$

- For $b=1,\dots,B$:
  - Draw a bootstrap sample $(x_1^{\ast b},y_1^{\ast b}),\dots,(x_n^{\ast b},y_n^{\ast b})$ by sampling $n$ observations **with replacement**
  - Grow a deep, unpruned tree $\hat{f}^{\ast b}$ on this bootstrap sample
- **Regression:** $\hat{f}_{\mathrm{bag}}(x) = \frac{1}{B}\sum_{b=1}^{B}\hat{f}^{\ast b}(x)$
- **Classification:** $\hat{G}_{\mathrm{bag}}(x) = \text{majority vote of }\hat{G}^{\ast 1}(x),\dots,\hat{G}^{\ast B}(x)$
</div>

### 2.3 Out-of-Bag (OOB) Error

Each bootstrap sample omits roughly a fraction $1-1/e\approx 36.8\%$ of the original observations. For each training point $x_i$, define its **OOB prediction** as the average (or majority vote) over only those trees whose bootstrap sample did not include $x_i$. The OOB error is a nearly unbiased estimate of the test error, and serves as a convenient alternative to cross-validation.

---

## 3. Random Forests

### 3.1 Decorrelating the Trees

Bagging reduces variance by averaging, but if a single strong predictor dominates, all $B$ trees tend to look similar — they are *correlated*. Averaging correlated quantities yields less variance reduction. Random forests address this by **decorrelating** the trees.

<div class="note-algo">
<div class="algo-title">Algorithm 4: Random Forest</div>

**Input:** Training data; number of trees $B$; number of candidate features per split $m$  
**Output:** Random forest predictor $\hat{f}_{\mathrm{RF}}$

- For $b=1,\dots,B$:
  - Draw a bootstrap sample of size $n$
  - Grow a deep tree $\hat{f}^{\ast b}$, but at **each split**:
    - Select $m$ features uniformly at random from $\{1,\dots,d\}$
    - Choose the best split *among these $m$ features only*
- Aggregate as in bagging (average for regression, majority vote for classification)
</div>

<div class="note-box note-rem">
<span class="note-label">Remark — Default choices for $m$</span>

Common defaults: $m = \lfloor d/3\rfloor$ for regression, $m=\lfloor\sqrt{d}\rfloor$ for classification. When $m=d$, random forests reduce to bagging.
</div>

### 3.2 Hyperparameters

The following can be tuned via cross-validation or OOB error: number of trees $B$; number of candidate features $m$; maximum tree depth; minimum leaf size; split criterion (MSE for regression; Gini or entropy for classification).

---

## 4. Boosting

Unlike bagging, which trains trees independently in parallel, **boosting** grows trees *sequentially*: each new tree is fitted to the errors of the current ensemble.

### 4.1 Boosting for Regression Trees

<div class="note-algo">
<div class="algo-title">Algorithm 5: Gradient Boosting for Regression (Squared Loss)</div>

**Input:** Training data $(x_i,y_i)_{i=1}^n$; iterations $B$; shrinkage $\lambda>0$; tree depth $J$  
**Output:** Boosted predictor $\hat{f}$

1. Initialise $\hat{f}(x)\gets 0$ and residuals $r_i\gets y_i$ for all $i$
2. For $b=1,\dots,B$:
   - Fit a regression tree $\hat{f}^{b}$ with $J$ terminal nodes to the data $(x_i, r_i)_{i=1}^n$
   - Update the model: $\hat{f}(x) \gets \hat{f}(x) + \lambda\,\hat{f}^{b}(x)$
   - Update residuals: $r_i \gets r_i - \lambda\,\hat{f}^{b}(x_i)$
3. **Return** $\hat{f}(x) = \sum_{b=1}^{B}\lambda\,\hat{f}^{b}(x)$
</div>

<div class="note-box note-rem">
<span class="note-label">Remark — Hyperparameters</span>

- **Number of trees $B$:** unlike bagging, boosting *can* overfit if $B$ is too large; select via cross-validation.
- **Shrinkage $\lambda$:** a small learning rate (typically $0.01$ or $0.001$) slows down learning; smaller $\lambda$ requires larger $B$.
- **Tree depth $J$:** controls the interaction order. Often $J=1$ (a *stump*) works well, giving a purely additive model.
</div>

### 4.2 Boosting for Classification: AdaBoost

For binary classification with $y_i\in\{-1,+1\}$, the ensemble takes the form $G(x) = \sign\bigl(\sum_{t=1}^{T}\alpha_t G_t(x)\bigr)$. AdaBoost fits this by iteratively reweighting the training data to focus on misclassified observations.

<div class="note-algo">
<div class="algo-title">Algorithm 6: AdaBoost</div>

**Input:** Training data $(x_i,y_i)_{i=1}^n$ with $y_i\in\{-1,+1\}$; iterations $T$; base learner class $\mathcal{F}$  
**Output:** Ensemble $G(x) = \sign\bigl(\sum_{t=1}^{T}\alpha_t G_t(x)\bigr)$

1. Initialise weights: $w_i \gets 1/n$ for $i=1,\dots,n$
2. For $t=1,\dots,T$:
   - Fit base learner to weighted data: $G_t = \argmin_{G\in\mathcal{F}}\;\sum_{i=1}^{n}w_i\,\ind\{G(x_i)\ne y_i\}$
   - Compute weighted error: $\varepsilon_t = \frac{\sum_{i=1}^{n}w_i\,\ind\{G_t(x_i)\ne y_i\}}{\sum_{i=1}^{n}w_i}$
   - Compute learner weight: $\alpha_t = \tfrac{1}{2}\ln\frac{1-\varepsilon_t}{\varepsilon_t}$
   - Update sample weights: $w_i \gets w_i\,\exp\bigl(-\alpha_t\,y_i\,G_t(x_i)\bigr)$ for all $i$
   - Normalise: $w_i \gets w_i\big/\sum_{j}w_j$
   </div>

<div class="note-box note-rem">
<span class="note-label">Remark — AdaBoost does not easily overfit</span>

In the absence of label noise, AdaBoost's test error often continues to decrease even after training error reaches zero. This is explained by *margin theory*: each additional iteration increases the margin of the ensemble, which controls generalisation. However, with noisy labels, mislabelled examples accumulate exponentially large weights.
</div>

### 4.3 Loss Functions for Boosting

All boosting variants can be understood as greedily minimising an empirical risk $\sum_i \ell(y_i\,f(x_i))$ with different surrogate loss functions $\ell$ that upper-bound the $0$-$1$ loss $\ind\{yf(x)<0\}$:

| Name | Loss $\ell(u)$ | Algorithm |
|------|----------------|-----------|
| Exponential | $e^{-u}$ | AdaBoost |
| Logistic | $\log_2(1+e^{-u})$ | LogitBoost |
| Quadratic | $(1-u)^2$ | GentleBoost |
| Hinge | $(1-u)_+$ | SVM |

All four are convex upper bounds on the $0$-$1$ loss. The logistic loss is more robust to label noise than the exponential loss.

### 4.4 Gradient Boosting: the General Framework

The key insight is that boosting performs **gradient descent in function space**: the vector of current predictions $\mathbf{f}=(f_1,\dots,f_n)\in\R^n$ is updated in the direction that most reduces $\sum_i\ell(f_i,y_i)$. Because we must stay within the class $\mathcal{F}$, we fit the next base learner to the *negative gradient*.

<div class="note-algo">
<div class="algo-title">Algorithm 7: Gradient Boosting (General Loss)</div>

**Input:** Training data; iterations $T$; differentiable loss $\ell(\cdot,\cdot)$; base learner class $\mathcal{F}$  
**Output:** Ensemble $\hat{f}(x) = \sum_{t=1}^{T}\alpha_t G_t(x)$

1. Initialise predictions: $f_i\gets 0$ for all $i$
2. For $t=1,\dots,T$:
   - Compute pseudo-residuals (negative gradient): $g_i = -\frac{\partial\,\ell(f_i,\,y_i)}{\partial f_i}$, $i=1,\dots,n$
   - Fit base learner to pseudo-residuals: $G_t = \argmin_{G\in\mathcal{F}}\;\sum_{i=1}^{n}\bigl(G(x_i) - g_i\bigr)^2$
   - Line search for step size: $\alpha_t = \argmin_{\alpha>0}\;\sum_{i=1}^{n}\ell\bigl(f_i + \alpha\,G_t(x_i),\;y_i\bigr)$
   - Update predictions: $f_i\gets f_i + \alpha_t\,G_t(x_i)$ for all $i$
   </div>

<div class="note-box note-rem">
<span class="note-label">Remark — Special cases</span>

With squared loss $\ell(f,y) = \frac{1}{2}(y-f)^2$, the pseudo-residuals become $g_i = y_i - f_i$ (actual residuals), and we recover Algorithm 5. With exponential loss $\ell(f,y) = e^{-yf}$, one recovers AdaBoost (Algorithm 6).
</div>

### 4.5 XGBoost (eXtreme Gradient Boosting)

XGBoost extends standard gradient boosting with explicit **regularisation**. At iteration $T$, it adds a tree $G_T(x) = \sum_{k=1}^{|T|}w_k\,\ind\{x\in\text{leaf }k\}$ by minimising

$$\widetilde{Q}_T \;=\; \sum_{i=1}^{n}\ell\bigl(f_i + G_T(x_i),\;y_i\bigr) \;+\; \gamma\,|T| \;+\; \frac{\lambda}{2}\sum_{k=1}^{|T|}w_k^2,$$

where $|T|$ is the number of leaves, $w_k$ is the prediction in leaf $k$, and $\gamma,\lambda\geq 0$ are regularisation parameters. The $\gamma|T|$ term acts like cost-complexity pruning ($\ell_0$-type), while $\frac{\lambda}{2}\sum w_k^2$ is $\ell_2$ shrinkage on leaf values.

XGBoost uses a **second-order Taylor expansion** of $\ell$ for more accurate optimisation, and includes practical enhancements: built-in handling of missing values, column subsampling (à la random forests), and GPU support.

---

## 5. Why Do Ensembles Work? A Bias–Variance Perspective

<div class="note-box note-prop">
<span class="note-label">Proposition — Bias–variance intuition for ensembles</span>

The two main ensemble strategies occupy opposite ends of the bias–variance spectrum:

1. **Bagging / Random Forests.** Base learners are deep trees (low bias, high variance). Averaging $B$ of them reduces variance while preserving the low bias.
2. **Boosting.** Base learners are shallow trees or stumps (high bias, low variance). Sequential addition of shrunken trees reduces bias, and with appropriate shrinkage, also controls variance.
</div>

<div class="note-box note-rem">
<span class="note-label">Remark — Practical comparison</span>

Bagging is better for small or noisy training sets (robust to label noise) and can train base learners in parallel. Boosting tends to achieve higher accuracy on clean data with complex decision boundaries, but is inherently sequential and sensitive to label noise. In practice: random forests are preferred for robustness and speed; XGBoost/LightGBM for maximum accuracy on tabular data.
</div>
