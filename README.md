# superquantile-opt

A semismooth-Newton-based proximal augmented Lagrangian (p-ALM) solver$`^{[1]}`$ for the sample average approximation (SAA) of the superquantile constrained problem
```math
\begin{equation}
\begin{array}{rl}
\displaystyle \underset{{x\in X\subseteq \mathbb{R}^n}}{\text{minimize}}\; & \displaystyle f(x) + \mathsf{CVaR}_{\tau_0}\bigl(G^{\,0}(x;\omega^{[m]})\bigr) \\[0.1in]
\mbox{subject to}
& \mathsf{CVaR}_{\tau_{\ell}}\bigl(G^{\,\ell}(x;\omega^{[m]})\bigr) \leq 0,\; \ell \in\{1,2,\ldots,L\},
\end{array}
\end{equation}
```
where
- $\mathsf{CVaR}_{\tau}\bigl(y\bigr)$ denotes the (empirical) superquantile of a (realization of a random) vector $y\in\mathbb{R}^m$ at confidence $\tau\in(0,1)$ with $`\mathsf{CVaR}_{\tau}\bigl(y\bigr) = k(\tau)^{-1}\mathsf{T}_{k(\tau)}(y)`$ where $k(\tau) := m\cdot(1-\tau)$ and where $`\mathsf{T}_k(\cdot)`$ denotes the top- $\\!\\!k$-sum operator$`^{[2]}`$
- $f$ is smooth and convex;
- for each $\ell\in\\{0,1,\ldots,L\\}$, $`G^\ell(x;\omega^{[m]}) := \bigl\{g^\ell(x;\omega^j)\bigr\}_{j=1}^m`$ is a collection of convex continuously differentiable scalar functions $`g^\ell(\,\cdot\,; \omega)`$ evaluated at $m$ SAA scenarios $`\{ \omega^j \}_{j=1}^m`$.
  

The repository contains two directories related to the paper https://arxiv.org/abs/2405.07965:
- `src/`: implementations of the p-ALM.
- `run/`: scripts for running the experiments.

---
$`^{[1]}`$ See https://arxiv.org/abs/2405.07965 for methodological details.<br>
$`^{[2]}`$ See https://arxiv.org/abs/2310.07224 and https://github.com/jacob-roth/top-k-sum.
