$$
\DeclareMathOperator{\Cov}{Cov}
\DeclareMathOperator{\Corr}{Corr}
\DeclareMathOperator{\Var}{Var}
\DeclareMathOperator{\prob}{\mathbb{P}}
\DeclareMathOperator{\qprob}{\mathbb{Q}}
\DeclareMathOperator{\E}{\mathbb{E}}
\newcommand{\set}[1]{\left\{#1\right\}}
\newcommand{\pa}[1]{\left(#1\right)}
\newcommand{\ang}[1]{\left<#1\right>}
\newcommand{\bra}[1]{\left[#1\right]}
\newcommand{\abs}[1]{\left|#1\right|}
\newcommand{\norm}[1]{\left\|#1\right\|}
$$

# Three unit rate model

## Model.

Neurons $\set{E_1,E_2,E_I}$ evolve with dynamics described by

$$
	\tau_i \dot{r}_i(t) = - r_i(t) + \bra{\mu_i + \sum_{j\in\set{E_1,E_2,I}}W_{ij}r_j(t)  + \sqrt{2\sigma_i}\xi_i(t)}_+, \qquad \text{for }i \in \set{E_1,E_2,I}
$$

^0898d9

where: 

- $\xi_i(t)$ are mutually independent white noise processes,
- $\tau_i$ are time constants
- $\sigma_i$ is the noise magnitude 
- $\mu_i$ is stimulus input to $i$ 
- $\bra{\cdot}_+ : \mathbb{R} \to \mathbb{R}$ is defined via $\bra{\cdot}\equiv \max\set{\cdot, 0}$ 

We may drop the $\bra{\cdot}_+$ under the assumption that the threshold is never sampled for sufficiently strong inputs. In vector form, these equations are 


$$\begin{aligned}
	\underbrace{\begin{bmatrix}
		\tau_{E_1} & 0 & 0 \\
		0 & \tau_{E_2} & 0 \\
		0 & 0 & \tau_I
	\end{bmatrix}}_{\boldsymbol{\tau}}	
	\begin{bmatrix}
		\dot{r}_{E_1} \\
		\dot{r}_{E_2} \\
		\dot{r}_{I}
	\end{bmatrix}
	=
	- \underbrace{\begin{bmatrix}
		r_{E_1} \\
		r_{E_2} \\
		r_{I}
		\end{bmatrix}}_{\mathbf{r}}
	+ 
		\underbrace{\begin{bmatrix}
			\mu_{E_1} \\
			\mu_{E_2} \\
			\mu_I
		\end{bmatrix}}_{\pmb{\mu}}
		& + 
		\underbrace{\begin{bmatrix}
			W_{E_1E_1} & W_{E_1E_2} & -W_{E_1I} \\
			W_{E_2E_1} & W_{E_2E_2} & - W_{E_2I}  \\
			W_{IE_1} & W_{IE_2} & -W_{II}
		\end{bmatrix}}_{\mathbf{W}}
		\begin{bmatrix}
			r_{E_1} \\
			r_{E_2}\\
			r_I
		\end{bmatrix} \\
		& \qquad +  
		\sqrt{2} 
		\underbrace{\begin{bmatrix}
			\sqrt{\sigma_{E_1}} & 0 & 0 \\
			0 & \sqrt{\sigma_{E_2}} &  0 \\
			0 & 0 & \sqrt{\sigma_{I}}
		\end{bmatrix}}_{\mathbf{S}}
		\underbrace{\begin{bmatrix}
			\xi_{E_1} \\
			\xi_{E_2}\\
			\xi_I
		\end{bmatrix}}_{\pmb{\xi}}
\end{aligned}$$

which can be written compactly as 
$$
\pmb{\tau} \mathbf{\dot{r}} = \pmb{\mu} + \pa{-\mathbf{1} + \mathbf{W}}\mathbf{r} + \sqrt{2}\mathbf{S} \pmb{\xi}.
$$
In the absence of noise, given $\mathbf{\bar{r}}$ we can solve for $\pmb{\mu}$ at the steady state (i.e., $\mathbf{\dot{r}}=0$)
$$
\pmb{\mu} = \pa{\mathbf{1} - \mathbf{W}}\mathbf{\bar{r}}.
$$
Let $\delta \mathbf{r} = \mathbf{r}-\mathbf{\bar{r}}$ denote the rates' fluctuations about their steady states. Then for sufficiently small noise, $\delta \mathbf{r}$ is an Ornstein-Uhlenbeck process given as the solution to the linear SDE 
$$
\frac{d}{dt}\delta \mathbf{r} = -\underbrace{\pmb{\tau}^{-1}\pa{\mathbf{1}-\mathbf{W}}}_{\mathbf{M}}\delta \mathbf{r} + \underbrace{\sqrt{2} \pmb{\tau}^{-1} \mathbf{S}}_{\mathbf{D}} \pmb{\xi} = -\mathbf{M} \delta \mathbf{r} + \mathbf{D} \pmb{\xi}.
$$

^e90dac

