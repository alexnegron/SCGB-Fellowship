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

# Long-time covariance calculation
In [[Three unit rate model|the note defining the rate model]], we wrote down a linear SDE characterizing the dynamics of the fluctuations $\delta \mathbf{r}.$ This SDE is given by: 

$$
\tau_i \dot{r}_i(t) = - r_i(t) + \bra{\mu_i + \sum_{j\in\set{E_1,E_2,I}}W_{ij}r_j(t)  + \sqrt{2\sigma_i}\xi_i(t)}_+, \qquad \text{for }i \in \set{E_1,E_2,I}
$$

The solution is obtained as

$$
	\delta \mathbf{r}(t) = \exp\left\{-\mathbf{M}t\right\}\delta\mathbf{r}(0) + \int_{-\infty}^{t}\exp\left\{-\mathbf{M}\pa{t-t'}\right\}\mathbf{D}\ \underbrace{d\mathbf{W}(t')}_{\boldsymbol{\xi}(t')}
$$

with mean $\ang{\delta \mathbf{r}(t)} = \exp\set{-\mathbf{M}t}\ang{\delta \mathbf{r}(0)}$ (expectation of integral is given by an [Ito integral property](marginnote3app://note/5CC92031-57BF-4688-9A49-DB3BCBB8BE43)).

## Long-time covariance function.

We now compute the $t\to \infty$ covariance function. Observe $\lim_{t \to \infty} \ang{\delta \mathbf{r}(t)} = 0$ (see the [[#^adf0fd|Remark]] below for the details about this calculation). As an operator, $\Cov(\cdot, \cdot)$ is shift-invariant in each argument, thus 

$$\Cov(\mathbf{r}(t), \mathbf{r}(s)) = \Cov(\mathbf{r}(t)-\mathbf{\bar{r}}, \mathbf{r}(s)-\mathbf{\bar{r}}) = \Cov(\delta \mathbf{r}(t), \delta\mathbf{r}(s)).$$

Hence in the following calculations, for notational simplicity, we will use $\mathbf{r}(t)$ in place of $\delta \mathbf{r}(t)$.

We make the notations

$$\Cov(\mathbf{r}(t), \mathbf{r}(s)) =: \ang{\mathbf{r}(t), \mathbf{r}(s)} 
	= 
	\E\bra{\pa{\mathbf{r}(t)-\E\bra{\mathbf{r}(t)}}\pa{\mathbf{r}(s)-\E\bra{\mathbf{r}(s)}}^\top} 
	= 
	\ang{\pa{\mathbf{r}(t) - \ang{\mathbf{r}(t)}}\pa{\mathbf{r}(s) - \ang{\mathbf{r}(s)}}^\top}.$$
	
For $s=t+\tau$ where $\tau > 0$, we have 

$$\begin{aligned}
	\ang{\pa{\mathbf{r}(t)-\ang{\mathbf{r}(t)}}\pa{\mathbf{r}(t+\tau)-\ang{\mathbf{r}(t+\tau)}}^\top}
	&=
	\ang{\mathbf{r}(t)\mathbf{r}^{\top}(t+\tau)-\ang{\mathbf{r}(t)}\mathbf{r}^{\top}(t+\tau) - \mathbf{r}(t)\ang{\mathbf{r}(t+\tau)}^{\top} + \ang{\mathbf{r}(t)}\ang{\mathbf{r}(t+\tau)}^{\top}}\\
	&=
	\ang{\mathbf{r}(t)\mathbf{r}^{\top}(t+\tau)} - 2\underbrace{\ang{\mathbf{r}(t)}\ang{\mathbf{r}^{\top}(t+\tau)}}_{\xrightarrow[t\to\infty]{}0} + \underbrace{\ang{\mathbf{r}(t)}\ang{\mathbf{r}(t+\tau)}^{\top}}_{\xrightarrow[t\to\infty]{}0},
\end{aligned}
$$

which implies that the long-time covariance function is the $t\to \infty$ limit of $\ang{\mathbf{r}(t)\mathbf{r}^{\top}(t+\tau)}$: 

$$
\begin{aligned}
\lim_{t\to \infty}
\ang{\mathbf{r}(t)\mathbf{r}^{\top}(t+\tau)} 
&= 
\lim_{t\to\infty}
\ang{
	\int_{-\infty}^{t}\int_{-\infty}^{t+\tau}
	\exp\set{-\mathbf{M}(t-t')}\mathbf{D} \pmb{\xi}(t')
	\pa{\exp\set{-\mathbf{M}\pa{t+\tau-t''}}\mathbf{D} \pmb{\xi}(t'')}^\top \ dt'' \ dt'
	}\\
&=
	\lim_{t\to \infty} 
	\int_{-\infty}^{t}\int_{-\infty}^{t+\tau}
	\ang{
		\exp\set{-\mathbf{M}(t-t')}\mathbf{D} \pmb{\xi}(t') 
		\pmb{\xi}^{\top}(t'')\mathbf{D}^{\top}\exp\set{-\mathbf{M}^{\top}(t+\tau-t'')}} \ dt'' \ dt'
		 \\
&=
	\lim_{t\to \infty} 
	\int_{-\infty}^{t}\int_{-\infty}^{t+\tau}
	\exp\set{-\mathbf{M}(t-t')}\mathbf{D} \ang{\pmb{\xi}(t') 
		\pmb{\xi}^{\top}(t'')}\mathbf{D}^{\top}\exp\set{-\mathbf{M}^{\top}(t+\tau-t'')} \ dt'' \ dt' \\
&= 
	\lim_{t\to\infty}
	\int_{-\infty}^{t}
	\exp\set{-\mathbf{M}(t-t')}\mathbf{D}\mathbf{D}^{\top}\exp\set{-\mathbf{M}^{\top}(t+\tau-t')} \ dt',
\end{aligned}
$$


where in the last line we used that the white noise process $\boldsymbol{\xi}(t)$ has a $\delta$ correlation in time. 

### The Lyapunov equation.

Let $\pmb{\Sigma}:=\ang{\mathbf{r}(t),\mathbf{r}^\top(t)}$ be the (instantaneous, at each $t$) stationary covariance matrix. We may evaluate this quantity algebraically as

$$
\begin{aligned}
	\mathbf{M}\pmb{\Sigma} + \pmb{\Sigma}\mathbf{M}^{\top}
	&=
	\int_{-\infty}^{t}\mathbf{M} 
	\exp\set{-\mathbf{M}\pa{t-t'}} 
	\mathbf{D}\mathbf{D}^{\top}
	\exp\set{-\mathbf{M}^{\top}(t-t')} \ dt' \\
	&\qquad \qquad + 
		\int_{-\infty}^{t} 
		\exp\set{-\mathbf{M}\pa{t-t'}} 
		\mathbf{D}\mathbf{D}^{\top}
		\exp\set{-\mathbf{M}^{\top}(t-t')}\mathbf{M}^{\top}  \ dt'\\
	&= 
	\int_{-\infty}^{t} \frac{d}{dt'}\bra{\exp\set{-\mathbf{M}\pa{t-t'}} 
	\mathbf{D}\mathbf{D}^{\top}
	\exp\set{-\mathbf{M}^{\top}(t-t')}} \ dt' \\
	&= \exp\set{-\mathbf{M}\pa{t-t'}} 
	\mathbf{D}\mathbf{D}^{\top}
	\exp\set{-\mathbf{M}^{\top}(t-t')}\Bigg|_{t'\to -\infty}^{t'=t} \\
	&=\exp\set{0}\mathbf{D}\mathbf{D}^{\top}\exp\set{0} - \lim_{t'\to -\infty} \pa{\exp\set{-\mathbf{M}\pa{t-t'}} 
	\mathbf{D}\mathbf{D}^{\top}
	\exp\set{-\mathbf{M}^{\top}(t-t')}}\\
	&=\mathbf{D}\mathbf{D}^\top.
\end{aligned}
$$

This is called the **Lyapunov equation**. 

> [!REMARK]- Remark
> In the last equality, we take limits 
> $$ 
> \lim_{t' \to -\infty} \exp\set{-\mathbf{M}(t-t')} \qquad \text{and} \qquad \lim_{t' \to -\infty} \exp\set{-\mathbf{M}^{\top}(t-t')}
> $$ 
> which is equivalent to looking at $\lim\limits_{t\to \infty} e^{\mathbf{M}t}$.  Suppose $\mathbf{M}$ is diagonalizable: 
> $$
> \mathbf{M} = \mathbf{U}^{-1}\boldsymbol{\Lambda}\mathbf{U}, \qquad \boldsymbol{\Lambda} = \text{diag}\pa{\lambda_1,\lambda_2,\ldots,\lambda_n}
> $$
> Then 
> $$ e^{t\mathbf{M}} 
> = \mathbf{U}^{-1}e^{t\boldsymbol{\Lambda}}\mathbf{U} 
> = \mathbf{U}^{-1}\pa{\text{diag}\pa{e^{\lambda_1t},e^{\lambda_2t},\ldots,e^{\lambda_n t}}}\mathbf{U}.$$
> This expression shows that the convergence as $t\to \infty$ depends on eigenvalues. In particular, we need their real parts to be $<0$. In this case, $e^{t\mathbf{A}} \xrightarrow[t\to\infty]{} 0$. 
> 
> Because we assume the network is stable, $\mathbf{M}$ must have eigenvalues with positive real part. 
^adf0fd

## Time correlation matrix in stationary state. 

As computed above and using the fact that $\mathbf{DD}^{\top}= \mathbf{M}\bf{\Sigma} + \bf{\Sigma}\bf{M}^\top$, if $t > s$, 

$$
\begin{aligned}
	\ang{\mathbf{r}(t)\mathbf{r}^{\top}(s)}
	&=
	\int_{-\infty}^{s}
	\exp\set{-\mathbf{M}(t-t')}\mathbf{D}\mathbf{D}^{\top}\exp\set{-\mathbf{M}^{\top}(s-t')} \ dt' \\
	&= 
	\int_{-\infty}^{s} \exp\set{-\mathbf{M}(t-t')} \mathbf{M}\boldsymbol{\Sigma} \exp\set{-\mathbf{M}^{\top}(s-t')} \ dt' \\ 
	&\qquad \qquad + \int_{-\infty}^{s}\exp\set{-\mathbf{M}(t-t')} \boldsymbol{\Sigma}\mathbf{M}^{\top} \exp\set{-\mathbf{M}^{\top}(s-t')} \ dt'\\
	&= 
	\int_{-\infty}^{s} \frac{d}{dt'} \bra{\exp\set{-\mathbf{M}(t-t')}\boldsymbol{\Sigma} \exp\set{-\mathbf{M}^{\top}(s-t')}} \ dt' \\
	&= 
	\exp\set{-\mathbf{M}\pa{t-s}}\boldsymbol{\Sigma} - \lim_{t'\to-\infty} \exp\set{-\mathbf{M}\pa{t-t'}} \boldsymbol{\Sigma} \exp\set{-\mathbf{M}^\top\pa{s-t'}}\\
	&= 
	\exp\set{-\mathbf{M}\pa{t-s}}\boldsymbol{\Sigma}.
\end{aligned}
$$

Likewise, if $s > t$, 

$$\begin{aligned}
	\ang{\mathbf{r}(t)\mathbf{r}^{\top}(s)}
	&=
	\int_{-\infty}^{t}
	\exp\set{-\mathbf{M}(t-t')}\mathbf{D}\mathbf{D}^{\top}\exp\set{-\mathbf{M}^{\top}(s-t')} \ dt' \\
	&= 
	\int_{-\infty}^{t} \frac{d}{dt'} \bra{\exp\set{-\mathbf{M}(t-t')}\boldsymbol{\Sigma} \exp\set{-\mathbf{M}^{\top}(s-t')}} \ dt' \\
	&= \boldsymbol{\Sigma}\exp\set{-\mathbf{M}^{\top}(s-t)}.
\end{aligned}$$

> [!NOTE]- Note
> As expected of a stationary solution, the correlation function depends only on time lag $h=|s-t|$. 

Define $\widetilde{\mathbf{C}}(t-s) = \ang{\mathbf{r}(t),\mathbf{r}^{\top}(s)}$. If $h = s-t$,  Then we may write 

$$
\widetilde{\mathbf{C}}(h) = 
\begin{cases}
	\exp\set{\mathbf{M}h}\boldsymbol{\Sigma}, & h <0\\[10pt]
	\boldsymbol{\Sigma}\exp\set{-\mathbf{M}^{\top}h}, & h > 0 
\end{cases}
$$


## Long-time covariance matrix.

The spectrum matrix is the Fourier transform of the correlation matrix, 

$$
\begin{aligned}
	\mathbf{S}(\omega) 
	= 
	\int_{-\infty}^{\infty} e^{-2\pi i\omega h} \widetilde{\mathbf{C}}(h)\ dh.
\end{aligned}
$$

Then $\mathbf{S}(0)=\int_{-\infty}^\infty \mathbf{\tilde{C}}(h) \ dh$. We may also compute this as follows:

$$\begin{aligned}
	\mathbf{S}(0) 
	&= 
	\int_0^{\infty} \boldsymbol{\Sigma}\exp\set{-\mathbf{M}^{\top}h} \ dh + \int_{-\infty}^{0}\exp\set{\mathbf{M}h}\boldsymbol{\Sigma} \ dh \\
	&= 
	\boldsymbol{\Sigma}\pa{\mathbf{M}^{\top}}^{-1} + \mathbf{M}^{-1}\boldsymbol{\Sigma},
\end{aligned}$$

and so 

$$\mathbf{M}\mathbf{S}(0)\mathbf{M}^{\top} = \mathbf{M}\boldsymbol{\Sigma} + \boldsymbol{\Sigma}\mathbf{M}^{\top}= \mathbf{DD}^\top$$

This gives

$$
	\int_{-\infty}^{\infty} \widetilde{\mathbf{C}}(h) \ dh
	=\mathbf{S}(0) 
	= \mathbf{M}^{-1}\mathbf{D}\mathbf{D}^\top\pa{\mathbf{M}^{\top}}^{-1}
	=
	\mathbf{M}^{-1}\mathbf{D}\pa{\mathbf{M}^{-1}\mathbf{D}}^\top.
$$

We define the **long-time covariance** precisely this quantity:

$$
	\mathbf{C}:=\int_{-\infty}^{\infty} \widetilde{\mathbf{C}}(h) \ dh = \mathbf{M}^{-1}\mathbf{D}\pa{\mathbf{M}^{-1}\mathbf{D}}^\top.
$$
