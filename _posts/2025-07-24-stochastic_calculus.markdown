---
layout: post
title:  "Stochastic Calculus"
date:   2025-07-24 13:37:00 +0100
categories: Stochastic
---

$\newcommand{\mbf}{\mathbf}$ 
$\newcommand{\und}{\underline}$
$\newcommand{\fdel}{\mathbf{\delta}}$
$\newcommand{\mbb}{\mathbb}$ 
$\newcommand{\bsym}{\boldsymbol}$ 
$\newcommand{\fder}{\mathcal{D}}$ 
$\newcommand{\trace}{\mathrm{Tr}}$ 
$\newcommand{\mcal}{\mathcal}$
$\newcommand{\msf}{\mathsf}$

$\newcommand{\sder}{\mathsf{d}}$


<style>
.wrapper {
max-width: 1500px !important;
}
</style>

{% include mathjax.html %}



The theory of stochastic calculus comes from the ideia of time varying random processes, such processes cannot be described traditionaly via the usual notion of tiny changes with respect to time, because random processes are weird. 

Let $X_t,Y_t$ be random processes that evolve stochastically over time. Consider $dX_t$ a symbolic differential that encodes how a stochastic process changes over an infinitesimally small time interval. Informaly we can describe it as an **infinitesimal random quantity**. Consider some arbitrary function $f:\mbb{R}\mapsto\mbb{R}$ then consider a taylor series expansion around $X_t$ and in the direction $dY_t$ then

$$
f(X_t+dY_t) = \sum_{k=0}^\infty \frac{dY_t^k}{k!}f^{(k)}(X_t) 
$$

We can interpet the value $dY_t^k$ in the sense of an expectation with respect to its probability distribution. This means that in some sense $dY_t^k$ will be a function of $dt$. Consider for instance $Y_t$ a wiener process, this type of process can be described by the distribution function

$$
f_{Y_t}(x) = \frac{1}{\sqrt{2\pi\sigma^2 t}}e^{-\frac{1}{2}\frac{x^2}{\sigma^2t}}
$$

Note that $Y_t$ has moments 

$$
\mbb{E}[Y_t^k] = \sigma^k\sqrt{t}^k
$$

for $k$ a positive odd integer. While for $k$ even it gives 0. Informally we can express the distribution function of $dY_t$ as a function of $dt$, then


$$
f_{dY_t}(x) = \frac{1}{\sqrt{2\pi\sigma^2 dt}}e^{-\frac{1}{2}\frac{x^2}{\sigma^2dt}}
$$

While $dY_t$ has the following moments:

$$
\mbb{E}[dY_t^{2k}] = (\sigma^2 dt)^k
$$

and $\mbb{E}[dY_t^{2k+1}]=0$ for any positive integer $k$. This sugest that $dY_t$ has second order moments of order $dt$ while higher order depend on higher powers of $dt$ which vanish since $dt$ is infinitesimally small. When considering wiener processes Ito's Calculus provides us the following rule 

$$
f(X_t+dY_t) = dY_tf'(X_t) + \frac{1}{2}dY_t^2 f''(X_t)
$$

which reduces to a more generalized chain rule when we try to compute the infinitesimal of $f$. Define 

$$
d(f(X_t)) \equiv f(X_t+dX_t) = dX_t f'(X_t) + \frac{1}{2}dX_t^2 f''(X_t)
$$

A particularly interesting result is when $f$ is an exponential function then 

$$
d(\exp(X_t)) = (dX_t + (1/2)dX_t^2)\exp(X_t)
$$

A question that now arises is why do we not care about higher powers of $dX_t$. I think the answer lies in the fact that we are considering processes such that higher powers of $dX_t$ are just higher orders of $dt$. What would happen, however, if we considered a stochastic process with $dX_t^2$ of order $dt$?? 


### Geometric Brownian Motion with Dividends

Let us start by consider the general linear time varying stochastic equation 

$$
dS(t) = dA(t)S(t) + dU(t)
$$

with $U(t)$ some random process. Our goal is to determine the solution to this time varying stochastic equation. Instead of going through the process of trying to determine the solution, let us start by computing itos derivatives using the Itos chain rule 

Let us consider the integral 

$$
I = \int_{t_0}^t \exp(\Phi(\tau,t)) U(\tau) d\tau
$$

It's derivative with respect to $t$ is 

$$
d I = \left(d\Phi + (1/2)(d\Phi^2)\right) \int_{t_0}^t \exp(\Phi(\tau,t)) dU(\tau) d\tau + dU(t)
$$


where we assumed that $\Phi(t,t)=0$ for all times t. This result sugests that the solution to the stochastic equation should be 

$$
S(t) = \exp(\Phi(t_0,t))S_0 + \int_{t_0}^t \exp(\Phi(\tau,t)) dU(\tau) d\tau
$$

with 

$$\boxed{
d\Phi + (1/2)(d\Phi)^2 = dA(t)
}\label{eq:itos:der:chain:rule}
$$

Of particular interest is the geometric brownian motion which can be expressed as 

$$
dA(t) = \mu(t)dt + \sigma(t)dW_t
$$

Then the function $\Phi(\tau,t)$ will have the particular form 

$$
\Phi(\tau,t) = \int_{\tau}^t \nu(\theta)d\theta + \int_{\tau}^t \sigma(\theta) dW_\theta \label{eq:Phi:from:integral}
$$

with $\nu(t)=\mu(t)-(1/2)\sigma^2(t)$. Taking the ito's derivative of $\eqref{eq:Phi:from:integral}$ we can imediatly write 

$$
d\Phi = \nu(t)dt +  \sigma(t) dW_t
$$

next compute $\eqref{eq:itos:der:chain:rule}$ to obatin

$$
\begin{split}
d\Phi + (1/2)(d\Phi)^2 &= \nu(t)dt +  \sigma(t) dW_t + (1/2)(\nu^2(t)dt^2 + \sigma^2(t)dW_t^2 + 2\nu(t)dt\sigma(t) dW_t)\\
                       &= (\mu(t)-(1/2)\sigma^2(t))dt + \sigma(t) dW_t + (1/2)\sigma^2(t)dW_t^2\\
                       &= (\mu(t)-(1/2)\sigma^2(t))dt  + \sigma(t) dW_t + (1/2)\sigma^2(t)dt 
\end{split}
$$

where we used $dW_t^2=dt$ and $dt^2=dW_tdt=0$ since $dW_t$ is of order $\sqrt{dt}$. Which makes us finally see that 

$$
dA = d\Phi + (1/2)(d\Phi)^2 = \mu(t)dt + \sigma(t)dW_t
$$

as required. An interesting particular case of the stochastic equation arises when we set $dU(t)=-dA(t)d(t)$, this particular case enables us to describe the paying of dividends at continuous points in time. With this choice the solution becomes of the form 

$$
S(t) = \exp(\Phi(t_0,t))S_0 - \int_{t_0}^t \exp(\Phi(\tau,t)) dA(\tau) d(\tau) d\tau\label{eq:divs}
$$


An even more trivial case is when we consider no dividends that is $d(t)=0$ for all time $t$. We can interpet the integral in $dW_t$ as infinitesimal sum of standard normal distributed functions this sugests the following

$$
\int_{\tau}^t \sigma(\theta)dW_\theta = \sqrt{\int_\tau^t \sigma^2(\theta)d\theta}\  Z(\tau,t)
$$

with $Z(\tau,t)$ a random normal variable. Note however that $Z(t,\tau)$ and $Z(t',\tau')$ are statistically independent only if the intervals $[\tau,t]$ and $[\tau',t']$ are disjoint. Let us define 

$$
\Sigma(\tau,t) = \sqrt{\int_\tau^t \sigma^2(\theta)d\theta}
$$

and 

$$
\bar{\nu}(\tau,t) = \int_{\tau}^t \nu(\theta)d\theta
$$


Then we can write the particular case where $d(t)=0$ in the form

$$
S(t) = \exp\left(\bar{\nu}(t_0,t)+\Sigma(t_0,t)Z(t_0,t)\right)S_0
$$

Then the estimated payoff value may be computed by determining 

$$
P(t) = \mbb{E}[(S(t)-K)_+]
$$

with $(X)_+=\max(0,X)$. Determining this expected value is straightforward since it only involves a single random variable $Z(t_0,t)$ which is normally distributed. However when we consider the solution with non-zero dividends i.e. equation $\eqref{eq:divs}$, computing the expected value then becomes rather difficult.

<!--Note however that $d[\exp(\Phi)]=\exp(\Phi)dA$, this suggest that we may rewrite the above equation. 
Let us consider a discretization of this-->

## Discretizing the Dividends

We might be able to provide an heuristic way to determining dividends by considering a discretization. In particular, assume that $d(t)$ is constant in the interval $t\in[t_k,t_{k+1}]$


We start by considering rewriting the solution to the stochastic differentia equation in the generic form

$$
S(t_1,t_2) = \exp(\Phi(t_1,t_2))S(t_0,t_1) +\int_{t_1}^{t_2} \exp(\Phi(\tau,t_2))dA(\tau)d(\tau)\ d\tau
$$

with $t_0<t_1<t_2$. This equation tells us how given some initial condition at time $t_1$ we can determine the value of the stock $S$ at time $t_2$. Note that using $\eqref{eq:itos:der:chain:rule}$ and Ito's derivative we have $d[\exp(\Phi)]=\exp(\Phi)dA$, this means that we can rewrite the above integral and thus write the solution in the form


$$
S(t_1,t_2) = \exp(\Phi(t_1,t_2))S(t_0,t_1) +\int_{t_1}^{t_2} \sder\left[\exp(\Phi(\tau,t_2))\right]d(\tau)\ d\tau\label{eq:stock:t1t2}
$$


with $\sder$ the stochastic derivative. Now consider that $d(t)$ is constant in the interval $t\in]t_1,t_2]$ then note that the integral becomes

$$
\int_{t_1}^{t_2} \sder\left[\exp(\Phi(\tau,t_2))\right]d(\tau)\ d\tau = \int_{t_1}^{t_2} \sder\left[\exp(\Phi(\tau,t_2))\right]\ d\tau\ d(t_2) = (\exp(\Phi(t_1,t_2)) - 1)d(t_2)
$$

where we used $\Phi(t_1,t_1)=0$. Then equation $\eqref{eq:stock:t1t2}$ takes the form 

$$
S(t_k,t_{k+1}) = \exp(\Phi(t_k,t_{k+1}))(S(t_{k-1},t_k) + d(t_{k+1})) - d(t_{k+1}) 
$$

for $0<t_0<t_1<t_2<\cdots<t_n<T$ positive real numbers. Next we define the following discrete time varying random variables

$$
\Pi_k = \exp(\Phi(t_k,t_{k+1})),\ \ \ \Omega_k = \Pi_k-1,\ \ \ S_k = S(t_{k-1},t_k),\ \ \ d_k = d(t_{k+1}),\ \ \ \Phi_k = \Phi(t_k,t_{k+1})
$$

as such we arrive at the following beautifull stochastic discrete equation 

$$
S_{k+1} = \Pi_kS_k+\Omega_k d_k
$$

Note that for each $k$ each $\Phi_k$ is computed from disjoint time intervals, which means that $\Phi_k$ is statistically independent from $\Phi_j$ for all $k\neq j$. Consequently $\Pi_k$ and $\Omega_k$ are statistically independent from $\Pi_j$ and $\Omega_j$ for all $k\neq j$. Statisticall independence means that the joint probability distribution function of $\Pi_0,\Pi_2,\dots,\Pi_n$ is decomposabled as the product of each probability. Even though the joint distribution function is decomposable the problem of computing expected expected payoffs from $S_k$ is quite complicated, to adress such a difuicult problem, let us start by considering a single step of the discrete equation, then write

 
$$
S = \exp(\nu +\sigma Z)(s+d) - d
$$

The expected payoff is 

$$
\mbb{E}[(S-K)_+] = (s+d)e^{\nu}e^{\frac{1}{2}\sigma}[1-F_Z((1/2)\sigma-\alpha)] - (K+d)[1-F_Z(\alpha)]
$$

with 

$$
\alpha= \frac{\log\left(\frac{d+K}{d+s}\right)-\nu}{\sigma}
$$

Note that this solution is deterministic, however when multiple time steps are involved $s$ will become a function of multiple random variable and, because $\alpha$ is a function of $s$ and we need to compute the $F_Z(\dot)$ the problem of determining the expected value becomes way more difficult in just the second time step.
