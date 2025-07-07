---
layout: post
title:  "Transition Functions and Probabilities"
date:   2025-07-06 13:37:00 +0100
categories: Functional Analysis
---


<style>
.wrapper {
  max-width: 1500px !important;
}
</style>

{% include mathjax.html %}



$\newcommand{\mbf}{\mathbf}$
$\newcommand{\und}{\underline}$


We explore the analog of markov chains and probability transition matrices for probability density functions. Let $\chi(x,y)$ be the transition function from point $x\in\mathbb{R}^n$ to point $y\in\mathbb{R}^n$. for this function to be a proper probability transition function all of its "colums" must sum to one this means that 

$$
\int_{\mathbb{R}^n} \chi(x,y)dy = 1,\ \forall x\in\mathbb{R}^n
$$

in matrix form this would read 

$$
\sum_j \chi_{ij} = 1, \forall i=1,2,\dots,N
$$

The transition function is aplied to probability density to give back another probability density, in other words, let $\rho(x)$ be a probability density, then 

$$
\theta(y) = \int_{\mathbb{R}^n}\rho(x) \chi(x,y) dx
$$

is a probability density function. To show that it is a density function we check if it integrates to one:

$$
\int_{\mathbb{R}^n} \theta(y)\ dy =  \int_{\mathbb{R}^n}\rho(x) \int_{\mathbb{R}^n}\chi(x,y)\ dy\ dx=\int_{\mathbb{R}^n}\rho(x) dx = 1
$$

For a transition matrix the transition equation is just matrix vector multiplication:

$$\theta_i = \sum_i \rho_i\chi_{ij}$$


What the transition function tells us is that if we are at some point $x$ what is the 'probability' of going to point $y$.

## Discrete time steps

Assume that our probability distribution evolves through the following discrete time equation

$$
\rho(y,k+1) = \int_{\mathbb{R}^n}\rho(x,k) \chi(x,y,k) dx\label{eq:discrte:eq:trans:func}
$$

where $\chi(x,y,k)$ is a transition function for time $k$. So for each time step $k$ the transition function transforms the density distribution into a distribution at time $k+1$. This is the Ananalog to a markov chain where we have the probability update 

$$
\rho(k+1) = \chi(k)\rho(k)
$$

with $\chi(k)\in\mathbb{R}^{N\times N}$ and $\rho(k)\in\mathbb{R}^N$.

## Examples of transition functions for probabilities

A simple example of a transition function is the following 

$$
\chi(x,y) = \mathcal{N}\{ \Sigma(x),\mu(x) \}(y) \equiv (2\pi)^{-n/2}\det(\Sigma(x))^{-1/2}\exp\left(-\tfrac{1}{2}(y-\mu(x))^\top\Sigma^{-1}(x)(y-\mu(x))  \right)
$$

It is easy to check that in fact this function satisfies the properties of a transition function, just note that it is a gaussian function on the variable $y$.


## Particular Case

For the particular case of $n=1$ and $\mathbb{R}^n=\mathbb{R}$ the transition function simplifies to 

$$
\chi(x,y,k) \equiv \frac{\exp\left(-\tfrac{1}{2}\sigma^{-2}(x,k)(y-x)^2  \right)}{\sqrt{2\pi\sigma^2(x,k)}}
$$

where $\sigma^2(x,k)$ is some arbitrary variance and where we set $\mu(x)=x$.


For computations we need to consider a discretization of the spatial domain. Let $X={x_1,x_2,\dots,x_N}$ be the discretization of space. Assume that $\rho(x_i+\Delta_x)=\rho(x_{i+1})$ for $\Delta_x<x_{i+1}-x_i$. That is $\rho$ is assumed constant in the interval $[x_i,x_{i+1}]$. Then the transition function reads

$$
\begin{split}
\rho(y,k+1) &= \int_{-\infty}^{+\infty} \rho(x,k)\chi(x,y,k)dx\\
&= \sum_i \int_{x_i}^{x_{i+1}} \rho(x,k)\chi(x,y,k)dx \\
&= \sum_i \rho(x_i,k)\int_{x_i}^{x_{i+1}}\chi(x,y,k)dx\\
&= \sum_{i}\rho_{ik} \chi_i(y,k)
\end{split}
$$

with $\int_{x_i}^{x_{i+1}}\chi(x,y,k)dx\equiv \chi_i(y,k)$. Evaluating this at some discrete points $y_j$ we arrive at the discrete time and discrete space equation 

$$
\rho(j,k+1) = \sum_i \rho(i,k)\chi(i,j,k)
$$

with $i,j,k$ appropriate integers. The discrete version of $\chi(x,y,k)$ is defined as the following integral

$$
\chi(i,j,k) \equiv \int_{x_i}^{x_{i+1}}\chi(x,y_j,k)dx\approx [\chi(x_{i+1},y_j,k)+\chi(x_{i},y_j,k)]\frac{\Delta_x}{2}
$$

Since the approximation does not ensure that the transition matrix is not a proper probability transition matrix we normalize the columns by defining 

$$
\hat{\chi}(i,j,k)\equiv\frac{\chi(i,j,k)}{\sum_{l}\chi(i,l,k)}
$$

and we use $\hat{\chi}$ instead!

## The integral as a linear operator

As we have seen there is a clear analog between transition functions and transition matrices. Indeed this analog is created by replacing indices $i,j\in\mathbb{N}$ with marks/arguments $x,y\in\mathbb{R}$ and sums over indices by integrals over arguments. Since we have a word for index objects $a_{ij}$ and we have defined how they operate (via matrix-matrix multiplication), we are impeled to define an operator that acts as an integral. We define it as the concatenation operation, concretely, for 'matrix'-'vector' multiplication we define 


$$
\chi \rho \equiv \int_{\mathbb{R}^n} \chi(x,y)\rho(y)\ dy
$$

where $\chi \rho$ is to be understood as $\chi$ acts on $\rho$ as an integral. Let us now call $\chi$ a linear operator and use the symbol $\mathcal{A}^2$ to represent rank-$2$ linear operators, where the $2$ is the number of arguments.  


Now that we have defined the integral as an operator we can do a lot cool shit with it. If we view $\chi$ as a matrix where composition of operators is just the repeated integration, then for $A,B\in\mathcal{A}^2$ the product $AB$ is well defined.

## Differential Equation for the Probability

Recall that the $k$-time update of the probability distribution was given by $\eqref{eq:discrte:eq:trans:func}$. To determine an equation in continuous time $t$ we rewrite $\eqref{eq:discrte:eq:trans:func}$ as

$$
\rho(x,t+\Delta t) = \int_{\mathbb{R}^n} \chi(x,y,t)\rho(y,t)\ dy
$$

subtracting $\rho(x,t)$ on both sides and identifying $\int_{\mathbb{R}^n}\delta(x-y) \rho(y,t)\ dy=\rho(x,t)$, we may write 

$$
(1/\Delta t)(\rho(x,t+\Delta t) - \rho(x,t)) = \frac{1}{\Delta t}\int_{\mathbb{R}^n} (\chi(x,y,t) - \delta(x-y))\rho(y,t)
$$

Now we approximate the left hand side as a time derivative, and we arrive at the equation

$$
\frac{d}{dt} \rho(x,t) = \alpha \frac{1}{\Delta t}\int_{\mathbb{R}^n} (\chi(x,y,t) - \delta(x-y))\rho(y,t) = \alpha(\chi-I)\rho(x,t)
$$


For the particular case where $\chi$ does not depend on time $t$ the equation can be solved by using the exponential function applied to the operator $\alpha(\chi-I)$, then 

$$
\rho(x,t) = \exp(\alpha(\chi-I)t)\rho_0(x)
$$

to show that this is in fact a solution, we take a derivative with respect to time. 

$$
\frac{d}{dt} \rho(x,t) = \alpha(\chi-I)\exp(\alpha(\chi-I)t)\rho_0(x) = \alpha(\chi-I)\rho(x,t)
$$

thus we obtain the desired differential equation

$$
\boxed{
\frac{\partial \rho}{\partial t} = \alpha(\chi-I)\rho
}
$$

Note that it is to be undersood that the exponential of and operator $A\in\mathcal{A}^2$ is defined as follows 

$$
\exp(A) = \sum_{k=0}^\infty \frac{A^k}{k!}
$$

and where $A^k$ are $k$-repeated integration of the function $A(x,y)$, in particular for $k=2$ we have 

$$
(A^2)(y,x) = \int_{\mathbb{R}^n} A(y,z)A(z,x) \ dz
$$

While $A^2$ is the integral of the above over $A(y,x)$ that is 

$$
A^3 = (A^2)A = \int_{\mathbb{R}^n} (A^2)(y,z)A(z,x) \ dz = \int_{\mathbb{R}^n} \int_{\mathbb{R}^n} A(y,z')A(z',z) A(z,x) \ dz'\ dz
$$

given this definition it is also quite trivial to see that for a scalar $\lambda \in\mathbb{R}$ we have $(\lambda A)^k=\lambda^k A^k$, as expected.

## Estimating the transition function from samples



Let us consider samples $x(k)$ from the probability distribution at time $k$, and that the probability transition function does not depend on time $k$. An approach that would yield almost good results would be to define the count function

$$
C(x,y) = \sum_{k=0}^{N-1}g(x(k)-x,x(k+1)-y)
$$

with $g$ an appropriate function.






