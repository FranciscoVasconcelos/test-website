---
layout: post
title:  "Functionals and Linear Operators"
date:   2025-07-06 13:37:00 +0100
categories: Functional Analysis
---


<!-- To not have a break on this part put all the newcommand's in the same line -->


$\newcommand{\mbf}{\mathbf}$ 
$\newcommand{\und}{\underline}$
$\newcommand{\fdel}{\mathbf{\delta}}$
$\newcommand{\mbb}{\mathbb}$ 
$\newcommand{\bsym}{\boldsymbol}$ 
$\newcommand{\fder}{\mathcal{D}}$ 
$\newcommand{\mcal}{\mathcal}$
$\newcommand{\msf}{\mathsf}$

$\newcommand{\matrixtrace}{\text{mtr}}$ <!-- Matrix Trace -->
$\newcommand{\ctrace}{\text{ctrace}}$ <!-- Complete trace: Operator and matrix trace -->
$\newcommand{\trace}{\mathrm{otrace}}$  <!-- Operator Trace -->

$\newcommand{\infd}{\mathbf{d}}$ <!-- Infinitesimal Operator -->
$\newcommand{\infD}{\mathrm{D}}$ <!-- Infinitesimal Volume Operator -->

<style>
.wrapper {
max-width: 1500px !important;
}
</style>

{% include mathjax.html %}


TODO:
- Clarify/prove Taylor series approximation for functionals
- Prove identity $\eqref{eq:der:of:der:identity}$
- Need to clarify notation between adjoint and differential operators of euclidean space and of functional space (Hilbert spaces and Banach spaces)

* TOC
{:toc}


The main goal of this post is to build a mathematical machinery that provides the user(i.e. the matematician) with clear and simple tools to do functional calculus. A key distinction of my approach is that I avoid the ideia of basis functions or orthogonal basis functions, however for some types of problems the basis function description can simplify computations, either on physical computer or analyticaly. One of my goals is to be able to relate ideas from linear algebra to the space of linear operators over an infinite dimensional vector space (which are the functions). The reader should start to get aquanted to the ideia of functions being interpeted as vectors with infinite basis elements indexed by infinitesimal components. 

## Finite Dimensional Euclidean Spaces ($\mbb{R}^n$)

For the reader, $n$-dimensional euclidean vector spaces $\mcal{E}^n$ are mostly familiar. When talking about euclidean spaces we like to think about a set of basis vector $\{e_1,\dots,e_n\}$ $e_i\cdot e_j=\delta_{ij}$ that span the entire vector space $\mcal{E}^n$ and any vector that belongs to this vector space can be expressed as a linear combination of its components, namely for $v\in\mcal{E}^n$ we have 

$$
v = v^1e_1 + v^2e_2 + \cdots + v^ne_n
$$

however the representation of this vector space through this basis is not unique, we can find a different set of basis from which we shall be able to describe the vector $v$. Essentialy the vector space satisfies the following trivialities

$$
\begin{align}
u+v &\in \mcal{E}^n,\\
\alpha u &\in \mcal{E}^n
\end{align}
$$

for any $u,v\in\mcal{E}^n$. There are a set of notions that are attached to the ideia of an Euclidean space. One is the notion of distance and length, but first we need to define an inner product $u\cdot v$ from which we can reclaim the notion of angles and distance. If we express the vectors $u$ and $v$ through the basis vectors. Then the inner product can be expressed as 

$$
u\cdot v = (u^1,u^2,\dots,u^n)\cdot (v^1,v^2,\dots,v^n) = \sum_k u^k v^k
$$

While the notion of angle can be expressed via the expression $u\cdot v = \|u\|\|v\|\cos(\theta)$ the notion of magnitude or distance can be expressed as 

$$
|u| = \sqrt{u\cdot u}
$$

The most useful tool in linear algebra are matrices. Matrices act on vector to produce other vectors, and they are linear on their argument, we can write a linear operation by the simple use of the inner product, let $a\in\mcal{E}^n$ then the following is an example of a linear operation by applying the matrix to the vector $v$:

$$
Av = au\cdot v
$$

$$
Av = au\cdot v
$$

with $A:\mcal{E}^n\mapsto \mcal{E}^n$ a matrix that takes $\mcal{E}^n$ and transforms it into $\mcal{E}^n$. Matrix vector multiplication can also be computed with the use of indexation in particular we may write 

$$
u_i = \sum_{j=1}^n A_{ij}v_j\label{mat:mult:components}
$$

while matrix-matrix multiplication can be written in component form as the following

$$
C_{ij} = \sum_{k=1}A_{ik}B_{kj}\label{eq:cij:mat:mul}
$$

We are being carefull to include all of this indexed operations because we want to make the bridge to the functional operators that will have continuous rather then integer indices.


A matrix $U\mbb{R}^{n\times n}$ is orthogonal when it satisfies the orthogonal condition $U^\top U=I$ with $I$ the identity matrix. In index notation the orthogonality condition can be expressed as 

$$
u_i^\top u_j = \sum_{k=1}^n U_{ik}U_{jk} = \delta_{ij}\label{eq:ortho;cond:lin:alg}
$$

where the $u_i$ are columns of $U$. 

## Function Spaces and Functional Operators

To be able to understand fully what a function space is we can relate it to our dear familiar finite dimensional euclidean space. A function space is the set of all objects that transform vector space $\mcal{V}^n$ into other vector spaces $\mcal{U}^m$ we denote this set a 

$$
\mcal{F}(\mcal{V}^n,\mcal{U}^m) = \{ f: \mcal{V}^n\mapsto\mcal{U}^m\} 
$$

Similarly to how functions act on vector spaces, functionals or functional operators ${F}\in\mcal{G}=\\{G:\mcal{F}(\mcal{V}^n,\mcal{U}^m)\mapsto \mcal{F}(\mcal{V}^n,\mcal{U}^m)\\}$ act on functions to produce other functions. The most simple kind of functionals are linear functionals or also called linear operators, as the name sugests they have the property that they transform functions linearly. To reduce notation and since for simplicity we want to focus on a particular case of function spaces let us consider some arbitrary vector space denoted $\Omega$ where the domain of the functions is defined, and the codomain to be $\mbb{R}$, then we shall use 

$$
\mcal{F}({\Omega}) =\mcal{F}_{\Omega} \equiv \{ f:\Omega\mapsto \mbb{R} \}
$$

That is $$\mcal{F}_{\Omega}$$ is the space of all functions that transform $\Omega$ into $\mbb{R}$. The space of linear functionals on $\mcal{F}_\Omega$ is denoted by $L(\Omega)\subset \mcal{G}(\Omega)$. Where $\mcal{G}(\Omega)$ is defined as 

$$
\mcal{G}(\Omega) \equiv \{ G:\mcal{F}_\Omega \mapsto \mcal{F}_\Omega \}
$$

- We use $\nabla$ to express the gradient with respect to arguments of functions
- We use $\fder$ to express the functional derivative
- We use greek letters for functions $\psi,\phi,\rho\in  \mcal{F}_{\Omega}$,
- Upper case and lower case latin letters for functionals e.g. $F,f,g,h,J\in\mcal{G}(\Omega)$. 
- Bold symbols for linear functionals or linear operators $\mbf{A,B,C,D}\in L(\Omega)$
- The last letters of the alphabet to represent arguments of functions $x,y,z,w\in\Omega$ or $x_1,x_2,\dots,x_m\in\Omega$ 
- The symbol ${}^\dagger$ is used to denote the transmutation operation over function spaces (is equivalent to matrix transpose)
- Kernel functions have the same letter as associated linear functionals, but non bold $A(x,y),B(x,y),C(x,y),D(x,y)\in \mcal{F}(\Omega\times \Omega)$


To a linear operator $\mbf{A}\in L(\Omega)$ we can always associate an integral operator, with an associated kernel function $A\in\mcal{F}(\Omega\times\Omega)$. Operating on some function with the linear operator $\mbf{A}$ is equivalent to integrating over the vector space $\Omega$, then 

$$
\phi(x) = \mbf{A}\psi = \int_{\Omega} A(x,y)\psi(y)\ dy 
$$

See the similarity with $\eqref{mat:mult:components}$? Just replace the index $i$ by $x$, the index $j$ by $y$ and the sum by the integral and voila we have the functional linear operation defined above. The transmutation ${}^\dagger$ operation can be defined for both function and linear operators, we may define it through its properties 

$$
\begin{align}
\psi^\dagger \mbf{A} &= (\mbf{A}^\dagger\psi)^\dagger\\
\psi^\dagger\phi &= \int_{\Omega} \psi(x)\phi(x) dx\in\mbb{R}\\ 
\phi\psi^\dagger &\in L(\Omega) \\
(\mbf{AB})^\dagger &= \mbf{B}^\dagger\mbf{A}^\dagger
\end{align}
$$

The last equation is similar to the outer product of vectors while the second is the inner product. Setting $\mbf{A}=\phi\psi^\dagger$ and applying $\mbf{A}$ to some function $\rho\in\mcal{F}(\Omega)$ we have 

$$
\mbf{A}\rho = \int_\Omega A(x,y)\rho(y)\ dy = \int_{\Omega} \phi(x)\psi(y)\rho(y)\ dy = \phi(x)\int_{\Omega} \psi(y)\rho(y)\ dy = \phi (\psi^\dagger\rho)
$$

So the outer product of functions produces a linear operator! The transmutation operation of a linear operators swaps arguments of the associated kernel function $A(x,y)$ we can express this by the following 

$$
\mbf{A}^\dagger\psi \equiv \int_{\Omega} A(y,x)\psi(y)\ dy 
$$

Multiplying linear operator is equivalent to integrating their respective kernel functions. Let $\mbf{C}=\mbf{AB}$ then a relation for the kernel functions happens to be 

$$
C(x,y) = \int_{\Omega} A(x,z)B(z,y)\ dz
$$

Look at how similar this expression is to $\eqref{eq:cij:mat:mul}$!!!


An orthogonal operator $\mbf{U}\in\mathcal{O}(\mcal{F}_\Omega)$ or orthogonal functional is a linear functional that satisfies the orthogonality condition

$$
\mbf{U^\dagger U} = \mbf{I}
$$

where $\mbf{I}$ is the identity operator, the associated kernel function is $I(x,y)=\delta(x-y)$ where $\delta(x)\in\mcal{F}_\Omega$ is the dirac-de-delta funcion defined over the vector space $\Omega$. In terms of integrals the orthogonality condition tells us that 

$$
\int_{\Omega} U(z,x)U(z,y)\ dz = \delta(x-y)
$$

Again, see the similarity with the linear algebra expression $\eqref{eq:ortho;cond:lin:alg}$!? Try to replace indices by components from $\Omega$ and the sum by an integral!

*Projection Operators* are symmetric and indepotent operators meaning for $\mbf{P}$ an projection operator, we have

$$
\begin{split}
\mbf{P}^\dagger &= \mbf{A}\\
\mbf{P}^2 &= \mbf{PP}^\dagger = \mbf{P}
\end{split}
$$

The most trivial projection operator is given as the outer product of two degree-one functions. 

$$
\mbf{P} = \psi\psi^\dagger
$$

with $\psi^\dagger\psi=1$. We can also form projection operators by integrating over some parameterization of degree-one functions $\psi_\lambda$

$$
\mbf{P} = \int_{\mbb{R}} \psi_x\psi_x^\dagger dx
$$

Note however in order for this to be a projector, each $\psi_x$ must satisfy 

$$
\begin{align}
\psi_x^\dagger\psi_{x'} &= \delta(x-x'), \ \forall x,x' \in {S}\\ 
\|\psi_x\| &= 0,\ \ \forall x\notin S
\end{align}
$$

with $S\subset \mbb{R}$.



$$
\mbf{P}^2 = \int_{S} \psi_x\psi_x^\dagger dx\int_{S} \psi_{x'}\psi_{x'}^\dagger dx' = \int_S\int_S \psi_x\psi_x^\dagger\psi_{x'}\psi_{x'}^\dagger dx dx' =  \int_S\int_S \psi_x\delta(x-x')\psi_{x'}^\dagger dx dx' = \mbf{P}
$$

to show that it is symmetric we do 

$$
\mbf{P}^\dagger = \int_{\mbb{R}}(\psi_x\psi_x^\dagger )^\dagger dx = \mbf{P}
$$

Any symmetric operator $\mbf{A}$ can be written as the integral 

$$
\mbf{A} = \int_{\mbb{R}}\lambda_x \psi_x\psi_{x}^\dagger dx
$$

we call the above the spectral decomposition of $\mbf{A}$, and we call $\psi_x$ the spectral functions of $\mbf{A}$. Note that each $\psi_x$ satisfies the eigenvalue problem 

$$
\mbf{A}\psi_x = \lambda_x \psi_x
$$

The determinant of the symmetric operator can be determined from the eigenvalues $\lambda_x$, in particular we have 

$$
\log\det(\mbf{A}) = \int_{\mbb{R}}\log(\trace(\lambda_x \psi_x\psi_{x}^\dagger))\ dx = \int_S \log(\lambda_x)\  dx
$$

there is another interesting way to determine the determinant of an operator. We use the following 

$$
\det(e^{\mbf{B}}) = e^{\trace({\mbf{B}})}\label{eq:det:from:trace}
$$

where $\trace(\cdot)$ is the operator trace (see [Operator Trace](#operator-trace)). 

#### Operator exponential

We can deterimine an operator $\mbf{B}$ such that $e^{\mbf{B}} = \mbf{A}$. First note that 

$$
\mbf{A}^2 = \int_{\mbb{R}}\int_{\mbb{R}}\lambda_x\lambda_{x'} \psi_x\psi_{x}^\dagger \psi_{x'}\psi_{x'}^\dagger\  dx\ dx' = \int_\mbb{R} \lambda_x^2\psi_x\psi_{x}^\dagger \ dx
$$

and more generaly we have 

$$
\mbf{A}^k = \int_\mbb{R} \lambda_x^k\psi_x\psi_{x}^\dagger \ dx
$$

for any positive integer $k$. Then assume that $\mbf{B}$ has the decomposition 

$$
\mbf{B} = \int_{\mbb{R}} \alpha_x\phi_{x}\phi_{x}^\dagger \ dx
$$

Now use the exponential power series to write 

$$
\begin{split}
\exp(\mbf{B}) &= \sum_{k=0}^\infty \frac{\mbf{B}^k}{k!} = \sum_{k=0}^\infty \int_{\mbb{R}} \frac{\alpha_x^k}{k!}\phi_{x}\phi_{x}^\dagger \ dx \\
&= \int_{\mbb{R}} \sum_{k=0}^\infty \frac{\alpha_x^k}{k!}\phi_{x}\phi_{x}^\dagger \ dx = \int_{\mbb{R}} \exp(\alpha_x) \phi_{x}\phi_{x}^\dagger dx
\end{split}
$$

Finally consider that 

$$
\mbf{A} = e^\mbf{B}
$$

which suggests that $\phi_x=\psi_x$ and $\lambda_x=\exp(\alpha_x)$. Now use $\eqref{eq:det:from:trace}$ to write 

$$
\begin{split}
\det(\mbf{A}) &= \exp(\trace(\mbf{B})) =  \exp\left(\trace\left(  \int_{\mbb{R}} \alpha_x\phi_{x}\phi_{x}^\dagger \ dx \right)\right) \\
&= \exp\left(\int_{\mbb{R}} \trace(\alpha_x\phi_{x}\phi_{x}^\dagger)\right) = \exp\left(\int_{R} \alpha_x\ dx \right) \\
&= \exp\left(\int_{S} \log(\lambda_x)\ dx \right) 
\end{split}
$$

where  $R=S$, because $\psi_x=\phi_x$ and where $\alpha_x=\log(\lambda_x)$. 


## The trace of operators 
{:#operator-trace}
We define the operator trace $\trace(\cdot)$ as 

$$
\trace(\mbf{A}) = \int_{\mbb{R}} A(x,x)\ dx
$$

the *Operator Trace* has the following properties

$$
\begin{align}
\trace(\alpha \mbf{A}) &= \alpha\!\ \trace(\mbf{A})\\
\trace(\mbf{A}+\mbf{B}) &= \trace(\mbf{A}) + \trace(\mbf{B})\\
\trace(\mbf{A}^\dagger) & = \trace(\mbf{A})\\
\trace(\mbf{AB}) &= \trace(\mbf{BA})\\
\trace(\alpha) &= \alpha 
\end{align}
$$

As a special case of the above we have

$$
\trace(\psi\phi^\dagger) = \phi^\dagger\psi = \int_{\mbb{R}} \phi(x)\psi(x)\ dx
$$

## Skew, Symmetric and Orthogonal Operators

Any linear operator $\mbf{A}$ can be wirtten as the summ of a symmmetric part $\mbf{A}\_+$ and a skew-symmetric part $\mbf{A}\_-$, *i.e.*

$$
\mbf{A} = \mbf{A}_++\mbf{A}_-
$$

with 

$$
\mbf{A}_+ = \frac{1}{2}(\mbf{A+A^\dagger}),\ \ \ \mbf{A}_-  =\frac{1}{2}(\mbf{A-A^\dagger})
$$

We say that an operator $\mbf{A}$ is *symmetric* if it satisfies $\mbf{A}^\dagger=\mbf{A}$, and we say that an operator $\mbf{A}$ is *skew-symmetric* if it satisfies $\mbf{A}^\dagger=\mbf{A}$. In terms of the kernel $A(x,y)$ of the operator $\mbf{A}$ this means that 

$$
\begin{split}
A(x,y)  = A(y,x) \quad \quad &(\text{symmetric operator})\\
A(x,y)  = -A(y,x) \quad \quad &(\text{skew-symmetric operator})
\end{split}
$$

Skew-symmetric operators generate orthogonal operators. For $\mbf{A}$ a skew-symmetric operator we can easily show that 

$$
\mbf{U} = \exp(\mbf{A})
$$

is an orthogonal linear operator. An orthogonal operator satisfies 

$$
\mbf{U}^\dagger\mbf{U} = \mbf{UU}^\dagger  = \mbf{I}
$$


To show that this is indeed true, recall that the exponential can be written as taylor series expansion and that $(\mbf{AB})^\dagger=\mbf{B}^\dagger\mbf{A}^\dagger$, then

$$
\mbf{U}^\dagger = \exp(\mbf{A})^\dagger = \sum_{k=0}^\infty \frac{\left(\mbf{A}^k\right)^\dagger}{k!} = 
\sum_{k=0}^\infty \frac{\left(\mbf{A}^\dagger\right)^k}{k!} = \exp(\mbf{A}^\dagger) = \exp(-\mbf{A}) = \mbf{U}^{-1}
$$


This suggests that we can generate any orthogonal linear operator by considering some linear operator $\mbf{A}$ with arbitrary kernel function, compute its skew-symmetric part and evaluate $\mbf{U}=\exp(\mbf{A}_-) = \exp\\!\left(\tfrac{1}{2}(\mbf{A-A}^\dagger)\right)$




## Probability Densities

In this section I explore the analog of markov chains and probability transition matrices for probability density functions. With the use of special types of functional operators we can consider linear functional transformations that transform probability density functions to other probability density functions. A probability density function $$\rho\in\mcal{F}_\Omega$$ has the property that is positive and that integrates to one. In order for the transition probability operator $\mbf{A}$ to preserve this quantity there are certain properties that the kernel function must satisfy. Let us say that a function is normalized with respect to some function $\mbf{a}\in\mcal{F}_\Omega$ when it satisfies 


$$
\mbf{a}^\dagger\rho = 1
$$

If this is satisfied we say that $\rho$ is $\mbf{a}$-normalized. In the context of probability densities the function $\mbf{a}$ is simply one, we will denote the one function as $\mbf{i}$ analogous to the vector of all ones. Note that we assume the following notation $\mbf{i}^\dagger \rho = \int_{\Omega} 1\rho(x)\ dx$. So we care to preserve $\mbf{i}$-normality, assume that $\rho$ is $\mbf{i}$-normalized, that is it satisfies $\mbf{i}^\dagger\rho=1$, then we want that $\theta=\mbf{A}\rho$ to be $\mbf{i}$-normalized, for any $\mbf{i}$-normalized function $\rho$ this implies the following:

$$
\mbf{i}^\dagger\theta = \mbf{i}^\dagger \mbf{A}\rho = 1
$$

which is satisfied if $\mbf{A}^\dagger\mbf{i}=\mbf{i}$. In other words the kernel function $A(x,y)$ must integrate to one, that is 

$$
\int_{\Omega} A(y,x) dy = 1 \ \text{for all } x\in\Omega. 
$$

## Discrete Time Updates

Probability distribution functions usually evolve through time $t\in\mbb{R}$, however we want to consider a discretization of time, let $k\in\mbb{N}$ denote the discrete time step, then a distribution that evolves through discrete time $k$ must satisfy the following difference equation 

$$
\rho(k+1) = \mbf{A}(k)\rho(k)
$$

with $\mbf{A}(k)\in L(\Omega)$, $\rho(k)\in \mcal{F}(\Omega)$ for all discrete time steps $k\in\mbb{N}$. This equation can be solved when given initial conditions $\rho(k_0)=\rho_0\in\mcal{F}_\Omega$, recalling that operator multiplication is well defined via the respective integration of the kernel functions we have the solution:

$$
\rho(k) = \mbf{A}(k-1)\cdots\mbf{A}(k_0)\rho_0
$$

Other than the $\mbf{i}$-normalization condition that $\rho(k)$ must satisfy to be a proper probability density we also have to ensure non-negativity. The trivial choice $A(x,y)\geqslant 0$ $\forall x,y\in\Omega$ does ensure non-negativity. And it is the one used for [Discrete Time Markov Chains](https://en.wikipedia.org/wiki/Markov_chain#Discrete-time_Markov_chain). An operator that satisfies $\mbf{i}^\dagger\mbf{A}=1$ and $\mbf{A}\geqslant 0$ is called a right stochastic linear operator.  

## Continuous time updates 

If we are expressing the evolution in continuous time we need to consider a differential equation on time $t\in\mbb{R}$. Such an equation is 

$$
\dot{\rho}(t) = \mbf{B}(t)\rho(t)
$$

for some linear functional $\mbf{B}\in L(\Omega)$. The solution to this equation is non-trivial but it can be expressed via the [Peano-Baker Series](https://en.wikipedia.org/wiki/State-transition_matrix#Peano%E2%80%93Baker_series). However we want to make sure that the solution is $\mbf{i}$-normalized, as such we want to determine what are the conditions on $\mbf{B}(t)$ in order for the probability $\rho(t)$ to stay normalized. One such condition is by considering that $\mbf{i}^\dagger\rho(t)$ is conserved, in other words 

$$
\mbf{i}^\dagger \rho(t) = 1\ \ \Leftrightarrow \ \ \mbf{i}^\dagger \dot{\rho} = 0
$$

Where we simply took the derivative on both sides with respect to time $t$. This implies that we must have 

$$
\mbf{i}^\dagger\mbf{B}(t)\rho(t) = 0\  \text{for all}\ t\in\mbb{R} 
$$

which can be achieved by requiring that  $\mbf{i}^\dagger\mbf{B}(t)=0$ for all $t\in\mbb{R}$. Note that while this is satisfied the exponential of $\mbf{B}t$ will satisfy 

$$
\mbf{i}^\dagger e^{\mbf{B}t} = \mbf{i}^\dagger \sum_{k=0}^\infty \frac{\mbf{B}^kt^k}{k!} = \mbf{i}^\dagger \left( \mbf{I} + \mbf{B}+\tfrac{1}{2}\mbf{B}^2+\cdots \right) = \mbf{i}^\dagger 
$$

since $\mbf{i}^\dagger$ multiplied with any non-zero power of $\mbf{B}$ yields the zero function. This means that for the transition matrix $\mbf{B}$ a constant function of time the solution to the differential equation is guaranteed to preserve $\mbf{i}$-normalization. Recall that for $\mbf{B}=\textrm{const}(t)$ the solution is 

$$
\rho(t) = e^{\mbf{B}t}\rho_0
$$

multiplying by $\mbf{i}^\dagger$ on both sides gives us 

$$
\mbf{i}^\dagger \rho(t) = \mbf{i}^\dagger e^{\mbf{B}t}\rho_0 = \mbf{i}^\dagger\rho_0 = 1
$$

when $\rho_0$ is $\mbf{i}$-normalized. The other condition that must be ensured is non-negativity of the probability density thus for the generator operator $\mbf{B}$ we must ensure that its kernel satisfies

$$
B(x,y) \geqslant 0,\ \forall x\neq y\in\Omega\label{eq:B(x,y):non-negative}
$$

this means that the 'diagonal' entries of $\mbf{B}$ may be negative while all other entries must non-negative. This non-negativity property ensures that the transition function $\mbf{P}(t)=e^{\mbf{B}t}$ has a non-negative kernel for all times $t\in\mbb{R}$. This continuous time description of a difusion process that evolves in time $t$ is similar to [Continuous Time Markov Chains](https://en.wikipedia.org/wiki/Markov_chain#Transition_probability_definition). This analogy between these two concepts, nicely interconnect ideas from two seperate worlds. Through this analogy we can deduce that if the transition matrix of continuous time markov chain is non-negative when the generator matrix has off-diagonal non-negative entries, then it also must be true that  $\mbf{P}(t)=e^{\mbf{B}t}\geqslant 0$ for the kernel of $\mbf{B}$ to satisfy $\eqref{eq:B(x,y):non-negative}$.

From now on when we write an inequality with respect to a function or operator we are refering to an inequality with respect to its kernel. In particular we are considering the following analogies

$$
\begin{align}
\psi > 0 \ \Leftrightarrow\ \ \psi(x)>0\ \forall x\in\Omega\notag \\ 
\mbf{A} \geqslant 0 \ \Leftrightarrow\ \ A(x,y) \geqslant 0\ \forall x,y\in\Omega  \notag 
\end{align}
$$


## Derivatives as Integrals

One of the interesting aspects of the formulation of linear operators via kernel functions is that it can also encode/represent derivatives. In particular with the use of the dirac-delta function and its derivatives we can compute kernel functions from derivatives, consider for instance the general derivative operator 


$$
\mbf{D} = \sum_{k=0}^n a_k(x)\frac{d^k}{dx^k}
$$

we shall show that the kernel of $\mbf{D}$ is 

$$
D(x,y) = \sum_{k=0}^n (-1)^ka_k(x) \delta^{(k)}(x-y)
$$

with $\delta^{(k)}$ the $k$-th derivative of the dirac-delta function. To show that $D(x,y)$ is indeed the kernel function of $\mbf{D}$ we will show that 

$$
\int_{\mathbb{R}} \delta^{(k)}(x-y)\psi(y)dy = (-1)^k\frac{\partial^{k}\psi}{\partial x^k}\label{eq:dirac:delta:derivative}
$$

To show how this holds recall the classic integration by parts formula

$$
\int_{I} \psi \phi'\ dx = [\psi\phi]_{\partial{I}} - \int_{I}\psi'\phi\ dx
$$

with $I\subset\mathbb{R}$ some interval and $\partial I$ its boundary. Applying this result to $\psi(x)\delta^{(k)}(x-y)$ yields the formula 

$$
\int_{\mathbb{R}} \delta^{(k)}(x-y)\psi(x)\ dx = [\delta^{(k)}(x-y)\psi(x)]_{-\infty}^{+\infty} - \int_{\mathbb{R}}\psi'(x)\delta^{(k-1)}(x-y)\ dx = - \int_{\mathbb{R}}\psi'(x)\delta^{(k-1)}(x-y)\ dx
$$

where $\delta^{(k)}$ plays the role of $\phi$. To obtain the final result we made $[\delta^{(k)}(x-y)\psi(x)]_{-\infty}^{+\infty}$ disapear, we argue that either $\psi$ vanishes at $\pm\infty$ or that $y$ will never attain the value of $\pm\infty$. Applying repeated integration by parts $k$-times will yield the following 


$$
\int_{\mathbb{R}} \delta^{(k)}(x-y)\psi(x)\ dx = (-1)^k \int_{\mathbb{R}} \delta(x-y) \frac{\partial^k\psi}{\partial{x}^k}\ dx
$$

which after evaluating the integral on the right will give $\eqref{eq:dirac:delta:derivative}$.


We have to be carefull when dealing with the kernel functions of derivative operators. Let us consider the particular case when we have an operator $\mbf{D}$ with kernel 

$$
D(x,y) = -\frac{d}{dx}\delta(x-y) = \frac{d}{dy}\delta(x-y)
$$

Now we use $\mbf{D}$ on a test function $\psi$, then 

$$
\begin{split}
[\mbf{D\psi}](x) &= \int_{\mbb{R}} \left(\frac{d}{dy}\delta(x-y) \right)\psi(y)\ dy\\
&= -\int_{\mbb{R}} \delta(x-y)\frac{d\psi}{dy}\ dy + [\delta(x-y)\psi(y)]_{-\infty}^{+\infty}\\
&= -\frac{d\psi}{dx}
\end{split}
$$

next we use the transmutated operator $\mbf{D}^\dagger$ to determine 

$$
\begin{split}
[\mbf{D^\dagger\psi}](y) &= \int_{\mbb{R}} \left(-\frac{d}{dx}\delta(x-y) \right)\psi(x)\ dx\\
&= \int_{\mbb{R}} \delta(x-y)\frac{d\psi}{dx}\ dx + [\delta(x-y)\psi(x)]_{-\infty}^{+\infty}\\
&= \frac{d\psi}{dy}
\end{split}
$$

This suggests that $\mbf{D}^\dagger=-\mbf{D}$ thus $\mbf{D}$ is a skew-symmetric operator! Note however that $k$-th derivatives are not all skew-symmetric, namely, consider the identity

$$
\frac{d^k}{dx^k}\delta(x-y) = (-1)^k\frac{d^k}{dy^k}\delta(x-y)
$$

thus $k$-the derivatives are symmetric for $k$ even and skew for $k$ an odd function. This is also obvious to see when we compute 

$$
(\mbf{D}^k)^\dagger = (\mbf{D}^\dagger)^k = (-1)^k \mbf{D}^k 
$$

## Derivative operator from Kernel Function

By the use of taylor series expansion of the kernel function we can find how the kernel function can give rise to a differential operator. We want to show that there exists coefficients $a_k(x)$ such that 


$$
\mbf{A} = \sum_{k=0}^\infty a_k(x)\frac{d^k}{dx^k}
$$


To show this we make a change of variables $y\mapsto x+y$ and taylor expand $\psi(x+y)=\sum_{k=0}^\infty \frac{y^k}{k!} \frac{d^k\psi(x)}{dx^k}$, then applying $\mbf{A}$ on the test function $\psi$ yields

$$
\begin{split}
\mbf{A}\psi = \int_\mathbb{R} A(x,y) \psi(y) \ dy &= \int_{\mathbb{R}} A(x,x+y)\psi(x+y) \ dy = \int_{\mathbb{R}} A(x,x+y)\sum_{k=0}^\infty\frac{y^k}{k!}{\psi^{(k)}(x)} \ dy\\
&= \sum_{k=0}^\infty {\psi^{(k)}(x)}\int_{\mathbb{R}} A(x,x+y)\sum_{k=0}^\infty\frac{y^k}{k!}\ dy = \sum_{k=0}^\infty a_k(x)\frac{d^k\psi}{dx^k}
\end{split}
$$

where we defined $a_k(x)\equiv\int_{\mathbb{R}} A(x,x+y){y^k}/{k!}\ dy$. Of course this assumes that the series converges as $k$ goes to inffinity. 

## The Domain Transformation Operator

Let $F\in\mathcal{G}$ be a functional that transforms the domain of functions, symbolically this is expressed as $\psi(x)\overset{F}{\longmapsto}\psi(f(x))$ where $f:\Omega\rightarrow \Omega$. $F$ is a linear functional, thus it has an associated linear operator $\mbf{F}$ and kernel function $K(x,y)$. The kernel $K$ of $\mbf{F}$ is:

$$
K(x,y) = \delta(y-f(x))
$$

by which is trivial to show that 

$$
\mbf{F}\psi = \int_{\Omega} \delta(y-f(x))\psi(y)\ dy  = \psi(f(x))
$$

To this kinds of linear operators we shall call them <ins>**Domain Transformation Operators**</ins>. 

**Functions that preserve $\mbf{i}$-normalization**

It is interesting to determine what types of transformations $f(x)$ will make the transformed function $\mbf{i}$-normalized. Consider $\psi$ to be an $\mbf{i}$-normalized function, then we want to find all $f's$ such that $\mbf{i}^\dagger\mbf{F}\psi=1$, then 

$$
\int_{\mbb{R}^n} \psi(f(x)) d^nx = \int_{\mbb{R}^n} \psi(y)\det(\underline{f})^{-1}d^ny
$$


where we used the chain rule 

$$
\mbf{d}y = \mbf{d}f(x) = \und{f}\mbf{d}x
$$

with $\mbf{d}$ the infinitesimal differential operator defined as $\mbf{d}f(x) \equiv \infd x\cdot\nabla_x f(x)$, with $\mbf{d}x$ an infinitesimal change. And consequently: $d^ny=\det(\und{f})d^nx$. The above integral suggests that in order to ensure $\mbf{F}$ keeps the $\mbf{i}$-normalization property we have to choose an $f(x)$ such that $\det(\underline{f})=1$.

To ensure $\mbf{i}$-normalization we can also consider satisfying the requirement $\mbf{i^\dagger F}=\mbf{i}$, thus using a change of variables we have 

$$
\int_{\mbb{R}^n} \delta(f(x)-y)d^nx = \int_{\mbb{R}^n} \delta(x-y)\det(\und{f})^{-1} d^nx = \left[\det(\und{f})^{-1}\right](y)
$$

This suggests that we should consider a kernel $F$ of $\mbf{F}$ as 

$$
F(x,y) = \delta(f(x)-y)[\det(\und{f})](y)
$$

which obviously ensures $\mbf{i}$-normalization. Note that we have assumed that $f(x)=y$ only has one solution for each $y$, however if we can find more then one solution, we may express it using the winding number $\mcal{k}(y)\in\mbb{N}$ the number of solutions to the equation $f(x)=y$, then $\mbf{i^\dagger F}$ takes the form 

$$
\int_{\mbb{R}^n} \delta(f(x)-y)d^nx  = \left[\det(\und{f})^{-1}\right]\!(y)\!\ k(y)
$$

as long as $\det(\und{f})\neq 0$ and $k(y)\neq 0$ we can allways redefine $\mbf{F}$ such that $\mbf{F}$ is $\mbf{i}$-normalized. Then the kernel $F$ of $\mbf{F}$ must be of the form 

$$
F(x,y) = \delta(f(x)-y)\frac{[\det(\und{f})](y)}{k(y)}
$$


### A More General Derivative Operator 

Recall that in section [Derivative Operator](# Derivative operator from Kernel Function) we showed how we could provide a differential operator representation of the linear operator by the taylor series of $\psi(x+y)$, however because of convergence properties this change of variables might not be the most appropriate, so instead we consider the change of variable $y\mapsto y+f(x)$. Then the operator will be the composition of a differential operator and a domain transformation operator. Mathematically this is expressed as

$$
\mbf{A}\psi = \sum_{k=0}^\infty a_k(x)\frac{d^k}{dx^k}\psi(f(x)) = \mbf{D}\mbf{F}\psi
$$

where $a_k(x)\equiv\int_{\mathbb{R}} A(x,f(x)+y){y^k}/{k!}\ dy$, the kernel of $\mbf{D}$ and of $\mbf{F}$ are


$$
\begin{align}
\mbf{D} &=  \sum_{k=0}^\infty a_k(x)\frac{d^k}{dx^k}\\
F(x,y) &= \delta(y-f(x)) 
\end{align}
$$


## The Functional Derivative 

In this section we explore a new concept in functional analysis, the functional derivative. In the functional analysis we want to study functionals $G\in\mcal{G}$ and how they change when their argument, which is a function $$\psi\in\mathcal{F}_{\Omega}$$, change. The functional derivative will help us answer the question: How does a functional change when we change its argument by a tiny amount. The first definition that will help us answer this question is the definition of the directional derivative. If we want to determine how a functional changes when its arguments changes in some particular "_direction_" $\phi$ we compute the directional derivative defined as follows

$$
\phi\cdot\mathcal{D} F(\psi) \equiv \lim_{\varepsilon\rightarrow 0} \frac{F(\psi+\varepsilon\phi)-F(\psi)}{\varepsilon}
$$


The directional functional derivative help us determine changes in some particular direction, but what if we do not care about any particular directions. Indeed the functional derivative $\mcal{D}$ will provide us with such an attribute, the functional derivative expresses the directional derivative in all possible different directions $\phi$. We can, indeed, consider a set of basis functions $\phi_1,\phi_2,\cdots$ and their reciprocal functions $\phi^1,\phi^2,\cdots$, that generate the entire functional space $$\mathcal{F}_{\Omega}$$, and write the derivative as

$$
\mcal{D} = \sum_\mu \phi^\mu\phi_\mu\cdot\mcal{D} 
$$

with $\phi_\mu^\dagger\phi^\nu=\delta_{\mu\nu}$. However we are not interested on a definition through some set of basis functions. As the reader may already be aware, the space of functions is infinite dimensional this means that the quantity of basis functions is infinite. To provide a basis independent definition of $\mcal{D}$ we consider an integral over the space of all possible functions, and write it as 


$$
\fder = \int d\phi\  \phi'\phi\cdot \fder = \int_{\mbb{R}^\infty} d^\infty z\ \phi'(z,x)\phi(z,x)\cdot\fder  \equiv \bsym{\Phi}' \bsym{\Phi}\cdot\fder
$$

where $\phi'$ is the operational inverse of $\phi$, that is 

$$
\bsym{\Phi}'^\dagger\bsym{\Phi} = \bsym{\Phi}^\dagger\bsym{\Phi}'=
\int_{\mbb{R}^\infty} \phi'(z,x)\phi(z,y) \ d^\infty z
= \delta(x-y)$$


We can think about $\phi(z,x)$ as being a parameterized function with parameter $z$ of the argument $x\in\Omega$. That is, for each different fixed value of $z$ we have a different function $\phi(z,x)$. We can also think on $\phi$ as the kernel function of a linear operator $\bsym{\Phi}$, with $\bsym{\Phi}'$ being its inverse with kernel $\phi'$. 

I am not satisfied with the previous two definitions since they either involve integration over an infinite dimensional space or a sum over an infinite number of functions, and when we are adding or integrating infinites amounts of something a lot of questions about convergence start to arise, and to be honest I want to avoid needing to answer such questions, let us leave it for the mathmatical purists to figure it out! So, I will provide another definition for the total functional derivative $\fder$ as the following:

$$
(\fder F)^\dagger{\phi} \equiv \int_{\mbb{R}} (\fder F)(y,x)\phi(x)\ dx = \lim_{\varepsilon\rightarrow 0} \frac{F(\psi+\varepsilon \phi)-f(\psi)}{\varepsilon} = \phi\cdot\fder F
$$

What this definition says is that applying the operator $(\fder f)^\dagger$ to some function $\phi$ yields the directional functional derivative, in the $\phi$ direction, evaluated at $\psi$. An important property of the defined derivatives is the following 

$$
\boxed{
\fder_{\phi} \phi\cdot\fder_\psi f(\psi) = \fder_\psi f(\psi)
}\label{eq:der:of:der:identity}
$$

note that $\phi\cdot\fder_\psi f(\psi)$ is just a linear function on $\phi$ thus the derivative transforms the linear equation $\phi\cdot\fder_\psi f(\psi)$ into the operator $\fder_\psi f(\psi)$.

### The Chain Rule 

We want to prove the general chain rule that is provided in the space of fuctions. In finite dimensional euclidean vector spaces, when we compose two functions and take their derivative we get the chain rule for the derivative. The same is true when we compose functionals, let $F(\psi) =  G(f(\psi))$, then the functional chain rule is:

$$\boxed{
\textbf{Chain rule:}\ \ \fder F = \fder_\psi G(f(\psi))  = \bar{f}(\fder_\phi) G(\phi)
}\label{eq:chain:rule}$$

with $\bar{f}\equiv (\fder f)^\dagger$ and $\phi=f(\psi)$. A particular case of the chain rule arises when we consider $f$ to be a linear functional, that is $\phi=f(\psi)=\mbf{A}\psi$ is 

$$
\fder_\psi G(\mbf{A}\psi) = \mbf{A}^\dagger(\fder_\phi) G(\phi) 
$$

To show that the functional chain rule holds we provide the following proof. We start by considering the directional derivative in the $\phi$ direction then

$$
\phi\cdot\mcal{D} F = \lim_{\varepsilon\rightarrow 0} \frac{F(\psi+\varepsilon\phi)-F(\psi)}{\varepsilon} = \lim_{\varepsilon\rightarrow 0} \frac{G(f(\psi+\varepsilon\phi))-G(f(\psi))}{\varepsilon}
$$

now taylor expand the functional $f$ as $f(\psi+\varepsilon\phi) = f(\psi)+\varepsilon \phi\cdot\mcal{D}f(\psi)+O(\varepsilon^2)$, since $\varepsilon\rightarrow 0$ we may ignore second order terms and then write

$$
\phi\cdot\mcal{D} F =\lim_{\varepsilon\rightarrow 0} \frac{G(f(\psi)+\varepsilon \phi\cdot\mcal{D}f(\psi)) - G(f(\psi))}{\varepsilon} =\lim_{\varepsilon\rightarrow 0} \frac{G(\theta+\varepsilon \rho) - G(\theta)}{\varepsilon}
$$

with $\theta=f(\psi)$ and $\rho=\phi\cdot\mcal{D}f(\psi)$. Notice how the leftmost limit is just the derivative of $G$ with respect to its argument thus:

$$
\lim_{\varepsilon\rightarrow 0} \frac{G(\theta+\varepsilon \rho) - G(\theta)}{\varepsilon} = \rho\cdot\mcal{D}_\theta G(\theta)
$$

But now recall that $\rho$ is itself a directional functional derivative of $f$, which I will write as

$$
\rho = \phi\cdot\mcal{D}f(\psi) = (\fder f)^\dagger\phi
$$

then 

$$
\rho\cdot\fder_\theta = \rho^\dagger\fder_\theta =  \phi^\dagger (\fder f)\fder_\theta = \phi\cdot \bar{f}(\fder_\theta)
$$

where 

$$
\boxed{
\begin{align}
\und{f}(\phi) &= \und{f}\phi = \phi\cdot \fder f \equiv (\fder f)^\dagger\phi\\
\bar{f}(\phi) &= \bar{f}\phi = \phi\cdot (\fder f)^\dagger \equiv (\fder f)\phi
\end{align}
}
$$

which will finally give the chain rule for the directional derivative 

$$
\phi\cdot\mcal{D} F = \phi\cdot \bar{f}(\fder_\theta) G(\theta)
$$

Taking the functional derivative with respect to $\phi$ of the above and using the identity $\eqref{eq:der:of:der:identity}$ will yield the chain rule $\eqref{eq:chain:rule}$.


First we need to understand the derivative $\fder$ operator as a function of $y$. In some sense the product $\fder F$ is an 'outer product' of $\fder(y)$ with $F(\psi,x)$, this implies that what the chain rule says: how is $\fder_\phi$ transformed by the linear operator $\bar{f} = \fder f$. Let us expand the integral operator of $\bar{f}(\fder_\phi)$ then

$$
\bar{f}(\fder_\phi) = \int_{\mbb{R}} (\fder f)(y,z)\fder_\phi(z)\ dz 
$$

The resulting product of $\bar{f}(\fder_\phi)$ with $G(\phi)$ will be the outer product, in other words it will result in a function of $x$ and $y$, that is, $\bar{f}(\fder_\phi)G(\phi) = \bar{f}(\fder_\phi)[y] G(\phi)[x]$. Let us express this product as

$$
\bar{f}(\fder_\phi)G(\phi) = \int_{\Omega}(\fder f)(y,z)(\fder_\phi G)(z,x)\  dz \equiv (\fder f)(\fder_\phi G)
$$

Which shows that the chain rule results in the product of two linear operators or in the product of linear operator with a function when $G$ does not explicitily depend on $x$.  

**Particular Case**

A particular case of the chain rule is when I assume that $G$ is simple on $\phi=f(\psi)$ and thus $\fder_\phi(z)G(\phi,x)=G'(\phi)\delta(x-z)$ then 

$$
\fder F = \int_{\mbb{R}} (\fder f)(y,z)\fder_\phi(z) G(\phi,x)\ dz =  \int_{\mbb{R}} (\fder f)(y,z) G'(\phi)\delta(x-z)\ dz = (\fder f)(y,x) G'(\phi,x)
$$

with $G'=\frac{\partial G}{\partial \phi}$ the classic derivative of $G$ with respect to $\phi$. Or more sucintily 

$$\boxed{\fder F=(\fder f) \frac{\partial G}{\partial \phi}}\label{eq:part:case:simple}$$

**Functional Divergence**

Note how we can also consider the functional divergence which we denote $\fder\cdot F$. The divergence is simply the trace of the gradient of $F$, in other words 

$$
\boxed{
\fder \cdot F = \trace(\fder F)= \int_{\mbb{R}} (\fder F)(x,x)\ dx 
}$$

### Simple Functionals and its Derivatives

We call a functional $$F:\mcal{F}_\Omega\mapsto\mcal{F}_\Omega$$ simple when it is a function of its arguments. For each value $x$ of the argument of the function $\psi$ it outputs a value which is a function of $\psi(x)$. Using some non-standard notation a simple functional $f$ satisfies:

$$\phi(x) = f(\psi(x)) = f(\psi)[x]$$

This means that evaluating $\phi = f(\psi)$ at the point $y$ gives us $f(\psi(y))$. If we relate this with functions defined on finite dimensional euclidean space $\mcal{E}^n$, we say that a function $\sigma:\mbb{R}^n\mapsto\mbb{R}^n$ is simple if it acts on vectors element wise. For instance $\sigma\left(\sum_{k=1}^n v_ne_n\right)=\sigma_0(v_n)e_n$. With $e_n$ the standard orthogonal basis vectors of $\mcal{E}^n$. Or with the use of vector notation $\sigma([v_1 \ v_2 \ \dots \ v_n])=[\sigma_0(v_1) \ \sigma_0(v_2) \ \dots \ \sigma_0(v_n)]$. 


The functional derivative of a simple functional $f$ results in a diagonal kernel function

$$
\fder f = \frac{\partial f}{\partial\psi}\delta(x-y)\label{eq:der:simp:func}
$$

To show that this holds we use the definiton of the directional derivative to show that 

$$
\begin{split}
\phi\cdot \fder f &= \lim_{\varepsilon\rightarrow 0} \frac{f(\psi+\varepsilon\phi)-f(\psi)}{\varepsilon}\\
&= \phi\lim_{\varepsilon\rightarrow 0} \frac{f(\psi+\varepsilon\phi)-f(\psi)}{\phi\varepsilon} \\
&= \phi \lim_{\delta\rightarrow 0} \frac{f(\psi+\delta)-f(\psi)}{\delta}\\
&= \phi\frac{\partial f}{\partial \psi}
\end{split}
$$

since at the limit when $\varepsilon$ goes to zero the ratio will not depend on $\phi$ or any of its derivatives, in other words we can take a change of variables $\delta=\varepsilon\phi$. The functional derivative of a functional $F\in\mcal{G}$ is the kernel that transforms $\phi$ in the linear operational equation $\phi\cdot\fder F$. Hence the kernel function of this linear transformation is a diagonal kernel and given by $\eqref{eq:der:simp:func}$. 






### Deriving Euler-Poisson From the Chain Rule

In this section we want to use some of the "mathematical tools" that we built to determine equations of variational calculus. In particular, of most importance we have the Euler-Lagrange equation and the Euler-Poisson equation. What this equations tells us is that their solution is the trajectory of a particle that is a stationary trajectory of the action integral. In other words solutions to the Euler-Lagrange equation aims to find stationary points of 

$$
J = \int_{\mbb{R}} L\left(\psi,\frac{d\psi}{dx}\right)\ dx\label{eq:act:int:eul:lag}
$$

with $L$ some function of both of its arguments $L:\mbb{R}\times \mbb{R}\rightarrow \mbb{R}$ (Note however that $L$ is a simple functional, that is it is a function in the classic sense). The Euler-Lagrange equations tells us that in order to find the stationary paths $\psi(x)$ for the action integral we must determine $\psi$ that solves the Euler-Lagrange equation

$$
\frac{dL}{d\psi} - \frac{d}{dt}\frac{dL}{d\dot{\psi}} = 0
$$

with $\dot{\psi}=\frac{d\psi}{dx}$. With the notation that we have introduced so far we can express the action integral as 

$$
J = \mbf{i}^\dagger L(\psi,\mbf{D}\psi)
$$

with $\mbf{D}=\frac{d}{dx}$. To derive the Euler Poisson equation we recall that the kernel of $\mbb{D}$ is $-\delta'(x-y)$, then using the chain rule and considering $L$ to be simple we have  

$$
\fder^\dagger J = (\fder L)^\dagger \mbf{i} = \left(\left(\frac{\partial L}{\partial \psi} + \mbf{D}^\dagger\frac{\partial L}{\partial \dot{\psi}}\right)\delta(x-y)\right)^\dagger \mbf{i} = \frac{\partial L}{\partial \psi} + \left(\frac{\partial L}{\partial \dot{\psi}}\mbf{D}^\dagger \delta(x-y)\right)^\dagger \mbf{i} 
$$


$\mbf{D}^\dagger$ is a linear operator that transforms elements with argument $x$. The second term is equal to the integral 

$$
\int_{\mbb{R}} \frac{\partial L}{\partial \dot{\psi}} \delta'(x-y)\ dy = -\frac{d}{dx}\frac{\partial L}{\partial \dot{\psi}}
$$

where we used $\mbf{D}^\dagger\delta(x-y)=\delta'(x-y)$. In a similar manner we can derive the Euler-Poisson equations, the key difference is that we will have a function $L$ of multiple arguments where each argument is a $k$ derivative of $\psi$, consequently a linear transformation (in the sense of functions) of $\psi$.  


## Fokker Plank equation

We want to give some sense to the Fokker Plank equation in terms of an evolving probability distribution and how it can relate with the operational equation $\dot{\psi}=\mbf{B}(t)\psi$. While the Fokker Plan is also a time evolving process which can be expressed as an operational equation $\dot{\psi}=\mbf{D}(t)\psi$ it is less general then the operational equation since it ignores higher order moments of the probability distribution and I am not sure if it satisfies the requisites for being a stochastic operator. 

For a one space dimensional process the focker Planck can be written in the form 

$$
\dot{\psi}(t,x) = \frac{d}{dx}(f(t,x)\psi(t,x)) + \frac{1}{2} \frac{d^2}{dx^2} (g(t,x)^2 \psi(t,x) )
$$

Now note that if we define the operator $\mbf{D}(t)=\frac{d}{dx} f(t,x)+\tfrac{1}{2} \frac{d^2}{dx^2} g(t,x)^2 $, then the above equation can be simply written as 

$$
\dot{\psi} = \mbf{D}(t)\psi  
$$

To really understand if this provides us with a proper generator operator that preserves $\mbf{i}$-normalization we have to check some properties of $\mbf{D}$ by looking at its kernel. These properties are 

$$
\mbf{i^\dagger D} = 0, \ \ D(x,y)\geqslant 0\ \forall x\neq y
$$

The kernel $D(x,y)$ of $\mbf{D}$ can be expressed as 

$$
D(x,y) = -f(t,x)\delta'(x-y) + \frac{df(t,x)}{dx}\delta(x-y) +  \frac{1}{2} \frac{d^2g(t,x)^2
}{dx^2} - \frac{1}{2} \frac{d g(t,x)^2}{dx}\delta'(x-y) + \frac{1}{2}g(t,x)^2\delta''(x-y)
$$

Note that since the dirac-delta and its derivatives are zero except at the point where its arguments are zero, then we see this is a diagonal kernel, which means that the inequality $D(x,y)\geqslant 0\ \forall x\neq y$ is straightforwardly satisfied. The other condition might take a little more work. 

$$
\mbf{i^\dagger D} = \mbf{Di} = \frac{d}{dx}(f(t,x)1) + \frac{1}{2} \frac{d^2}{dx^2} (g(t,x)^2 1 ) = \frac{df(t,x)}{dx} + \frac{1}{2}\frac{d^2g(t,x)^2
}{dx^2}
$$

which must be equal to zero. However using the fundamental theorem of calculus we also see that 

$$
\mbf{i^\dagger D}\psi = \int_{\mbb{R}} \frac{d}{dx}(f(t,x)\psi(t,x)) + \frac{1}{2} \frac{d^2}{dx^2} (g(t,x)^2 \psi(t,x) )\ dx = 0
$$

assuming boundary terms vanish i.e. $f(t,\pm\infty)\psi(t,\pm\infty)= \frac{d}{dx}g(t,\pm\infty)^2 \psi(t,\pm\infty) = 0$. This means that, even though $\mbf{i^\dagger D}$ does not vanish the total probability is preserved preserved since $\mbf{i^\dagger D}\psi=\mbf{i}^\dagger\dot{\psi} = \frac{d}{dt}\mbf{i}^\dagger\psi = 0$, which means that $\mbf{i^\dagger}\psi$ is constant on $t$. 


One of the most used differential equations for the description of space-time evolution of density distributions is the fokker plank equation. We can derive a differential equation by truncating the taylor series of $\rho$ at some integer value. To see how this approximation gives rise to a differential equation we write 

Note that both $A(x,y)$ and $\rho(x)$could be made time $t$ independent without any change being made to the above. Then an equation of the form 

$$
\boxed{
\frac{\partial \rho}{\partial t} = A\rho
}
$$

will become of the form 

$$
\boxed{
\frac{\partial \rho}{\partial t} = \sum_{k=0}^N a_k(x,t)\frac{\partial^k \rho}{\partial x^k}
}
$$

which turns into the Fokker plank equation when we set $N=2$ and $a_0(x,t)=0$.



The integral is equivalent to taking spatial $x$ derivatives of the input function $\rho$ and then multiplying each derivative by some function of $x$ and $t$. 


We now are interested in the other problem: Given some coeficients $a_k(x)$ of the differential equation how do we find $A(x,y)$? The answer to this question can be answered by proposing the following transition function 

$$
\boxed{
A(x,y) = \sum_k \delta^{(k)}(x-y) a_k(x)
}
$$

with $\delta^{(k)}(x-y)=\frac{\partial^{k}\delta}{\partial x^k}$ the $k$-th derivative of the delta-de-dirac distribution function. To see if this indeed is equivalent to the differential equation we first need to recall that 
where we assume that $\psi$ and all its derivatives vanish at $-\infty$ and $+\infty$. Then the application of $A$ to $\rho$ gives 

$$
A\rho = \int_{\mathbb{R}}\sum_k \delta^{(k)}(x-y) c_k(x)\rho(y)\ dy = \sum_k a_k(x) \int_{\mathbb{R}}\delta^{(k)}(x-y)\rho(y)\ dy = \sum_k (-1)^k a_k(x)\frac{\partial^{k}\rho}{\partial x^k}
$$

We can consider an aproximate solution by considering a gaussian distribution function instead of the delta-de-dirac, making it computationaly feasible in terms of integral operations. To show how the relationship between the derivatives of the delta-de-dirac function as an operator yield the signed derivative, first 

## Functional directional and total derivatives 

In this section we understand how we can evaluate derivatives with respect to functions of functions and what does that mean in terms of operations. We start by clearly defining the directional functional derivative in the following manner 

$$
\delta_\phi f = \phi\cdot\delta f = \lim_{\varepsilon\rightarrow 0} \frac{f(\psi+\varepsilon\phi) - f(\psi)}{\varepsilon} = \frac{\partial}{\partial t} f(\psi+t\phi) = \phi f'(\psi)
$$

with $$ f'(\psi)=\frac{\partial f}{\partial \psi} = \lim_{\rightarrow 0} \frac{f(a+b)-f(a)}{b} \\|_{a=\psi} $$.


An orthonormal set $e(x,y)$ is defined as the following integral equation

$$
\boxed{
\int_{\mathbb{R}} e(z,y)e(z,x)\ dz = \delta(x-y)
}
$$

for each $x\in\mathbb{R}$ the function $e(x,y)$ is an orthogonal vector (Each column of $e(x,y)$ it's orthogonal to the other columns). 

The total derivative can be thought as a linear operator, in particular applying it to some function $\phi$ yields 

$$
\phi\cdot\delta f = (\delta f)\cdot \phi = \int_{\mathbb{R}} (\delta f)(y,x) \phi(x)\ dx

$$

The differential operator $\delta$ can be regarded as a vector while the product $\delta f$ can be understood as the outer product $\delta(y)f(x)$, this symbol must not be confused with the dirac-de-delta function. If we consider $f=f(\psi)$ to be a simple function of $\psi$ then the functional gradient reduces to 

$$
\delta f = \frac{\partial f}{\partial \psi} \delta(x-y)
$$

Note that when we say a simple function we mean a function that acts element wise on the components of $\psi$, for example the following function is not simple:

$$
f(\psi) = \int_{\mathbb{R}} \int_{\mathbb{R}} \psi(x)\psi(y) A(x,y)\ dx\ dy = \psi\cdot A\cdot\psi
$$

where the values of $\psi$ of diferent components get multiplied together. The directional derivative of this quadratic equation yields 

$$
\phi\cdot\delta f(\psi) = \psi\cdot A\cdot \phi + \phi\cdot A \cdot \psi = \phi\cdot(A+A^\top)\cdot\psi
$$

now taking the derivative with respect to $\phi$ yields the vector 

$$
\delta f = \delta_\phi \phi\cdot\delta f = (A+A^\top)\cdot \psi
$$

while the second derivarive yields the operator 

$$
\delta^2 f = A+A^\top = A(x,y)+A(y,x)
$$

Thus the hessian operator of a scalar quadratic function results on a linear operator. If $A$ is symmetric then we can recover $f$ from the hessian just by computing $f=\frac{1}{2}\psi\cdot (\delta^2 f)\cdot \psi$.


Let us now consider a linear function 

$$
f(\psi) = A\cdot \psi \equiv \int_{\mathbb{R}} A(y,x) \psi(y)\ dy
$$

Using the definition of the directional derivative we easily find that 

$$
\phi\cdot\delta f = A\cdot \phi
$$

while the functional gradient is just the linear operator $A$, that is 

$$
\delta f = \delta_\phi\phi\cdot\delta f = A
$$

Thus for the cases of linear transformations and quadratic equations we recover exactly the operators that we were expecting. And of course we deduce that for $f$ a linear 'transformation' $f(\psi)=(\delta f)\cdot \psi$.


To see that this result aplied as a linear operator to some function $\phi$ does indeed recover the directional derivative we compute 

$$
(\delta f)\cdot \phi = \int_{\mathbb{R}} \frac{\partial f}{\partial \psi} \delta(x-y) \phi(x)\ dx = \frac{\partial f}{\partial \psi}\phi(x) = \lim_{\varepsilon\rightarrow 0} \frac{f(\psi+\varepsilon\phi) - f(\psi)}{\varepsilon\phi}\phi(x) = \phi\cdot\delta f
$$

With $\frac{\partial f}{\partial \psi} = \lim_{\varepsilon\rightarrow 0} \frac{f(\psi+\varepsilon\phi) - f(\psi)}{\varepsilon\phi}$


### A new attempt at defining this crapity crap

Let us think about the total derivative $\fder$ in a different way. Recall that the gradient vector $\nabla$ can be defined as the sum of directional derivatives, thus

$$
\nabla = \gamma^\mu\partial_\mu = e^\mu e_\mu\cdot\nabla
$$

where $\gamma^\mu$ and $e_\mu$ are arbitrary linearly independent vectors with $e_\mu\cdot e^\nu = \delta_{\mu\nu}$. In a similar fashion we may express the total functional derivative $\fder$.

$$
\fder = \sum_\mu \phi^{\mu}\phi_\mu\cdot \fder 
$$

with $\phi_\mu\cdot \phi^\nu = \int_{\mbb{R}} \phi_\mu(x) \phi^\nu(x)\ dx=\delta_{\mu\nu}$. Note however that the space of functions is infinite dimensional, meaning that the sum will go up to infinity. We can provide another representation by considering an integral over all possible functions 

Another property is that if $f(\psi)$ is of rank-$k$ the  its derivative will be of rank-$(k+1)$, in other words in terms of number of arguments of functions, if $f$ is a function of $k$ arguments, i.e. $f=f(\psi,x_1,x_2,\dots,x_k)$ then $\fder f = (\fder f)(\psi,x_1,\dots,x_{k+1})$  





with $\phi=\mbf{A}\psi$. The derivative of a simple function $f$ (a function that only operates on functions element wise) is a diagonal operator. An example of such a function is $f(\psi)=(\psi(x))^\alpha$ for some $\alpha\in\mbb{R}$. 

$$
\fder f = f'(\psi(x))\delta(x-y)
$$

with $f:\mbb{R}\mapsto \mbb{R}$ and its $f'$ is its derivative. $f'(\alpha) = \lim_{\varepsilon \rightarrow 0} \frac{f(\alpha+\varepsilon)-f(\alpha)}{\varepsilon}$ with $\alpha\in\mbb{R}$.  

However if $\phi\in L(\mbb{R}^n,\mbb{R}^{m})$ then 

$$\boxed{
\fder f = \nabla_\psi f(\psi(x))\delta(x-y)
}$$

where $\nabla$ is the conventional $m$-dimensional gradient operator. 


**Notes on the chain rule**

It is not that much clear what the chain rule expressed in $\eqref{eq:chain:rule}$ can be interpeted both computationaly and philosophically. First we need to understand the derivative $\fder$ operator as a vector aka a function. In some sense the product $\fder F$ is a product of $\fder(y)$ with $F(\psi,x)$, this implies that what the chain rule says is how the $\fder_\phi$ gets transformed by the jacobian transformation $\bar{f} = (\fder f)^\top$. Let us expand the integral operator of $\bar{f}(\fder_\phi)$

$$
\bar{f}(\fder_\phi) = \int_{\mbb{R}} (\fder f)(x,y)\fder_\phi(y)\ dy 
$$




In variational problems of mechanichs and physics we usually enconter action integrals, which are just the integral of the Lagrangian function. Let $F$ be some Lagrangian function, then define 

$$
J(F)= \int_{\mbb{R}} F(\psi,x) a(x)\ dx \equiv \mbf{a}\cdot F
$$

we want to determine the functional derivative of $J$ with respect to $\psi$. In order to do that we will make use of the chain rule, then 

$$
\fder J(F) = \bar{F}(\fder_\phi) J(\phi) = \bar{F}(\fder_\phi)\phi\cdot \mbf{a} = \fder{F}(\mbf{a}) = (\fder F)^\top \mbf{a} = \int_{\mbb{R}} (\fder F)(y,x)a(x)\ dx\label{eq:chain:rule:part:case:int}
$$


#### Gradients as Linear Integral Operators 

The gradient vector $\nabla$ can be expressed as an integral, namely reacalling that $\nabla=e^\mu e_\mu\cdot\nabla$ and taking $e_\mu\cdot\nabla = -\int_{\mbb{R}^n}dy\ e_\mu\cdot\nabla \delta(x-y)$ 

$$
\nabla = \mbf{D} = -\int_{\mbb{R}^n}dy\ e^\mu e_\mu\cdot\nabla \delta(x-y) = -\int_{\mbb{R}^n}dy\ \nabla\delta(x-y)
$$

Let $\psi:\mbb{R}^n\rightarrow\mbb{R}$ and define $F=F(\nabla \psi) = F(\mbf{D}\psi)$ then the derivative of $F$ may be computed using the chain rule, thus

$$
\fder F = \mbf{D}^\top(\nabla_y)F(y) = -\nabla_x\cdot\nabla F(\psi(x))
$$

with $y=\nabla \psi$. With $\nabla$ acting on $F$.

#### Deriving Euler-Poisson

To derive the Euler-Poisson equations from the machinery that we have set down so far we will make use of the chain rule. In particular we will consider $J(F)=F\cdot\mbf{i}=\int_{\mbb{R}}F\ dx$ and $F=G(\mbf{D}\psi)$ with $\mbf{D}$ a differential operator and $G$ simple. Equation $\eqref{eq:part:case:simple}$ together with $\eqref{eq:chain:rule:part:case:int}$ gives us the equation 

$$
\fder J = \int_{\mbb{R}}(\fder F)(y,x)\ dx = \int_{\mbb{R}}D(y,x)G'(\phi,x) dx
$$

But if we set $D(y,x)=(-1)^k\delta^{(k)}(x-y)$ the above integral reduces to the derivative with respect to $y$, then 

$$
\fder J  = \int_{\mbb{R}}(-1)^k\delta^{(k)}(x-y)G'(\phi,x) dx = \frac{d^k}{dy^k} G'(\phi,y) = \frac{d^k}{dy^k} \frac{d}{d\phi}G(\phi,y)
$$

with $\phi=A\psi$.








# **IS THE FUNCTIONAL DERIVATIVE PROPERLY DEFINED???**


### Important Identities


The generic chain rule is, let $F(\psi)=G(f(\psi))$ 

$$
\fder_\psi G(f(\psi))= \bar{f}(\fder_\phi) G(\phi)
$$

with $\phi=f(\psi)$

$$
\bar{f}(\phi) = \bar{f}\phi = (\mcal{D}f)^{\ddagger}\phi
$$

Now let $G$ be a simple function, then 

$$
\begin{split}
[\bar{f}(\fder_\phi) G(\phi)](x,y) &= \int_{\mbb{R}^n} (\fder f)^\top(z,x) \fder_\phi(z)[G(\phi)](y)d^nz\\ 
&= \int_{\mbb{R}^n} (\fder f)^\top(z,x)[G'(\phi)](y)\delta(z-y)d^nz\\
&= [G'(\phi)](y) \left[(\fder f)^\ddagger\right](x,y)
\end{split}\label{eq:chain-rule:complete:der}
$$


where we used $$\fder_\phi(z)[G(\phi)](y)=[G'(\phi)](y)\delta(z-y)$$. Another important special case is the derivative of the function $J(\psi)=\mbf{i}^\dagger G(f(\psi))=\mbf{i}^\dagger F(\psi)$, to evaluate this special case first note that 

$$
\fder^\ddagger J = \fder^\ddagger \mbf{i}^\ddagger F = (\fder F)^\ddagger \mbf{i}
$$

Now use $\eqref{eq:chain-rule:complete:der}$ to write 

$$
\fder^\ddagger J = \int_{\mbb{R}}[G'(\phi)](y) \left[(\fder f)^\ddagger\right](x,y)\ dx =   (\fder f)^\ddagger G'(\phi) = \left(G'(\phi)^\ddagger (\fder f) \right)^\ddagger
$$

Finally we arrive at the required result:

$$
\boxed{
\fder J = \fder\!\ \mbf{i}^\ddagger G(f(\psi))= G'(\phi)^\ddagger (\fder f) 
}\label{eq:fder:J:simple}
$$

where the last expression is to be understood as the linear operator $(\fder f)^\ddagger$ applied to the function $G'(\phi)$.


**A trivial proof**

Let $\theta:\mbb{R}^n\mapsto\mbb{R}^n$ and $A\in\mbb{R}^{n\times n}$ some arbitrary matrix function, then 

$$
\und{\theta}\cdot \mathsf{A} = \matrixtrace(\underline{\theta}\mathsf{A}) = \matrixtrace(\theta\nabla^\top\mathsf{A}) = (\mathsf{A}^\top\nabla)^\top \theta = (\mathsf{A}^\top\nabla) \cdot\theta\label{eq:matrix:differential:trace} 
$$

**Gradients** 




Let us now consider the case when we have gradients $\nabla$ of functions involved when computing derivatives. We want to show that the Euler-Lagrange equations only involving $\und{\psi}$ (the jacobian matrix of $\psi$) can be written as 

$$
\fder J = - \nabla_x^\top\nabla_{\und{\psi}} \mcal{L}(\und{\psi}) = 0
$$

Let 

$$
J(\psi) = \int_{\mbb{R}^n} \mcal{L}(\und{\psi}) d^nx
$$

Next compute the directional derivative of $J$ in the $\phi$ direction 

$$
\phi\cdot\fder J(\psi) = \lim_{\varepsilon\rightarrow 0}  \int_{\mbb{R}^n} \frac{\mcal{L}(\und{\psi}+\varepsilon\und{\phi})-\mcal{L}(\und{\psi})}{\varepsilon} d^nx
$$

which simplifies into 

$$
\phi\cdot\fder J(\psi)  =   \int_{\mbb{R}^n} \und{\phi}\cdot\nabla_{\und{\psi}} \mcal{L}(\und{\psi})\ d^nx
$$

Now use $\eqref{eq:matrix:differential:trace}$ and note that $\nabla_{\und{\psi}}$ is a matrix of derivatives, then write:

$$
\und{\phi}\cdot\nabla_{\und{\psi}} = (\nabla_\und{\psi}^\top\nabla_x)\cdot\phi
$$

Now express $\phi\cdot\fder J(\psi) =F(\phi)=\mbf{i}^\dagger G\left(\left(\nabla_\und{\psi}^\top\nabla_x\right)\cdot\phi\right)=\mbf{i}^\dagger G(\mbf{A}\cdot\phi)$ as a function of $\phi$. With $G\left(\left(\nabla_\und{\psi}^\top\nabla_x\right)\cdot\phi\right)=\und{\phi}\cdot\nabla_{\und{\psi}} \mcal{L}(\und{\psi})$ a linear function of the argument $\und{\phi}\cdot\nabla_{\und{\psi}}$. 

Then use $\eqref{eq:fder:J:simple}$ to write 

$$
\fder_\phi \phi\cdot\fder J(\psi) = \fder_\phi F(\phi) = G'(\theta)^\ddagger \mbf{A} = \left(-\nabla_y^\top\nabla_{\und{\psi}} \mcal{L}(\und{\psi})\right)^\ddagger
$$

where we used $G'=\mcal{L}(\und{\psi})$. And also note that $\mbf{A}^\ddagger = (\nabla_y^\top\nabla_{\und{\psi}})^\ddagger = -\nabla_{\und{\psi}}^\top \nabla_x$, when $\nabla_{\und{\psi}}^\top \nabla_x$ interpreted as a linear operator. ($\nabla_x^\dagger = \nabla_y$).

#### Multivariate Euler-Poisson Equation

Now consider the case when $\mcal{L}$ depends on both $\psi$ and $\und{\psi}$ then it can be shown (left as an exorcism for the reader) that 

$$
\fder_\psi J(\psi) = \nabla_{\psi}\mcal{L}-\nabla_x^\top\nabla_{\und{\psi}}\mcal{L} \label{eq:mult:euler:poisson}
$$

where 

$$
J(\psi) = \int_{\mbb{R}^n} \mcal{L}(\psi,\und{\psi})d^nx
$$





## Manifold Projection Operators 

Let us consider operators $\mbf{P}$ that project to $m$-dimensional manifolds of $\mbb{R}^n$.Manifold projection operators take some arbitrary function $\psi$ with domain in $\mbb{R}^n$ and project that function into $\mcal{M}$, in particular, for $\mbf{P}$ an $\mcal{M}$-projection operator we want 

$$
\mbf{P}\psi = 
\begin{cases}
\psi(x)&\text{if}\ x\in\mcal{M}\\
0&\text{otherwise}
\end{cases}
$$

Or in other words, for all points $x$ such that $x\notin \mcal{M}$ we have $$[\mbf{P}\psi](x)=0$$. We can express manifold projection operators by defining 

$$
\mbf{P}\psi = \int_{\mcal{M}} \delta(x-y)\psi(y)d^my 
$$

The question on wether $\mbf{P}$ is a projection operator or not can be answered by determining the transmutation of $\mbf{P}$, and by computing $\mbf{P}^2$. 

To determine the transmutation of $\mbf{P}$ we start by defining $\phi=\mbf{P}\psi$ and then computing 

$$
\omega^\dagger \mbf{P}\psi = \omega^\dagger \phi = \int_{\mbb{R}^n} \omega(x)\phi(x) d^nx = \int_{\mbb{R}^n}\int_{\mcal{M}} \delta(x-y)\psi(y)\ d^my\ d^nx =  \int_{\mcal{M}} \omega(y)\psi(y)\ d^my
$$

where we changed the order of integration. Now note that the last integral can also be written as 

$$
\int_{\mcal{M}} \omega(y)\psi(y)\ d^my = \int_{\mbb{R}^n}\int_{\mcal{M}}\delta(x-y)\omega(y)\ d^my\ \psi(x)\ d^nx = \psi^\dagger\mbf{P}\omega 
$$

Thus $\psi^\dagger \mbf{P}\omega = \omega^\dagger \mbf{P}\psi$ for any functions $\psi,\omega\in\mcal{F}(\mbb{R}^n)$ imply that 

$$
\mbf{P}^\dagger=\mbf{P}
$$

thus manifold projection operators are symmetric operators! To show that $\mbf{P}$ is an indempotent operator we compute 

$$
\mbf{P}^2\psi = \mbf{P^\dagger P}\psi = \int_{\mcal{M}} \delta(x-z)\int_{\mcal{M}} \delta(z-y)\psi(y)\ d^my \ d^mz = \int_{\mcal{M}} \delta(x-z)\psi(z) \int_{\mcal{M}}  \delta(z-y)\ d^my \ d^mz
$$

Now notice that 

$$\boxed{
\int_{\mcal{M}} \delta(z-y)\ d^my = 1, \ \ \forall z\in\mcal{M}
}$$

which means that the above integral will reduce to 

$$
\mbf{P}^2\psi = \int_{\mcal{M}} \delta(x-z)\psi(z)d^mz = \mbf{P}\psi
$$

thus proving that $\mbf{P}^2=\mbf{P}$ is an indempotent operator! By which we conclude that $\mbf{P}$ is a projection operator. Manifold projection operators can be expressed with respect to some parameterization. Consider $\mcal{M}$ is given by a parameterization $$\mcal{M}=\{y\in\mbb{R}^n \ \\| \ y=\sigma(z), \forall z\in\mbb{R}^m\}$$ with $\sigma:\mbb{R}^m\mapsto \mbb{R}^n$ such a parameterization of $\mcal{M}$, then we can show that 

$$
\mbf{P}\psi = \int_{\mbb{R}^m} \delta(x-\sigma(y))\psi(x)\det(\und{\sigma}) d^my  \label{eq:proj:manif:param}
$$

with $\und{\sigma}$ the jacobian matrix of $\sigma$. Note that we have written $\det(\und{\sigma})$, however $\und{\sigma}$ is not a square matrix, this is why in [Change of Variables](#chage-of-variables) we have extended the notion of determinant to non-square matrices.

One of the applications for which I aim to used manifold projection operators is when determining signal that belong to some $m$-dimensional manifold. Think for instance that we have some dataset that expresses the digitalization of some process. Assume that the dataset does not have any type of structure. Now assume that the real process is defined only on an $m$-dimensional manifold, this means that the reconstructed process must be projected to the manifold. A more concrete example is the problem of 3D image reconstruction when the dataset are point clouds, i.e. discrete measuraments of the 3D space. By the nature of the sensures used we know that the dataset must come from discrete samples of a 2-dimensional manifold. However we do not know what is that manifold. Which means that we need to somehow provide interesting algorithms for which we can determine $\mcal{M}$.  

### A Problem on Manifold Projections

Let us consider the problem of determining a manifold $\mcal{M}$ projection. Consider two signals $\phi,\psi\in\mcal{F}(\mbb{R}^n\mapsto \mbb{R})$ then assume that 

$$
\phi \approx \mbf{P}\psi
$$

holds approximatly for some $\mbf{P}\in\mcal{P}$ ($\mcal{P}$ is the space of manifold projection operators). The goal is to find the manifold projection operator $\mbf{P}$ that makes $\phi$ very close to $\mbf{P}\psi$. Consider the following optimization problem 

$$
\underset{\mbf{P}\in\mcal{P}}{\text{minimize}} \ \|\mbf{P}\psi-\phi\|^2\label{eq:p:opt:proj:manif}
$$

We can pose this problem in terms of some parameterization function $\sigma$, namely recall $\eqref{eq:proj:manif:param}$, now express $\mbf{P}=\mbf{P}(\sigma)$ as a functional of the parameterization, then we can rewrite $\eqref{eq:p:opt:proj:manif}$ as the following optimization problem 

$$
\underset{\sigma\in\mcal{F}{(\mbb{R}^m\rightarrow \mbb{R}^n)}}{\text{minimize}} \ \|\mbf{P}(\sigma)\psi-\phi\|^2
$$

Now write the cost function as the quadratic equation 

$$
\|\mbf{P}\psi - \phi\|^2 = \psi^\dagger\mbf{P}^\dagger\mbf{P} \psi + \phi^\dagger\phi - 2\phi^\dagger \mbf{P}\psi = (\psi^\dagger - 2\phi^\dagger)\mbf{P}\psi  + \phi^\dagger\phi 
$$

where we used $\mbf{P}^\dagger\mbf{P}=\mbf{P}$. Now define $\theta=\psi-2\phi$, then the problem will reduce to 

$$
\underset{\sigma\in\mcal{F}{(\mbb{R}^m\rightarrow \mbb{R}^n)}}{\text{minimize}} J(\sigma)\equiv\theta^\dagger\mbf{P}(\sigma)\psi 
$$

We start by computing the derivative with respect to $\sigma$, then

$$
\fder_{\sigma}^\ddagger J(\sigma) = (\fder_\sigma \mbf{P}(\sigma)\psi)^\ddagger\theta
$$

Now write 

$$
\mbf{P}(\sigma)\psi = \int_{\mbb{R}^n}\mcal{L}(\sigma,\und{\sigma}) d^ny
$$

with $\mcal{L}(\sigma,\und{\sigma})=\delta(x-\sigma(y))\psi(x)\det(\und{\sigma})$. Now use the Euler-Poisson chain rule $\eqref{eq:mult:euler:poisson}$ to write 

$$
\begin{split}
\fder_\sigma \mbf{P}(\sigma)\psi &= \nabla_{\sigma}\mcal{L} - \nabla_y^\top\nabla_{\und{\sigma}}\mcal{L} \\
&= \psi(x)\left(\det(\und{\sigma})\nabla_\sigma^\top \delta(x-\sigma(y))-[\nabla_y^\top\delta(x-\sigma(y))]\nabla_{\und{\sigma}}\det(\und{\sigma}) - \delta(x-\sigma(y))[\nabla_y^\top\nabla_{\und{\sigma}} \det(\und{\sigma})]\right)\\
&= \psi(x)\left(-\det(\und{\sigma})\nabla_x^\top \delta(x-\sigma(y)) + (\bar{\sigma}\nabla_x)^\top \delta(x-\sigma(y)) \nabla_{\und{\sigma}}\det(\und{\sigma}) -\delta(x-\sigma(y))[\nabla_y^\top\nabla_{\und{\sigma}} \det(\und{\sigma})]\right)
\end{split}
$$

Next we compute 

$$
\begin{split}
(\fder_\sigma \mbf{P}(\sigma)\psi)^\ddagger\theta &= \int_{\mbb{R}^n}  \psi(x)\theta(x) \left(-\det(\und{\sigma})\nabla_x^\top \delta(x-\sigma(y)) + (\bar{\sigma}\nabla_x)^\top \delta(x-\sigma(y)) \nabla_{\und{\sigma}}\det(\und{\sigma}) -\delta(x-\sigma(y))[\nabla_y^\top\nabla_{\und{\sigma}} \det(\und{\sigma})]\right)\ d^nx\\
&= -\nabla_{\sigma}^\top (\psi(\sigma)\theta(\sigma)) \left[\und{\sigma} \nabla_{\und{\sigma}} \det(\und{\sigma}) -\mathsf{I}\det(\und{\sigma}) \right] - \psi(\sigma)\theta(\sigma) [\nabla_y^\top\nabla_{\und{\sigma}} \det(\und{\sigma})]
\end{split}
$$



<!-- $$ -->
<!-- \theta^\dagger\mbf{P}(\sigma)\psi  = \int_{\mbb{R}^n} \theta(x)[\mbf{P}(\sigma)\psi](x) \ d^nx -->
<!-- $$ -->


### A note on the change of variable for non-square jacobians
{: #chage-of-variables}

We want to show the definition of the determinant by showing how a parameterization of the manifold gives rise to an infinitesimal volume change. In particular we want to show that 

$$
\begin{split}
|\det(\und{f})| &= |\und{f}(I_m)|=|\und{f}(e_1)\wedge \und{f}(e_2)\wedge\cdots\wedge \und{f}(e_m)|\\
&= |(\partial_1f)\wedge(\partial_2f)\wedge\cdots\wedge (\partial_m f)|  
\end{split}
$$

with $\partial_i=e_i\cdot\nabla$, the derivatives at all possible directions $e_i$, with $e_i$ spanning the entire vector space $\mbb{R}^m$. 

We start with some important definitions. The infinitesimal change operator tells us how a function $f$ changes when we change $x$ by an infinitesimal amount, that is $x+\infd x$. The infinitesimal change operator $\infd$ is defined as

$$\boxed{\infd f(x) \equiv f(x+\infd x)-f(x)= \infd x\cdot\nabla f(x) = \und{f}\infd x}.$$

Next we define the infinitesimal directed volume change operator $\infD$ as 

$$\boxed{\infD^m f(x) \equiv (\infd^1 f(x))\wedge(\infd^2 f(x))\wedge\cdots\wedge (\infd^m f(x)) }\label{eq:def:inf:vol:el}$$

with $\infd^if(x)=\infd x^i\cdot\nabla f(x) = \und{f}\infd x^i$, the change with respect to the $i$-th component of $x$. We call $\infD^m f(x)$ the infinitesimal directed volume element at the point $x$. From the definition of the infinitesimal directed volume element we can use the outermorphism property to compute

$$
\infD^m f(x) = (\und{f}\infd x^1)\wedge (\und{f}\infd x^2)\wedge\cdots\wedge (\und{f}\infd x^m) = \und{f}(\infd x^1\wedge \infd x^2\wedge\cdots\wedge \infd x^m) = \und{f}(\infD^m x)
$$

where of course, from the definition of $\infD$ in $\eqref{eq:def:inf:vol:el}$ we have  

$$
\infD^m x = \infd x^1\wedge \infd x^2\wedge\cdots\wedge \infd x^m
$$

Now note that however we can further express $\infD^m x$ in terms of the infinitesimal volume components $dx_i$ by setting $\infd x^i = e_i dx_i $, then 

$$
\infD^m x = e_1\wedge e_2\wedge \cdots\wedge e_m dx_1dx_2\cdots dx_m = I_m d^mx\label{eq:infD:vol:elem:Im}
$$

with $I_m = e_1\wedge e_2\wedge \cdots\wedge e_m$ the unit pseudoscalar of $\mbb{R}^m$ and $d^mx=dx_1dx_2\cdots dx_m$ the infinitesimal volume element. Now note that we can take the magnitude on both sides of $\eqref{eq:infD:vol:elem:Im}$ to obtain 

$$
|\infD^m x| = |I_m|d^mx = d^mx
$$

This means that also by a change of variables $y=f(x)$ we have

$$
d^m y = |\infD^m y| = |\infD^m f(x)| = |\und{f}(\infD^m x)| = |\und{f}(I_m)|d^mx = \det(\und{f})d^mx
$$

where we set $\|\und{f}(I_m)\|\equiv\det(\und{f})$. 


### Taylor Series expansion 

We can consider taylor expanding any arbitrary function, recall the classic taylor series expansion, and apply it to functions of functions to yield the operational equations

$$
f(\psi + \phi) = \sum_{k=0}^\infty \frac{(\phi\cdot\delta)^kf(\psi)
}{k!} $$

For the case where we have at most quadratic terms this reduces to 

$$
f(\psi + \phi) = f(\psi) + \phi\cdot(\delta f(\psi))+(1/2)\phi\cdot (\delta^2 f(\psi))\cdot \phi
$$

with $\delta f$ a vector and $\delta^2 f$ a matrix when $f$ is a scalar function.


## Euler Poisson Equations

We can derive the Euler-Poisson equations by transforming derivatives into linear integral operators and then using the chain rule while evaluating the derivative. First we consider the following property for the composition of functions, let $F(\psi)=G(f(\psi))$ then 

$$
\delta F = \delta_\psi G(f(\psi))  = \bar{f}\delta_\phi G(\phi)|_{\phi=f(\psi)}
$$

where $\bar{f}$ is a linear integral operator and is the Jacobian transpose of $f$, namely $\bar{f}\equiv (\delta f)^\top$, which on itself is a linear integral operator. A generalization of this chain rule can be attained when we consider $G$ with multiple arguments, specifically $F(\psi)=G(f_1(\psi),\dots,f_n(\psi))$ then 

$$
\delta F = \sum_k \bar{f}_k \delta_{\phi_k} G(\phi_1,\dots,\phi_n)|_{\phi_k=f_k(\psi)}
$$

While for the particular case when each $f_k$ is a linear operator, by which $f_k(\psi)=A_k\psi$ we have the following 

$$
\delta F = \sum_k A_k^\top \delta_{\phi_k} G(\phi_1,\dots,\phi_n)|_{\phi_k=A_k\psi}
$$

Letting the $A_k$ be derivative operators $A_k=\frac{\partial^k}{\partial x^k}=(-1)^k\int_\mathbb{R}dy\ \delta^{(k)}(x-y)$ expressed as linear integral operators. Next use the notation $\psi^{(k)}=A_k\psi=\frac{\partial^k\psi}{\partial x^{k}}$ and recall that for a simple function $f=f(\psi)$ $\delta f = \frac{\partial f}{\partial \psi} \delta(x-y)$, since we assume that $G$ is simple(acts 'element' wise) on each of its arguments then $\delta_{\phi_k}G=\frac{\partial G}{\partial\phi_k}\delta(x-y)$, then the functional derivative becomes 

$$
\delta F = \sum_k A_k  \frac{\partial G}{\partial\phi_k}\delta(x-y) = \sum_k (-1)^k \frac{\partial^{k}}{\partial y^k}\frac{\partial G}{\partial\phi_k}\delta(x-y)
$$

Integrating the above expression over $y$ yields 

$$
\delta J({F}) = \int_{\mathbb{R}} \delta F(\psi(x))\ dy = \int_{\mathbb{R}}\sum_k (-1)^k \frac{\partial^{k}}{\partial y^k}\frac{\partial G}{\partial\phi_k}\delta(x-y)\ dy =  \sum_{k} (-1)^k \frac{\partial^k}{\partial x^k} \frac{\partial F}{\partial \psi^{(k)}}
$$

as required from the Euler-Poisson equations.



We are expected to arrive at the expression 

$$
\delta J({F}) = \sum_{k} (-1)^k \frac{\partial^k}{\partial x^k} \frac{\partial F}{\partial \psi^{(k)}}
$$

where $J(F) = \int_{\mathbb{R}} F(\psi(x))\ dx = i\cdot F(\psi)$, with $i$ a "vector of all ones", i.e $i \cdot = \int_\mathbb{R}\ dx$. First note that we can rewrite the derivative of $J$ as 

$$
\delta J = \delta (i\cdot F) = (\delta F)^\top\cdot i = \int_{\mathbb{R}} \left(\delta F(\psi(x))\right)^\top \ dx
$$


The transpose ${}^{\top}$ changes the arguments of functions for example $A\cdot = \int_{\mathbb{R}}\ dx A(x,y)$ and $A^\top\cdot = \int_{\mathbb{R}}\ dx A(y,x)$, thus $A\rightarrow A(x,y)$ and $A^\top\rightarrow A(y,x)$.



## The operation equation as an eigenvalue equation 

We are interested in solving the following linear operational equation 

$$
\dot{\psi} = A\psi
$$

with initial conditions $\psi(t_0,x)=\psi_0(x)$. In terms of operators we have seen that the straightforward solution is simply 

$$
\psi(t,x) = e^{At}\psi_0(x)
$$

But as we have also seen evaluating $e^{At}$ involves the composition of multiple integrals. Another interesting approach is to determine the eigendecomposition of the operator $A$. Consider that we find operators $U,\Lambda$ such that 

$$
A = U\Lambda U^{-1}
$$

with $UU^{-1}=I=\int_{\mbb{R}}dy\ \delta(x-y)$ with $I$ the identity operator, and $\Lambda=\int_{\mbb{R}}dy\ \lambda(x)\delta(x-y)$ a diagonal operator, that is, it acts element wise. Let $U(x,z)$ label an eigenvector, that is, for each different $z$ the function $U(x,z)$ is the $z$-th eigenvector of $A$. If we write it instead as $u_z(x)=U(z,x)$, then we may write the eigenvalue equation in the form

$$
\boxed{
Au_z = \lambda_z u_z
}
$$

So instead of indexing the eigenpair with integers $(\lambda_i,u_i)$, $i\in\mbb{N}$ we index with a real number $z\in\mbb{R}$, thus the eigenpair $(\lambda_z,u_z)$ not only is a function of $x$ but also of its indexing $z$. Assuming that we found the eigendecomposition of the operator $A$ then the exponential is trivial to evaluate, namely

$$
e^{At} = \sum_{k=0}^\infty {(U\Lambda U^{-1}t)^k}{k!} = I + tU\Lambda U^{-1}+(t^2/2)U\Lambda U^{-1}U\Lambda U^{-1} + \cdots = I  + tU\Lambda U^{-1}+(t^2/2)U\Lambda^2U^{-1} + \cdots = Ue^{\Lambda t}U^{-1}
$$

The operator $e^{At}$ shares the same eigenvectors of $A$ but has eigenvalues $e^{\Lambda t}$, the eigenpair of $e^{At}$ is $(e^{\lambda_z t},u_z)$. This can also be stated as the eigenvalue equation 

$$
e^{At}u_z = e^{\lambda_z t}u_z
$$

Then to determine the general solution given the initial conditions, we must express $\psi_0(x)$ as a 'linear combination' of the eigenvectors, for that we assume the existence of $a_0=a_0(x)$ such that 

$$
\psi_0 = Ua_0 = \int_{\mbb{R}}U(x,y)a_0(y)\ dy
$$

then 

$$
\psi = e^{At}\psi_0 = Ue^{\Lambda t}U^{-1}Ua_0 = Ue^{\Lambda t} a_0 = U(e^{\lambda(x)t}a_0(x))
$$

The operator $U$ is just a change of basis operator! Defining $\phi=U^{-1}\psi$ we find 

$$
\phi(x,t) = e^{\lambda(x)t}a_0(x)
$$

In the different basis we find that the solution $\phi=\phi(x,t)$ is very simple. 


## Estimating the transition function from samples



Let us consider samples $x(k)$ from the probability distribution at time $k$, and that the probability transition function does not depend on time $k$. An approach that would yield almost good results would be to define the count function

$$
C(x,y) = \sum_{k=0}^{N-1}g(x(k)-x,x(k+1)-y)
$$

with $g$ an appropriate function.



## A note on Functional Neural Networks

One of the goals of this formalism is to enable the realization of functional neural networks. A functional neural network is indeed a functional $F\in\mcal{G}$ that transforms input functions of $\mcal{F}(\Omega)$ into output functions $\mcal{F}(\Theta)$, formally this can be expressed as 

$$
F:\mcal{F}(\Omega)\mapsto\mcal{F}(\Theta)
$$

while in classical neural networks the weights are matrices or vector of $\mbb{R}^n$ in functional neural networks weights are linear operators or functions. By composing linear functional layers with non-linear layers we can provide a complex/complicated neural network that will express some generalized functional. 

An example application is the problem of classification. The classification problem aims to find to which class or cathegory does an object belong to. Usually the input is just a set of matrices or vectors that represent the object, while the output is discrete probability that tells the probability of that object being in that class. We want a more geometrically meaningfull approach, so, instead of encoding the object data as ordered arrays of floating points we express it as a function, this could be a density function in space and time, a colored function which tells at each point $x\in\mbb{R}^3$ what the color is.  The input is a function that maps $\mbb{R}^n$ into $\mbb{R}^m$. If I have RGB data in 3D space and time then the input function would take space time coordinates $x\in\mbb{R}^4$ as their input and as output a color coordinate $\mbb{R}^3$. The output function of the neural network will then provide a probability density distribution, it can be a discrete probability density for each point in space-time ${x}\in\mbb{R}^4$ that tells us what is the probability of the point to belong to that object. Another possibility is to consider as output a continuous distribution function for which we assume a continuous cathegorization of data (consider continuous labels). Or even a continuous probability density for each point $\rho(x,y)$ where $x$ is the point in space-time and $y$ is the point in cathegorical space.

In summary, input data is described as a function of space and time, and the output describes a cathegorical probability associated with that input data.  

In order to train functional neural networks we need a new set of tools that provides with a way to update weights that will decrease the value of the cost function. This is one of the reasons why we introduced the functional derivative and the directional functional derivative, so that we can think about stationary points and local/global minimum with respect to the weight functions. A simple example of a functional neural network can be expressed as 

$$
\psi_{k+1} = \sigma\left(\mbf{W}_k\psi_{k} + \mbf{c}_k \right)
$$

with $\sigma$ a simple functional of its argument, $\mbf{W}_k$ are linear functionals and $\mbf{c}_k$ are trainable functions and $\psi_0$ is the input. $k$ is the neural-network layer. Then the goal is to minimize the cost functional $J(\mbf{W}_0,\dots,\mbf{W}_N)$ with respect to the weight functionals $\mbf{W}_0,\dots,\mbf{W}_N$ for all input functions that represent different datasets. 


Functional neural networks naturaly provides with point permutation invariance, since data is codified in a manner that does not order objects with coordinates (or does not order coordinates).  

### Minimizing Cost Functionals

To be able to optimize cost functionals we will use the functional taylor series expansion. This can then be used to approximate a functional around a function. At some point in time me or someone else will prove that any functional $F$ can be expressed as the taylor series

$$
F(\psi+\phi) = \sum_{k=0}^\infty \frac{(\phi\cdot\fder)^kF(\psi)}{k!}
$$

Newton's method and gradient descent can be established by approximating the taylor series to second and first order respectively. For gradient descent we assume the approximation 

$$
F(\psi+\phi)\approx F(\psi) + \phi\cdot\fder F(\psi) = F(\psi) + \phi^\dagger \fder F(\psi)
$$

The direction of greatest descent is $\phi = -\eta\fder F$ with $\eta>0$. Then the gradient descent update is 

$$\boxed{
\psi_{k+1} = \psi_k - \eta_k \fder F(\psi_k)
}$$

Newton's method considers a second order approximation then 

$$
F(\psi+\phi)\approx  F(\psi) + \phi^\dagger \fder F(\psi)+\tfrac{1}{2}\phi^\dagger \fder^2 F(\psi)\phi
$$

where $\fder^2 F$ is the hessian operator of $F$ (equivalent to the hessian matrix). Minimizing the approximate functional yields newtown's method

$$
\boxed{
\psi_{k+1} = \psi_k - \alpha_k (\fder^2 F(\psi_k))^{-1}\fder F(\psi_k)
}
$$

while Newtons' method involves computing the inverse of the linear operator $\fder^2 F$ gradient descent only involves the computation of the functional gradient of $F$. 

### Equivariant and Invariant Neural Networks


$\newcommand{\RR}{\msf{R}}$

The other goal that we also aim to achieve is Neural Networks with equivariance and invariance properties. To achieve such goal we start by defining two types of functional, the covariant functional, which transforms rotations of the domain of functions into the counter-domain, this can be expressed as 

$$
F(\psi(\mathsf{R}x)) = \mathsf{R}^\top F(\psi(x))
$$

with $\mathsf{R}\in SO(p,q)$ a special pseudo orthogonal transformation matrix. A very trivial functional with such a property is 

$$
F(\psi) = \int_{\mbb{R}^n} x \rho(\|x\|,y) \psi(x) dx 
$$

with $\rho(\|x\|,y)$ some arbitrary function of its arguments. The functional $F:\mcal{F}(\Omega)\mapsto\mcal{F}(\Theta)$. Note however that we may also consider non-linear functionals, that transform $\psi$ non-linearly. The other functional that we are interested in is the invariant functional $g:\mcal{F}(\Theta)\mapsto\mcal{F}(\Theta)$. Consider $\Psi\in\mcal{F}(\Theta)$ then $g$ must satisfy

$$
g(\RR\Psi(y)) = g(\Psi(y))
$$


A trivial example of such a functional is 

$$
g(\Psi) = \|\Psi(y)\|^2_{\mcal{V}}[y]
$$

where the norm is to be taken element wise, that is, the norm in the vector space $\mcal{V}$. This functional is a simple functional and satisfies the equivariance property since for any $\Phi$ and $\Psi$ we have

$$\RR(\Phi)\cdot\RR(\Psi) = \Phi\cdot\Psi$$

where $\cdot$ is the inner product on $\mathcal{V}$. However we can construct other functionals, another trivial one is 

$$
g(\Psi) = \Psi(x)\cdot\Psi(y)
$$

This results on a kernel of a linear operator and can be interpeted as the inner product operator with respect to the function $\Psi$. 

## Estimating Rotations from covariant transformations 

Let us start by defining an arbitrary function $\Psi(y):\mbb{R}\rightarrow \mcal{V}$ which is a functional of $\psi$, in other words let $F\in\mcal{F}(\mcal{V})$ some covariant functional, then define

$$
\Psi \equiv F(\psi)
$$

Now consider a function $\phi$ that relates with $\psi$ via a change of coordinates i.e. :

$$
\phi(x) = \psi(\RR x)
$$

then because of the equivariance property of the covariant functional $F$ we find that 

$$
\Phi \equiv F(\phi) = F(\psi(\RR x)) = \RR^\top F(\psi) = \RR^\top \Psi
$$

this sugests that we may be able to compute $\RR$ by minimizing the cost function 

$$
J(\RR) = \| \RR\Phi - \Psi\|^2
$$

subject to the constraint $\RR^\top\RR = \msf{I}$ and with $\\|\cdot\\|^2$ the functional norm given by 

$$
\|A\|^2 = \int_{\mbb{R}} |A(x)|^2 dx
$$

where $\|\cdot\|^2 $ is the norm on the field $\mcal{V}$. To provide a way to weigh more some points we can also provide a weighing operator $\mbf{W}$ with scalar kernel, this way we can rewrite the cost function in a weight dependent manner as

$$
J(\RR) = \|\mbf{W}(\RR\Phi-\Psi)\|^2
$$

where the kernel of $\mbf{W}$ is $W:\mbb{R}\times \mbb{R}\mapsto \mbb{R}$ a scalar function with two scalar arguments.


## Functional Operators on Fields


Now we generalize the ideia of functional to the idea of functional operators on some arbitrary field $\mcal{V}$. 




To distinguish between vector spaces and function spaces, we will call vectors to elements from some arbitrary finite dimensional pseudo-Euclidean space $\mcal{V}$. $\mcal{V}$ is so general that it can represent any associative hypercomplex system or just some arbitrary vector space. $\mcal{V}$ is a field, a geometric algebra $\mbb{G}$ can be interpeted as a field. Transformations of fields can be achieved by the use of transformation functions $\msf{F}\in\mcal{F}(\mcal{V}\mapsto\mcal{V})=\\{\msf{G} : \mcal{V}\mapsto\mcal{V}\\}$. 



  - Scalar functions of fields: $\psi,\theta,\phi\in\mcal{F}(\mcal{V}) = \\{\psi : \mcal{V}\mapsto\mbb{R}\\}$  
  - Functions that transform fields $\msf{F},\msf{G},\msf{H}\in\mcal{H}(\mcal{V})=\\{\msf{F}:\mcal{V}\mapsto\mcal{V}\\}$ 
  - Non scalar functions of fields $\Psi,\Theta,\Phi,\Omega\in\mcal{F}(\mcal{V}\mapsto\mcal{V})=\\{\Psi : \mcal{V}\mapsto\mcal{V}\\}$
  - Linear transformations $\msf{A}$ of fields can be denoted $\msf{A}\in\mcal{L}(\mcal{V})\equiv\mcal{L}_{\mcal{V}}$ that transforms $\mcal{V}$ into $\mcal{V}$ linearly
  - While ${}^\dagger$ swaps the order of the arguments of functions, the operator ${}^\top$ swaps the order of indexation (This means that $$[\mbf{A}^\top]_{ij}=[\mbf{A}]_{ji}$$ while $$[\mbf{A}^\dagger](x,y)=[\mbf{A}](y,x)$$).
  - The degree of a function tells us how many arguments the function has and rank means how many indices it has, for example a rank-2 degree-3 function $A$ is $A_{ij}(x,y,z)$ with $x,y,z\in\mcal{V}$.

A field $\mcal{V}$ is a pseudo-euclidean vector space with some important multiplicative properties. When we talk about fields we can also talk about transformations of fields, of upmost importance are the linear transformations of fields. When we talk about a linear transformation we are just considering a 'matrix' that transforms fields $\mcal{V}$, let $x\in\mcal{V}$ then $\msf{A}x\in\mcal{V}$ is a linear transformation of $x$. Of course we can also talk about general transformations that transform fields, let $\msf{A}\in\mcal{F}(\mcal{V}\mapsto\mcal{V})$ then $\msf{A}(x)\in\mcal{V}$. 



To generalize linear operators, we have to consider non-scalar kernel functions, thus in general the kernel of $\mbf{A}$ is either itself a field $\mcal{V}$ or a general linear transformation 'matrix'. Since the prior is a generalization of the previous we shall use $\msf{A}(x,y)$ to denote non-scalar kernels of linear operators. This means that for a function $\Psi(x)$ the operator $\mbf{A}$ will act as the integral 

$$
\mbf{A}\Psi = \int_{\mcal{V}} \msf{A}(x,y)\Psi(y)\ dy
$$

where $\msf{A}(x,y)\in\mcal{L}(\mcal{V})$ or $\msf{A}\in\mcal{K}(\mcal{V})$ meaning the kernel $\msf{A}$ of $\mbf{A}$ is a 'matrix' valued function of two arguments $x,y\in\mcal{V}$. For each value of $x$ and $y$ we are computing inside of the integral the matrix vector product $\msf{A}(x,y)\Psi(y)$, so its a linear transformation inside of a linear operator, it is the inception of linear operators and transformations (a transformation inside a transformation inside a transformation inside a transformation ....).  


# Valuing Options from Stochastic Models



To value an option we need to estimate future stock values, the problem is that stock values are not deterministic, they are stochastic! The value of an option is the average payoff value adjusted to today's value. This means that to determine the option value we need to estimate the probability distribution of the stock at time $t+T$ knowing that the value of the stock today is $s_t$. This can be expressed a linear operator equation 

$$
\dot{\psi}(\tau) = \mbf{A}(\tau)\psi(\tau)
$$

The Peano Baker series ensures that we can compute an operator $\bsym{\Pi}=\bsym{\Pi}(\tau,\tau+T)$ such that 

$$
\phi_t = \psi_t(T) = \bsym{\Pi}(t,t+T)\psi_t
$$

$\phi_t$ is the probability distribution function of a stock at time $t+T$ given the probability distribution $\psi_t$ at today's time $t$. Let $S_t\sim\phi_t$ be the random stock variable distrubuted as the function $\phi_t$. The random payoff variable $P_t$ for a call option is 

$$
P_t = \max(0,S_t-K)
$$

The expected payoff can be determined using the definition of expected value:

$$
\mbb{E}[P_t] = \int_{\mbb{R}}\max(0,s-K)\phi_t(s)\ ds = \int_{K}^\infty (s-K)\phi_t(s)\ ds
$$

As important as the expected payoff we also care about the standard deviation of $P_t$, which can be computed as

$$
\text{Var}[P_t] = \mbb{E}\left[(P_t-\mbb{E}[P_t])^2\right] =   \mbb{E}[P_t^2] - \mbb{E}[P_t]^2 = \int_{K}^\infty (s-K)^2\phi_t(s)\ ds - \left(\int_{K}^\infty (s-K)\phi_t(s)\ ds\right)^2

$$

The value $v$ of an option is an explicit function of both the strike value $K$ and time $t$. Let us define it clearly as


$$
v(t,K) = e^{-r(t)T}\mbb{E}[P_t] = e^{-r(t)T}\int_{K}^\infty (s-K)\phi_t(s)\ ds
$$

while $\phi_t$ is a functional of $\psi_t$. In other words assuming $\bsym{\Pi}(t,t+T)\equiv\bsym{\Pi}_{t}^T$ a known linear operator we write $v$ as also a function of $\psi_t$:

$$
v(t,K,\psi_t) = e^{-r(t)T}\int_{K}^\infty (s-K)\left[\bsym{\Pi}_{t}^T\psi_t\right](s)\ ds 
$$

Note that of particular importance we care about a known initial value of the stock price, this means that the distribution function $\psi_t$ simplifies to 

$$
\psi_t(s) = \delta(s-s_t)
$$

where $s_t$ is the know stock value at time $t$, this suggests that $v$ is in fact also a function of the today's stock value. 



## Shift Operators and Generators

Let us consider the following operator equation 

$$
\dot{\psi} = d\mbf{H}\psi
$$

with $\mbf{H}$ a constant operator of time. The solution to this equation is trivial and given by


$$
\psi(t) = e^{td\mbf{H}}\psi_0
$$

where $\psi_0$ is the initial distribution function. We care about shift operators $\mbf{S}$ that shift the domain of functions i.e, 

$$
\psi(x) \overset{\mbf{S}}{\longmapsto}\psi(x+y)
$$

Interestingly if we taylor expand $\psi(x+y)$ around $\mbf{y}$ we find that 

$$
\psi(x+y) = \sum_{k=0}^\infty \frac{y^k}{k!}\frac{d^k}{dx^k}\psi(x)
$$

Note that if we define $\mbf{H}\equiv \frac{d}{dx}$ the above reduces to 

$$
\psi(x+y) = \sum_{k=0}^\infty \frac{y^k}{k!}\mbf{H}^k\psi = e^{y\mbf{H}}\psi
$$

which means that $\mbf{H}$ is the generator of shift/translations of the domain of functions. $\mbf{H}$ generates the lie-algebra of shift operators. Another simple way to express the shift operator is through its kernel, in simple notation we may write that 

$$
\mbf{S}\psi = \int_{\mbb{R}}\delta(x+y-z)\psi(z)\ dz
$$

thus the kernel $S$ of $\mbf{S}$ is the function 

$$
S(x,z) = \delta(x+y-z)
$$

The key ideia between a discrete and continuous time models is that the continuous solution generates the lie-algebra for the discrete solution. 

In a discrete description we would write a shift operator equation in the form

$$
\psi_{k} = \mbf{S}(k)\psi_{k+1}
$$

with the kernel of $\mbf{S}(t)$ given by $S(t,x,y)= \delta(x-y+d_k\Delta t)$. While the continuous time operator equation is 

$$
\dot{\psi}(t) = d(t)\mbf{H}\psi(t)
$$

where we can relate the  shifts in discrete and continuous time domains by

$$ d_k \Delta t= \int_{t}^{t+\Delta t} d(t)\ dt$$


This suggests that we may introduce a dividend shift in the continuous time equation by pre or post multiplying by $d(t)\mbf{H}$ then the differential equation takes the form

$$
\dot{\psi}(t) = d(t)\mbf{H}\mbf{A}(t)\psi(t)
$$


### Space Dependent Shifts

We May also consider shifts that also depend on space $x$, in particular the Kernel function $S$ of a shift operator $\mbf{S}$ that shifts by a $\mu(x)$ amount is 

$$
S(x,y) = \delta(x+\mu(x)-y)
$$

Similarly to what we did previously we can taylor expand $\psi(x+\mu(x))$ thus

$$
e^{\mu(x)\mbf{H}}\psi = \psi(x+\mu(x)) = \sum_{k=0}^\infty \frac{\mu(x)^k}{k!} \frac{d^k}{dx^k}\psi(x)
$$

thus $\mbf{H}=\frac{d}{dx}$. And this means that the differential equation which dictates a space shift $\mu(x)$ is 

$$
\dot{\psi} = \mu(x)\mbf{H}\psi
$$

We can also introduce shift that vary through time by making $\mu$ and $x$ and $t$ dependent, then

$$
\dot{\psi} = \mu(x,t)\mbf{H}\psi
$$


To check if the operator $\mu(x,t)\mbf{H}$ preserves $\mbf{i}$-normalization we compute 

$$
\mu\mbf{Hi} = \int_{\mbb{R}} \mu(x,t)\delta'(x-y)\ dx = -\left[ \frac{d\mu}{dx}\right]_{-\infty}^{+\infty}
$$

which must vanish at the endpoints $\pm\infty$. Thus $\frac{d\mu}{dx}$  must go to zero at the endpoints for all times $t$.

### Dividends Paying

A reasonable model that includes dividends paying at discrete time intervals can be expressed as the following 

$$
\boxed{
\psi_{k+1} = e^{d_k\mbf{D}}\mbf{B}_k \psi_k
}
$$

recall that the operator for geometric brownian motion $\mbf{K}=\mbf{K}(t)$ is a function of time, this means that the solution to 

$$
\dot{\psi}(t) = \mbf{K}(t)\psi(t)
$$

is not given by an exponential of $\mbf{K}$. Instead it must be computed from a product integral, the general transition matrix is given by:

$$
\bsym{\Phi}(\tau,t) = \mbf{I} + \int_{\tau}^t \mbf{K}(\sigma_1)\ d\sigma_1 + \dots + \int_{\tau}^t \mbf{K}(\sigma_1) \int_{\tau}^{\sigma_1} \mbf{K}(\sigma_2) \int_{\tau}^{\sigma_2} \mbf{K}(\sigma_3) \   d\sigma_3 d\sigma_2d\sigma_1 + \cdots
$$

This means that $$\mbf{B}_k=\bsym{\Phi}(t_k,t_{k+1})$$ for $$t_0<t_2<\cdots<t_n$$. Note however that we can express the series using a [time-ordering](https://en.wikipedia.org/wiki/Time-ordering) operator $\mcal{T}$:

$$
\mbf{\Phi}(\tau,t) = \exp\mcal{T}\int_{\tau}^t \mbf{K}(\sigma)\ d\sigma
$$

We can also express the solution as the exponential of some time varying operator $\bsym{\Omega}(t,t_0)$ this way we can write the solution as 

$$
\psi(t) = \exp(\bsym{\Omega}(t,t_0))\psi(t_0) 
$$

Note however that the [derivative of the exponential map](https://en.wikipedia.org/wiki/Derivative_of_the_exponential_map) is not as simple as when considering operators as linear functions of time. The [Magnus Expansion](https://en.wikipedia.org/wiki/Magnus_expansion) provides us with a way to compute the operator $\mbf{\Omega}$ from $\mbf{K}$. The magnus expansion provides us with good approximations for the solution, in particular first and second order terms are 

$$
\begin{align}
\bsym{\Omega}_1(\tau,t) &= \int_\tau^t \mbf{K}(t_1)dt_1\\
\bsym{\Omega}_2(\tau,t) &= \frac{1}{2}\int_\tau^t dt_1\int_{\tau}^{t_1}dt_2\ [\mbf{K}(t_1),\mbf{K}(t_2)]\\
\end{align}
$$

with $[\mbf{A},\mbf{B}]\equiv \mbf{AB-BA}$ the operator commutator of $\mbf{A}$ and $\mbf{B}$. Admiting a small time step $t-\tau=\Delta t$ we consider a first order approximation of the Magnus Expansion, as such we write the discrete equation in the form 

$$
\psi_{k+1} = e^{d_k\mbf{D}}e^{\int_{t_k}^{t_{k+1}} \mbf{K}(\sigma)d\sigma }\psi_k
$$

with $\mbf{D}$ the first derivative operator. Under the assumption that the time interval $\Delta t_k = t_{k+1}-t_k$ are quite small we can assume that this discrete equation models the continuous version pretty well. Let $$\bsym{\Theta}_k=\int_{t_k}^{t_{k+1}} \mbf{K}(\sigma)d\sigma$$ then rewrite the above assume

$$
\psi_{k+1}=e^{d_k \mbf{D}} e^{\bsym{\Theta}_k}\psi_k
$$


Computationaly solving this discrete equation is quite straightforward when we discretize the operators $\mbf{D}$ and $\bsym{\Theta}_k$ by computing an appropriate matrices $\mathrm{D},\mathrm{\Theta}\in\mathbb{R}^{N\times N}$ with $N$ the number of discretization points, then tyhe operator exponential becomes matrix exponential which is straightforward to compute in standard linear algebra libraries.



### Changing the Standard Deviation

Next we want to show that the equation 

$$
\dot{\psi} = \tfrac{1}{2}\sigma^{2}(x)\mbf{G}\psi
$$

with $\mbf{G}=\frac{d^2}{dx^2}$ will shift the standard deviation by a $t\sigma^2(x)$ amount. First we consider $\sigma^2$ constant and the equation reduces to 

$$
\dot{\psi} = \tfrac{1}{2}\sigma^2\mbf{G}\psi
$$

the Kernel of the generating operator $\mbf{K} = \exp({\tfrac{t}{2}\sigma^2\mbf{G}})$ is 

$$
K(x,y,t) = \frac{1}{\sqrt{2\pi\sigma^2 t}}\exp\left(-\frac{1}{2}\frac{(x-y)^2}{\sigma^2t}\right)
$$

I want to determine the kernel of the generating operator i.e.


$$
\mbf{K}(t) = e^{\frac{1}{2}\sigma^2(x)\mbf{G}} = \exp\left(\frac{t}{2}x^2\sigma_0^2\mbf{G}\right) = ???
$$


I want to consider the cases $\sigma(x)=\sigma_0x$ . This will lead to the [Log-normal distribution](https://en.wikipedia.org/wiki/Log-normal_distribution). The goal is to find the discretization operator for the *log-normal diffusion PDE* by detemining the kernel of $\mbf{K}(t)$.


The kernel for geometric brownian motion is (Shreve Vol. II page 119. Ex. 3.6) 

$$
K(\tau,x,y) = \frac{1}{\sigma x \sqrt{2\pi\tau}}\exp\left( -\frac{1}{2}\frac{(\log(x/y)-\nu\tau)^2}{\sigma^2\tau} \right)
$$

with $\nu=\mu-\tfrac{1}{2}\sigma^2$.


The probability density function is 

$$
\psi(\tau,x) = \frac{1}{\sigma x\sqrt{2\pi \tau}}\exp\left( -\frac{1}{2}\frac{(\log(x/x_0)-\nu\tau)^2}{\sigma^2\tau} \right)
$$

when the initial distribution function is $\psi_0(x)=\delta(x-x_0)$. 


The discrete time update is determined by considering a discrete operator $\mbf{L}(k)$ with the following kernel


$$
L(k,x,y) = \frac{1}{\sigma_kx\sqrt{2\pi \Delta t}}\exp\left( -\frac{1}{2}\frac{(\log(x/y)-\nu_k\Delta t)^2}{\sigma_k^2\Delta t} \right)
$$

with $\Delta t$ the discretization period and with $\nu_k=\mu_k-\tfrac{1}{2}\sigma_k^2$. Where we assume that $\sigma_k$ and $\mu_k$ are constant inside of the time interval $t\in[k\Delta t,(k+1)\Delta t]$. The discretized equation then is of the following form

$$
\psi(k+1) = \mbf{L}(k)\psi(k)
$$

Consider a initial known value of $\psi$ at time $k$ and denote it by $\psi_k$, now determine the value of the probability in the next $T$ steps, then obtain 

$$
\phi_k = \bsym{\Pi}(k,k+T)\psi_k = \mbf{L}(k+T)\mbf{L}(k+T-1)\cdots\mbf{L}(k)\psi_k
$$

thus the transition functional operator from step $k$ to step $k+T$ is 

$$
\bsym{\Pi}_k^T\equiv\bsym{\Pi}(k,k+T) = \mbf{L}(k+T)\mbf{L}(k+T-1)\cdots\mbf{L}(k)
$$

Note how the matrix $\bsym{\Pi}_k^T$ satisfies the recursive equation

$$
\bsym{\Pi}_{k+1}^T\mbf{L}(k)=\mbf{L}(k+T+1)\bsym{\Pi}_k^T
$$

which can be used programatically to update the value of $\bsym{\Pi}_k^T$ for each new step $k$.


### Geometric Brownian Motion (GBM)

We can wirite the stochastic differential equation for the geometric brownian motion as 

$$
dS_t = \mu S_tdt+\sigma S_t dW_t
$$

with $W_t$ a [Wiener Process](https://en.wikipedia.org/wiki/Wiener_process). The solution to this equation is 

$$
S(t) = S(0)\exp(\nu t+\sigma W_t)
$$

with $\nu=\mu-(1/2)\sigma^2$. The wiener process $W_t$ satisfies 

$$
W_{t_2} - W_{t_1} = \sqrt{t_2-t_1}Z
$$

with $Z$ is an independent standard normal variable. This enables us to write the following relation

$$
\frac{S(t_2)}{S(t_1)} = \frac{\exp(\nu t_1+\sigma W_{t_1})}{\exp(\nu t_2+\sigma W_{t_2})} = \exp(\nu(t_2-t_1)+\sigma(W_{t_2} - W_{t_1})) = \exp(\nu\Delta t + \sigma\sqrt{\Delta t}Z)
$$

with $\Delta t=t_2-t_1$. Setting $t=t_1$ and $t_2=t+\Delta t$ we can establish the following result 

$$
S(t+\Delta t) = S(t) \exp(\nu\Delta t + \sigma\sqrt{\Delta t}Z)
$$

This results suggests that we can consider the discretization by writting

$$
S_{k+1} = S_k\exp(\nu_k\Delta t + \sigma_k\sqrt{\Delta t}Z_k)\label{eq:disc:update:GMB}
$$

We can solve this equation backwards up to step $k=0$, the solution is given by 

$$
S_{k} = S_0 \exp(\nu_0\Delta t + \sigma_0\sqrt{\Delta t}Z_0)\cdots \exp(\nu_{k-1}\Delta t + \sigma_{k-1}\sqrt{\Delta t}Z_{k-1})
$$

which simplifies using the properties of the exponential function

$$
S_{k} = S_0\exp\left(\sum_{j=0}^{k-1} \nu_j\Delta  t+\sigma_j\sqrt{\Delta t} Z_j  \right)
$$

since the random variables $Z_j$ are statistically independent their sum is also a normal random variable with standard deviation $$\bar{\sigma}^2_k=\sum_{j=0}^{k-1}\sigma_j^2$$. Defining $$\bar{\nu}_k=\sum_{j=0}^{k-1} \nu_j$$ we may write the above solution as the following 

$$
S_k = S_0\exp(\bar{\nu}_k\Delta t+ \bar{\sigma}_k\sqrt{\Delta t}Z)
$$

**The Expected Payoff** of a call option at time $k$ can then be computed by determining the expected value of $P_k = \max(0,S_k-K)\equiv(S_k-K)_+$ thus

$$
\mbb{E}[P_k] = \mbb{E}[(S_k-K)_+] = \frac{1}{\sqrt{2\pi}}\int_{\mbb{R}} (s_0\exp(\bar{\nu}_k\Delta t+ \bar{\sigma}_k\sqrt{\Delta t}x) -K)_+e^{-x^2/2}dx
$$

with $s_0$ a known value. Note that the integrand is non-zero only for $s_0\exp(\bar{\nu}_k\Delta t+ \bar{\sigma}_k\sqrt{\Delta t}x)-K>0$, thus we care about integrating only when the variable of integration is greater then


$$
\alpha_k\equiv\frac{\log(K/s_0) - \bar{\nu}_k \Delta t}{\bar{\sigma}_k\sqrt{\Delta t}}
$$

this means that we may rewrite the integration as

$$
\mbb{E}[P_k]  = \frac{1}{\sqrt{2\pi}}\int_{\alpha_k}^\infty (s_0\exp(\bar{\nu}_k\Delta t+ \bar{\sigma}_k\sqrt{\Delta t}x) -K)e^{-x^2/2}dx
$$

This can be evaluated by considering the cumulative distribution function $F_Z(x)$ of the standard normal random variable $Z$

$$
\mbb{E}[P_k] =   s_0e^{\bar{\nu}_k\Delta t}e^{\tfrac{1}{2}\bar{\sigma}_k\sqrt{\Delta t}}[1-F_Z(\tfrac{1}{2}\bar{\sigma}_k\sqrt{\Delta t}-\alpha_k)] -K[1-F_Z(\alpha_k)]
$$



### Estimating Parameters for GBM

Rewriting $\eqref{eq:disc:update:GMB}$ we find that 

$$
\log(S_{k+1}/S_k) = \nu_k\Delta t+\sigma_k\sqrt{\Delta t}Z_k
$$

since $Z_k$ is an i.i.d. standard normal variable we can easily deduce the parameters $\nu_k$ and $\sigma_k$. In order to provide multiple samples for each time we will consider a running model of the form 

$$
A_j = \alpha_kZ_j + \beta_k,\ \text{for }j=k,\dots,k+N
$$

with $A_k=\log(S_{k+1}/S_k)$, $\alpha_k=\sigma_k\sqrt{\Delta t}$ and $\beta_k=\nu_k\Delta t$. Let $a_k$ be samples from $A_k$ then the parameters $\alpha_k$ and $\beta_k$ can be computed from the moving average of $a_j = \log(s_{j+1}/s_j)$ and from its standard deviation, thus

$$
\hat{\beta}_k = \frac{1}{N} \sum_{j=k}^{k+N}\log(s_{j+1}/s_j), \ \quad \hat{\alpha}_k = \sqrt{\frac{1}{N} \sum_{j=k}^{k+N}(\log(s_{j+1}/s_j)-\hat{\beta}_k)^2}
$$

while the parameters $(\sigma_k,\nu_k,\mu_k)$ are 

$$
\begin{align}
\sigma_k& = \frac{\alpha_k}{\sqrt{\Delta t}}, \quad \nu_k = \frac{\beta_k}{\Delta t}\\
\mu_k& = (1/2)\sigma_k^2 + \nu_k = (1/\Delta t) ((1/2)\alpha_k^2+\beta_k)
\end{align}
$$

The chosen model works given the assumption that the values of $\alpha_k$ and $\beta_k$ do not change much in the next $N$ time steps.


### Adjusting Dividends to Present Value

Assume a continuous dividend $d(t)$ and a continuous risk-free rate $r(t)$. The value of the dividend adjusted to time $t+T$ is 

$$
\hat{d}(t) = d(t)\exp\left(\int_{t}^{t+T} r(\tau)d\tau \right)\label{eq:forward:dividend}
$$

To compute the total dividend we have to add all of the dividend values. Consider that we want to determine the total dividend value inside the time interval $t\in[t_1,t_2]$ then we compute the integral

$$
\boxed{
D(t_1,t_2) = \int_{t_1}^{t_2} \hat{d}(\theta)d\theta = \int_{t_1}^{t_2} d(\theta)\exp\left(\int_{\theta}^{\theta+T} r(\tau)d\tau \right) d\theta
}
$$


Of particular interest is the cumulative value of the dividend from $t$ to $t+T$ this sugests that we define the present value dividend $D(t)$ as 

$$
D_f(t+T) \equiv \int_{t}^{t+T} d(\theta)\exp\left(\int_{\theta}^{\theta+T} r(\tau)d\tau \right) d\theta
$$


<!-- To determine the total amount of dividends up to some value $t+T$ we simply integrate from -->


**A backward model** considers adjusting $d(t)$ to the time $t-T$, thus

$$
b(t) = d(t)\exp\left(-\int_{t-T}^{t} r(\tau)d\tau \right)
$$

note that this is the same as setting $T\rightarrow -T$ in equation $\eqref{eq:forward:dividend}$. The total amount of dividend between $t-T$ and $t$ is 

$$
D_b(t-T) \equiv \int_{t-T}^{t} d(\theta)\exp\left(-\int_{\theta-T}^{\theta} r(\tau)d\tau \right) d\theta
$$

We can define a cumulative dividend function that can determine the cumulative dividends both in the forward and backward directions, as such we define

$$
D(t+T) \equiv\left| \int_{t}^{t+T} d(\theta)\exp\left(\int_{\theta}^{\theta+T} r(\tau)d\tau \right) d\theta \right|
$$


## Functional Optimization Problems 

The goal of this section is to provide a generalization of optimization problems in [Function spaces](https://en.wikipedia.org/wiki/Function_space) such [Hilbert spaces](https://en.wikipedia.org/wiki/Hilbert_space) and [Banach spaces](https://en.wikipedia.org/wiki/Banach_space) and [$L^p$ spaces](https://en.wikipedia.org/wiki/Lp_space). A generalized functional optimization problem aims to determine the minimum of some static simple functional $F:\mcal{F}(\mcal{V})\rightarrow \mbb{R}$ under some equality or inequality constraints. 

In general we may consider the following form to express the optimization problem 

$$
\begin{split}
\underset{\psi\in\mcal{F}_\mcal{V\rightarrow U}}{\text{minimize}} &\ F(\psi)\\
\text{subject to}\ &G(\psi) = 0\\
&H(\psi) \leqslant 0
\end{split}
$$

with $H,G\in\mcal{G}_\mcal{V\rightarrow U}=\{\mcal{F}(\mcal{V})\rightarrow \mcal{F}(\mcal{U})\}$ some appropriatly chosen functionals. We may consider lagrange functionals associated to this optimization problem, indeed we can show that 

$$
\mcal{L}(\psi,\mu,\lambda) = F(\psi) + \lambda^\dagger G(\psi) + \mu^\dagger H(\psi)
$$


where $\mu,\lambda\in\mcal{F}(\mcal{U})$ are functions with domain $\mcal{U}$. We will eventually want to show that local optimum functions $\psi^*$ of the optimization problem satisfy the [Karush-Kuhn-Tucker conditions](https://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions) that is 

$$
\begin{align}
\fder \mcal{L}(\psi*) &= 0\\
G(\psi^*) &= 0, \ \  H(\psi^*) \geqslant 0\\
\mu &\geqslant 0,\ \ \mu^\dagger H(\psi^*) = 0
\end{align}
$$

The last condition is sometimes written in the equivalent form $\mu(x) H(\psi^*(x))=0,\ \ \forall x\in\mcal{V}$

A note about equality and inequalities of functions: inequalities and equalities $=,<,>,\leqslant,\geqslant$ are to be interpeted as 'element wise' for instance if we write 

$$
\psi < \phi 
$$

this is equivalent to writing 

$$
\psi(x) < \phi(x),\ \forall x\in\mcal{V}
$$

## Basis Functions and Reciprocals 

One way that we like to do computations on real hardware is by considdering a set of basis functions $\omega_1,\omega_2,\dots,\omega_N$ and their reciprocals $\omega^1,\omega^2,\dots,\omega^N$. These two sets of functions satisfy 

$$
\omega_i^\dagger \omega_j = \delta_{ij}
$$

There are certain types of functions which can be written with respect to this basis vectors in particular for $\psi$ belonging to some function space we have 

$$
\psi=\sum_i \psi_i\omega^i = \sum_i \psi^i \omega_i
$$

where 

$$
\begin{align}
\psi_i &= \omega_i^\dagger \psi \\
\psi^i &= \omega^{i\dagger}\psi
\end{align}
$$

## Discretizing Linear Operators

When we consider that functions can be expressed with respect to basis functions differential operators can be simplified greatly. Consider the linear equation 

$$
\phi = \mbf{A}\psi = \mbf{A}\omega^i \psi_i=\mbf{A}\omega_i \psi^i
$$

Now take the $j$-th component of $\phi$ and its $j$-th reciprocal component then 

$$
\begin{align}
\phi_i &= \omega_i^\dagger \phi = \omega_i^\dagger\mbf{A}\omega^j \psi_j=\omega_i^\dagger\mbf{A}\omega_j \psi^j\\
\phi^i &= \omega^{i\dagger} \phi = \omega^{i\dagger}\mbf{A}\omega^j \psi_j=\omega^{i\dagger}\mbf{A}\omega_j \psi^j
\end{align}
$$

this suggests that we have four ways to compute the matrix of $\mbf{A}$:

| $\omega_i^\dagger\mbf{A}\omega^j$|$\omega_i^\dagger\mbf{A}\omega_j$|
| $\omega^{i\dagger}\mbf{A}\omega^j$ | $\omega^{i\dagger}\mbf{A}\omega_j$ |

### Different Types of Basis Functions 

Let us consider step basis functions, let the discretization step be $\Delta x\in\mbb{R}$ then define 

$$
u(x) = 
\begin{cases}
1&\text{if}\ |x|_1 < \Delta x\\
0&\text{otherwise}
\end{cases}
$$

next define the following basis functions 

$$
\omega_i(x) = \frac{u(x-x_i)}{\sqrt{(\Delta x)^n}}
$$

And next compute 

$$
\omega_i^\dagger \omega_j = \frac{1}{(\Delta x)^n}\int_{\mbb{R}^n}  u(x-x_i)u(x-x_j) d^nx = \delta_{ij}\frac{1}{(\Delta x)^n} \int_{\mbb{R}^n} u(x)d^nx = \delta_{ij}
$$

since $\|x_i-x_j\|=\Delta x$ for $i\neq j$, and where $\int_{\mbb{R}^n} u(x)d^nx=(\Delta x)^n$ is the volume of the $n$-dimensional unit square with side $\Delta x$. Thus the functions $\omega_i$ form an orthonormal basis set and their reciprocals are $\omega^i=\omega_i$. This set of basis functions can be used to discretize functional equations!!  

**Sampling Basis Functions and Operators**

Let us start by considering sampled signal at fixed distances, namely assume a sampling domain $\mcal{S}=\\{x_1,x_2,\dots,x_N\in\mbb{R}^n\\}$ where each $x_i$ satisfies $\|x_i-x_j\|=\Delta x(1-\delta_{ij})$. Sampling can be achieved with the use of a sampling operator $\mbf{S}=\mbf{S}(\mcal{S})$, as we will see the kernel of this sampling operator can be written as 

$$
S(x,y) = \sum_{k=1}^N \delta(y-x_k)\delta(x-y)
$$

A sampled signal $\psi_d$ is expressed as the following

$$
\psi_d(x) = \sum_{y\in\mcal{S}} \delta(x-y) \psi(y)     
$$

While using the operator $\mbf{S}$ we have

$$
\psi_d = \mbf{S}\psi = \int_{\mbb{R}^n} \sum_{k=1}^N \delta(y-x_k)\delta(x-y)\psi(y)\ d^ny = \sum_{k=1}^N \delta(x-x_k)\psi_k
$$

with $\psi_k=\psi(x_k)$. The components $\psi_d^i$ of the discrete signal can also be computed by the use of the basis functions 

$$
\omega_i(x) =\delta(x-x_i) 
$$

then we have in component form that 

$$
\psi^i_d = \omega_i^\dagger \psi
$$

with each $\psi_d^i\in\mbb{R}$ and where 

$$\psi_d = \sum_i \omega_i\psi_d^i = \sum_i \omega_i\omega_i^\dagger \psi$$

Thus we can write the sampling operator in the form 

$$
\mbf{S} = \sum_i \omega_i\omega_i^\dagger
$$

Note that while $\psi$ is some arbitrary continuous function, the discrete counterpart $\psi_d$ is a sum of dirac delta functions. 



### Functional Equations with Basis Functions

A digital computer cannot solve functional equations in raw form. In order for us to provide algorithms that a computer can use to solve functional equations we have to consider some basis functions $\omega_i$. Consider the basic functional equation:

$$
\boxed{F(\psi)=0}
$$

assume that $F:\mcal{F}\mapsto\mcal{F}$. When we want to solve this equation we write $\psi$ as a linear combination of the basis functions $\omega_i$ then write 

$$
F\left(\sum_i\psi^i\omega_i\right) = 0
$$

however note that the output is also a function, which means that we must also compute components of $F(\psi)$, since $F$ can in general be a very complicated functional equation of $\psi$ there is no guarante that $F(\psi)$ can also be expressible as a linear combination of the basis functions $\omega_i$, which sugests that we consider some other basis function set $\alpha_k$, not necessarily orthogonal, then we try to solve 

$$
\alpha_j^\dagger F\left(\sum_i\psi^i\omega_i\right) = 0,\ \text{for all}\ j
$$

by finding the components $\psi^i\in\mbb{R}$ that solve this equation. 

#### Aproximating Cost Functionals with Basis Functions

Let us consider the particular case of a cost functional, that is a function that receives as input another function and outputs a real number. Let $J:\mcal{F}\mapsto \mbb{R}$. We are particularly interested in finding stationary points of $J$, in order to achieve that we will consider a taylor series approximation 

$$
J(\psi+\phi) \approx J(\psi) + \phi\cdot \fder J + \tfrac{1}{2}\phi^\dagger \fder^2 J \phi 
$$

Now consider $\phi$ expressed in terms of some basis function $\omega_i$, that is $\phi=\phi^i\omega_i$, then 

$$
J(\psi+\phi)\approx J(\psi) + \phi^i\omega_i\cdot \fder J + \tfrac{1}{2}\phi^i\phi^j \omega_i^\dagger\fder^2 J \omega_j    
$$

Notice how the linear term becomes vector-vector product of $\bsym{\phi}=(\phi^1,\dots,\phi^N)$ with $\mathsf{a}=(\omega_1\cdot \fder J,\dots,\omega_N\cdot \fder J)$, while the quadratic term becomes matrix-vector and vector-vector product, letting $[\mathsf{A}]_{ij}=\omega_i^\dagger\fder^2 J \omega_j$ we can write the cost function as 

$$
J(\psi+\phi) \approx g(\bsym{\phi})\equiv J(\psi) + \mathsf{a}^\top\bsym{\phi} + \tfrac{1}{2}\bsym{\phi}^\top \mathsf{A}\bsym{\phi}
$$

This suggests that we may minimize the approximate function with respect to the components $\phi^i$, then we compute 

$$
\nabla_{\bsym{\phi}}g(\bsym{\phi}) = 0\ \ \Leftrightarrow \  \ \bsym{\phi} = -\mathsf{A}^{-1}\mathsf{a} \label{eq:sol:min:from:components}
$$

Functional quadratic equations turn into quadratic vector equations. 

#### Zeros of functionals

Consider now the problem of finding zeros of functionals. To find zeros we can consider a first order approximation, let $F:\mcal{F}\mapsto \mcal{F}$, then 

$$
F(\psi+\phi) \approx F(\psi) + (\fder F)^\dagger\phi = 0
$$

Note how solving this for $\psi$ is equivalent to solving a simple linear functional equation, in particular we can determine $\phi$ by computing 

$$
\phi = - (\fder F)^{-\dagger} F(\psi)
$$

However we can consider a solution involving matrices instead of functions, let $\phi=\phi^i\omega_i$ and extract the components of $F(\psi+\phi)$ to find 

$$
\alpha_j^{\dagger } F(\psi+\phi^i\omega_i) \approx \alpha_j^{\dagger } F(\psi) + \alpha_j^{\dagger } (\fder F)^\dagger\omega_i \phi^i
$$

Now set $\mathsf{b}=(\alpha_j^{\dagger } F(\psi),\dots,\alpha_M^{\dagger } F(\psi))$ and $[\mathsf{B}]_{ij}=\alpha_j^{\dagger } (\fder F)^\dagger\omega_i$ to arrive at the equation 

$$
f(\bsym{\phi}) \equiv \mathsf{b} + \mathsf{B}\bsym{\phi} = 0
$$

which solving for $\bsym{\phi}$ yields 

$$
\bsym{\phi} = -\mathsf{B}^{-1}\mathsf{b}
$$

Now consider the particular case that $F$ is the functional derivative of the cost function $J$, that is 

$$
F(\psi) = \fder J(\psi)\label{eq:F:fder:J}
$$


Then, given $\eqref{eq:F:fder:J}$ we can write 

$$
F(\psi + \phi) \approx  F(\psi) + (\fder F)^\dagger\phi  = \fder J(\psi) + \fder^2 J(\psi)\phi
$$


Now compare the two ways we tried to determine solutions to the minimum of the approximate quadratic equation, on one hand we considered writing the cost functional in terms of components of $\phi$ and then minimize it with respect to those components, while on the other hand we try to solve for the functional optimality conditions $F(\psi) = \fder J(\psi)=0$, by extracting components from $F$ and requiring that $\phi$ be expressed through basis functions. The question that arises is: The solutions to the optimization problem on both cases are fundamentally different??





<!-- Then the Newtown step can be expressed in terms of a functional equation as  -->
<!---->
<!-- $$ -->
<!-- F(\psi) + \phi^\dagger(\fder F) = 0 -->
<!-- $$ -->
<!---->

## Image Sequencing Estimation

Consider that we have a sequence of RGB images encoded in some functions $\psi_0,\psi_1,\dots,\psi_N\in\mcal{F}(\mbb{R}^3\mapsto\mbb{R}^3)$ the domain of $\psi_k$ describes the position in the image, while the co-domain is the RGB color domain $[0,1]^3$. We want to estimate a linear functional $\mbf{A}\in\mcal{G}(\mbb{R}^{3\times 3})$ that transforms each element in the sequence to the next element, mathematically this can be expressed as 

$$
\psi_{k+1} \approx \mbf{A}\psi_k
$$

in order to determine the operator $\mbf{A}$ that best transforms $\psi_k$ to $\psi_{k+1}$ we have to formulate the following functional optimization problem 

$$
\underset{\mbf{A}\in\mcal{G}(\mbb{R}^{3\times 3})}{\text{minimize}} \sum_k \|\mbf{A}\psi_k-\psi_k\|^2 
$$

The norm squared can be written in the following form 

$$
\|\psi\|^2 = \int_{\mbb{R}^3} |\psi(x)|^2\ d^m {x} 
$$


$$
|\psi(x)|^2 = \psi^\top(x)\psi(x)
$$

$$
\|\mbf{A}\psi_k-\psi_k\|^2  = \int_{\mbb{R}^3} \left|[\mbf{A}\psi_k](y) - \psi_k(y)\right|^2\ d^m {y} 
$$



$$
\left|[\mbf{A}\psi_k](y) - \psi_k(y)\right|^2 = |\psi_k(y)|^2 + |[\mbf{A}\psi_k](y)|^2 - \psi_k^\top(y)[\mbf{A}\psi_k](y)
$$


$$
\psi_k^\top(y)[\mbf{A}\psi_k](y) = \psi_k^\top(y)\int_{\mbb{R}^n} A(y,z)\psi_k(z) d^mz = \int_{\mbb{R}^n}\matrixtrace(A(y,z)\psi_k(z)\psi_k(y)^\top) d^mz = \matrixtrace\left(  \int_{\mbb{R}^n}A(y,z)\psi_k(z)\psi_k(y)^\top d^mz\right)
$$

Where we used the linearity of the matrix trace function to take $\matrixtrace$ out of the integral. Now set $B(z,y)=\psi_k(z)\psi_k(y)^\top$ and then notice that the integral evaluates to 
\omega\phi^\dagger\psi = \phi\ \ \Leftrightarrow \ \  \phi^\dagger\psi = \frac{\omega^\dagger\phi}{\omega^\dagger\omega} \ \ \Leftrightarrow\ \ \psi = \frac{\phi\ (\omega^\dagger\phi)}{(\phi^\dagger\phi)(\omega^\dagger\omega)}
$$

a particular solution for $\psi$, the general solution is obtained by adding solutions of $\psi_p^\dagger\phi=0$ to the above particular solution. Given the constraint $\dot{\psi} = \mbf{A}(t)\psi$ we can provide a differential equation for $\omega$ of the form 

$$
\frac{d}{dt}\left( \frac{\phi\ (\omega^\dagger\phi)}{(\phi^\dagger\phi)(\omega^\dagger\omega)}
 \right) = \frac{ \omega^\dagger\phi}{(\phi^\dagger\phi)(\omega^\dagger\omega)}\mbf{A}\phi
$$

which certainly is not trivial to solve for $\omega$. Note however that this result is only valid when we have a finite integrable function $\phi$, note that for the case $\phi=\delta(x-s(t))$ we have $\phi^\dagger\phi = \delta(0)=+\infty$, then we have to reconsider $\eqref{eq:der=0}$, recall that $\psi(t)^\dagger\phi(t)=\psi(t,s(t))$ then $\eqref{eq:der=0}$ takes the form

$$
\frac{\delta(x-s(t))}{\psi(t,s(t))} =  \mbf{A}^\dagger(t)\lambda(t) + \dot{\lambda}(t) \label{eq:der:lambda:eq}
$$

let us rewrite the above in the form 

$$
\rho(t,x)\delta(x-s(t)) = \omega(t,x)
$$

with $\omega$ defined as above. And with $\rho(t,x)=1/\psi(t,x)$. This suggests that we can write instead 

$$
\delta(x-s(t)) = [\mbf{A}^\dagger(t)\lambda(t) + \dot{\lambda}(t)](x)\psi(t,x)
$$

Integrating on $x$ on both sides yields 

$$
1 = \psi^\dagger\mbf{A}^\dagger\lambda + \psi^\dagger\dot{\lambda} = \dot{\psi}^\dagger\lambda + \psi^\dagger\dot{\lambda}
$$

then it is straightforward to conclude that 

$$
\frac{d}{dt}(\psi^\dagger\lambda) = \text{const}
$$

thus 

$$
\psi^\dagger\lambda = t-t_0
$$


Now recall that the solution to the equation $\eqref{eq:der:lambda:eq}$ is given by 

$$
\lambda(t) = \Phi(t_0,t)\lambda_0 + \int_{t_0}^t \Phi(t_0,\tau)\theta(\tau)\ d\tau  
$$

with $\Phi$ an appropriate transition function that solely depends on $\mbf{A}(t)$, and with $\theta(t,x) = -\frac{\delta(x-s(t))}{\psi(t,s(t))}$.


## Need to motivate better!!!


consider that we measured a sample function $\phi$, the probability that this sample function comes from  the probability density $\psi$ is given by

$$
 \log P(\psi,\phi)\equiv \int_{\mbb{R}^n} \log(\psi(x)\phi(x))\ d^nx
$$

This is a generalization for when we have multiple continuous samples at each point in time??? 

