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
$\newcommand{\trace}{\mathrm{Tr}}$ 
$\newcommand{\mcal}{\mathcal}$

<style>
.wrapper {
  max-width: 1500px !important;
}
</style>

{% include mathjax.html %}






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

Similarly to how functions act on vector spaces, functionals or functional operators $\mbf{F}\in\mcal{G}=\\{G:\mcal{F}(\mcal{V}^n,\mcal{U}^m)\mapsto \mcal{F}(\mcal{V}^n,\mcal{U}^m)\\}$ act on functions to produce other functions. The most simple kind of functionals are linear functionals or also called linear operators, as the name sugests they have the property that they transform functions linearly. To reduce notation and since for simplicity we want to focus on a particular case of function spaces let us consider some arbitrary vector space denoted $\Omega$ where the domain of the functions is defined, and the codomain to be $\mbb{R}$, then we shall use 

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
\phi\psi^\dagger &\in L(\Omega) 
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

where $\mbf{I}$ is the identity operator, the associated kernel function is $I(x,y)=\delta(x-y)$ where $\delta(x)\in\mcal{F}_\omega$ is the dirac-de-delta funcion defined over the vector space $\Omega$. In terms of integrals the orthogonality condition tells us that 

$$
\int_{\Omega} U(z,x)U(z,y)\ dz = \delta(x-y)
$$

Again, see the similarity with the linear algebra expression $\eqref{eq:ortho;cond:lin:alg}$!? Try to replace indices by components from $\Omega$ and the sum by an integral!




## Probability Densities

In this section I explore the analog of markov chains and probability transition matrices for probability density functions. With the use of special types of functional operators we can consider linear functional transformations that transform probability density functions to other probability density functions. A probability density function $$\rho\in\mcal{F}_\Omega$$ has the property that is positive and that integrates to one. In order for the transition probability operator $\mbf{A}$ to preserve this quantity there are certain properties that the kernel function must satisfy. Let us say that a function is normalized with respect to some function $\mbf{a}\in\mcal{F}_\omega$ when it satisfies 


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

assuming boundary terms vanish i.e. $f(t,\pm\infty)\psi(t,\pm\infty)= \frac{d}{dx}g(t,\pm\infty)^2 \psi(t,\pm\infty) = 0$.


One of the most used differential equations for the description of space-time evolution of density distributions is the fokker plank equation. We can derive a differential equation by truncating the taylor series of $\rho$ at some integer value. To see how this approximation gives rise to a differential equation we write 


$$
\begin{split}
\int_\mathbb{R} A(x,y) \rho(y) \ dy &= \int_{\mathbb{R}} A(x,x+y)\rho(x+y) \ dy = \int_{\mathbb{R}} A(x,x+y)\sum_{k=0}^\infty\frac{y^k}{k!}{\rho^{(k)}(x)} \ dy\\
&= \sum_{k=0}^\infty {\rho^{(k)}(x)}\int_{\mathbb{R}} A(x,x+y)\sum_{k=0}^\infty\frac{y^k}{k!}\ dy\approx \sum_{k=0}^N a_k(x){\rho^{(k)}(x)}
\end{split}
$$


where we made a change of variables $y\mapsto x+y$ in the first step and where we defined $a_k(x)\equiv\int_{\mathbb{R}} A(x,x+y){y^k}/{k!}\ dy$. Note that both $A(x,y)$ and $\rho(x)$could be made time $t$ independent without any change being made to the above. Then an equation of the form 

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

$$
\int_{\mathbb{R}} \delta^{(k)}(x-y)\psi(y)dy = (-1)^k\frac{\partial^{k}\psi}{\partial x^k}\label{eq:dirac:delta:derivative}
$$

where we assume that $\psi$ and all its derivatives vanish at $-\infty$ and $+\infty$. Then the application of $A$ to $\rho$ gives 

$$
A\rho = \int_{\mathbb{R}}\sum_k \delta^{(k)}(x-y) c_k(x)\rho(y)\ dy = \sum_k a_k(x) \int_{\mathbb{R}}\delta^{(k)}(x-y)\rho(y)\ dy = \sum_k (-1)^k a_k(x)\frac{\partial^{k}\rho}{\partial x^k}
$$

We can consider an aproximate solution by considering a gaussian distribution function instead of the delta-de-dirac, making it computationaly feasible in terms of integral operations. To show how the relationship between the derivatives of the delta-de-dirac function as an operator yield the signed derivative, first recall the classic integration by parts formula

$$
\int_{I} \psi \phi'\ dx = [\psi\phi]_{\partial{I}} - \int_{I}\psi'\phi\ dx
$$

with $I\in\mathbb{R}$ some interval and $\partial I$ its boundary. Applying this result to $\psi(x)\delta^{(k)}(x-y)$ yields the formula 

$$
\int_{\mathbb{R}} \delta^{(k)}(x-y)\psi(x)\ dx = [\delta^{(k)}(x-y)\psi(x)]_{-\infty}^{+\infty} - \int_{\mathbb{R}}\psi'(x)\delta^{(k-1)}(x-y)\ dx = - \int_{\mathbb{R}}\psi'(x)\delta^{(k-1)}(x-y)\ dx
$$

where $\delta^{(k)}$ plays the role of $\phi$. To obtain the final result we made $[\delta^{(k)}(x-y)\psi(x)]_{-\infty}^{+\infty}$ disapear, we argue that either $\psi$ vanishes at $\pm\infty$ or that $y$ will never attain the value of $\pm\infty$. Applying repeated integration by parts will yield the following 


$$
\int_{\mathbb{R}} \delta^{(k)}(x-y)\psi(x)\ dx = (-1)^k \int_{\mathbb{R}} \delta(x-y) \frac{\partial^k\psi}{\partial{x}^k}\ dx
$$

which after evaluating the integral on the right will give $\eqref{eq:dirac:delta:derivative}$.


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

$$
\fder = \int d\phi\  \phi'\phi\cdot \fder = \int_{\mbb{R}^\infty} d^\infty z\ \phi'(z,x)\phi(z,x)\cdot\fder  \equiv \bsym{\Phi}' \bsym{\Phi}\cdot\fder
$$

where $\phi'$ is the operational inverse of $\phi$, that is 

$$
\bsym{\Phi}'\cdot\bsym{\Phi} = 
\int_{\mbb{R}^\infty} \phi'(z,x)\phi(z,y) \ d^\infty z
= \delta(x-y)$$

Let us define the total functional derivative $\fder$ as the following functional equation:

$$
(\fder f)^\top\bsym{\phi} \equiv \int_{\mbb{R}} (\fder f)(y,x)\phi(x)\ dx = \lim_{\varepsilon\rightarrow 0} \frac{f(\psi+\varepsilon \phi)-f(\psi)}{\varepsilon} = \phi\cdot\fder f
$$

What this definition says is that applying $(\fder f)^\top$ to some function $\phi$ yields the directional functional derivative, in the $\phi$ direction, evaluated at $\psi$. 

An important property of the defined derivatives is the following 

$$
\boxed{
\fder_{\phi} \phi\cdot\fder_\psi f(\psi) = \fder_\psi f(\psi)
}
$$

note that $\phi\cdot\fder_\psi f(\psi)$ is just a linear function on $\phi$ thus the derivative transforms the linear equation $\phi\cdot\fder_\psi f(\psi)$ into the operator $\fder_\psi f(\psi)$. Another property is that if $f(\psi)$ is of rank-$k$ the  its derivative will be of rank-$(k+1)$, in other words in terms of number of arguments of functions, if $f$ is a function of $k$ arguments, i.e. $f=f(\psi,x_1,x_2,\dots,x_k)$ then $\fder f = (\fder f)(\psi,x_1,\dots,x_{k+1})$  




$$\boxed{
\textbf{Chain rule:}\ \ \fder F = \fder_\psi G(f(\psi))  = \bar{f}(\fder_\phi) G(\phi)
}\label{eq:chain:rule}$$

with $\bar{f}\equiv (\fder f)^\top$ and $\phi=f(\psi)$. And a particular case is 

$$
\fder_\psi G(\mbf{A}\psi) = \mbf{A}^\top(\fder_\phi) G(\phi) 
$$

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


A particular case of the chain rule is when I assume that $G$ is simple on $\phi$ and thus $\fder_\phi(z)G(\phi,x)=G'(\phi)\delta(x-z)$ then 

$$
\fder F = \int_{\mbb{R}} (\fder f)(y,z)\fder_\phi(z) G(\phi,x)\ dz =  \int_{\mbb{R}} (\fder f)(y,z) G'(\phi)\delta(x-z)\ dz = (\fder f)(y,x) G'(\phi,x)
$$

with $G'=\frac{\partial G}{\partial \phi}$ the classic derivative of $G$ with respect to $\phi$. Or more sucintily 

$$\boxed{\fder F=(\fder f)(y,x) G'(\phi,x)}\label{eq:part:case:simple}$$

Note how we can also consider the functional divergence which we denote $\fder\cdot F$. The divergence is simply the trace of the gradient of $F$, in other words 

$$
\boxed{
\fder \cdot F = \trace(\fder F)= \int_{\mbb{R}} (\fder F)(x,x)\ dx 
}$$

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

### Manifold Signals

We want to study signals/functions/fields that are defined only over some $m$-dimensional manifold $\mathcal{M}\subset \mbb{R}^n$. In particular we want to consider the following 

$$
\psi(x) = 
\begin{cases}
0 & \text{if}\ x\notin \mathcal{M}\\
\psi(x) & \text{otherwise}
\end{cases}
$$

This type of signal can be expressed as the application of a manifold operator $\mbb{A}=\mbb{A}(\mathcal{M})$ to some arbitrary field $\phi$. In particuar we may choose $\mbb{A}$ to be a convolutional operator such that its domain of integration is $\mathcal{M}$. Thus

$$
\mbf{A}\bsym{\phi} = \int_{\mathcal{M}} A(x,y) \phi(x) d^mx 
$$

where $A(x)$ is some arbitrary function of $x\in\mathcal{M}$. In order for $\mbf{A}$ to be a proper $\mathcal{M}$-projection operator we must ensure that $A(x,y)$ is zero when the argument is non-zero, we can easily find such a family of functions, in particular the choice $A(x,y)=a(x)\delta(x-y)$ is the proper one. And the integral then becomes

$$
\int_{\mathcal{M}} a(x)\delta(x-y) \phi(x) d^mx 
$$

Note how this simplifies when we use the properties of the $\delta$-function 

$$
\mbf{A}\bsym{\phi}  = a(y)\phi(y) \int_{\mathcal{M}} \delta(x-y) d^mx 
$$

this simplification makes now think about a new class of operator that we will call convolutional $\mathcal{M}$-manifold operators, contrary to the $\mathcal{M}$-projection operator, these operators do not project to $\mathcal{M}$ rather they are only defined over $\mathcal{M}$. So let us consider 

$$
\mbf{A}=\mbf{A}(\mathcal{M}) = \int_{\mathcal{M}}\ d^mx A(x-y)
$$

with $A:\mathbb{R}^n\mapsto \mbb{R}$. 

To distinguish between integrals and operators we will provide a different way of expressing the integrals. When we have an integral we will write the $d^mx$ at the rightmost part of the equation, but when we are just expressing an integral operator we will write $d^mx$ right after the integral $\int$ symbol. So the above expression means the integral operator and not the integral itself. 

A common problem that we encounter is determining manifolds for operators. This problem seems very abstract but it can be made very concrete by considering parameterizations of functions. In particular we shall consider that the manifold $\mathcal{M}$ can be parameterized by a function $\rho:\mbb{R}^m\mapsto\mbb{R}^n$. Then the problem of determining $\mathcal{M}$ reduces to the problem of determining $\rho$.  

By taking the change of variables $x=\rho(z)$ the differential will transform into $d^mx = \det(\underline{\rho})d^mz$ where $\underline{\rho}$ is the jacobian matrix of the transformation $x=\rho(z)$ and the determinant of the transformation has to be computed as a determinant in $m$-D??. In geometric algebra terms it is $d^mz = \|\underline{\rho}(dz^1)\wedge\cdots\wedge \underline{\rho}(dz^m)\|=\|\underline{\rho}(I_m)\|d^mz$ with $dz^k=e_kdz_k$ the directional measure, and with $I_m$ the unit pseudo-scalar of $\mbb{R}^m$. Under this change of variables the operator takes the form 

$$\boxed{
\mbf{A} = \int_{\mbb{R}^m} d^mz \det(\underline{\rho}) A(\rho(z)-y)
}$$

This is clearly non-linear equation on the function $\rho=\rho(z)$. In order to evaluate the functional derivative of this operator we will make use of the directional derivative. In particular we want to compute $\phi\cdot\fder \mbf{A}$, which using the definition is 


$$\begin{split}
\phi\cdot\fder \mbf{A}&= \lim_{\varepsilon\rightarrow 0} \int_{\mbb{R}^m} d^mz\ \frac{\det(\underline{\rho}+\varepsilon \underline{\phi}) A(\rho(z)+\varepsilon \phi(z)-y) - \det(\underline{\rho}) A(\rho(z)-y)}{\varepsilon}\\
&= \int_{\mbb{R}^m} d^mz\ \underline{\phi}\cdot\nabla_{\underline{\rho}}\det(\underline{\rho})  A(\rho(z)-y) + \det(\underline{\rho}) \phi\cdot\nabla_\rho A(\rho(z)-y) 
\end{split}
$$

where we used the product rule.



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








## LIXO 





A total derivative should be comparable to the gradient of some vector valued function. In particular if $f:\mathbb{R}^n\mapsto\mathbb{R}^m$ then the derivative is 

$$
J_{ij} = \frac{\partial f_i(\psi)}{\partial \psi_j} = e_j^\top\nabla f^\top e_i
$$

the functional derivative can then be understood by taking the $i$ and $j$ and replacing by the arguments 


The $i$-th component of the 'vector' $f$ becomes the $y$-th component $f_i=e_i^\top f\mapsto e(y)\cdot f = \delta(x-y)f$ and the $j$-th component of $\psi$ becomes the $x$-th component thus $\psi_j=e_j^\top\psi = e(x)\cdot \psi=\psi(x)$


$$
J(x,y) = e(x)\cdot\delta f\cdot e(y)
$$


with $e(x)$ an orthonormal basis set 

$$
\boxed{
e(x)\cdot e(y) = \delta(x-y)
}
$$



The unit vector $e(x)$ extracts the $x$-th component from any function of $y$. In particular 

$$
e(x)\cdot \psi(y) = \int_{\mathbb{R}} \delta(x-y)\psi(y) =   \psi(x)
$$

this result in a vector with a component at 'index' $x$.


$$
e_x\cdot \delta f\cdot e_y = \left(\lim_{\varepsilon\rightarrow 0} \frac{f(\psi+\varepsilon e_x) - f(\psi)}{\varepsilon}\right)\cdot e_y
$$



## Mais lixo

It is not that much clear what the chain rule expressed in $\eqref{eq:chain:rule}$ can be interpeted both computationaly and philosophically. First we need to understand the derivative $\fder$ operator as a vector aka a function. In some sense the product $\fder F$ is a product of $\fder(y)$ with $F(\psi,x)$, this implies that what the chain rule says is how the $\fder_\phi$ gets transformed by the jacobian transformation $\bar{f} = (\fder f)^\top$. Let us expand the integral operator of $\bar{f}(\fder_\phi)$

$$
\bar{f}(\fder_\phi) = \int_{\mbb{R}} (\fder f)(x,y)\fder_\phi(y)\ dy 
$$

Now assume that $\fder f$ is a simple differential operator of the form $\\fder f = (-1)^k \delta^{(k)}(x-y)$ where $\delta^{(k)}$ is the $k$-the derivative of the dirac function. And use the chain rule to write 

$$
\fder F = \int_{\mbb{R}} (\fder f)(y,z)\fder_\phi(z)\ dz\ G(\phi,x) = \int_{\mbb{R}} (-1)^k \delta^{(k)}(y-z) \fder_\phi(z)G(\phi,x)\ dz\label{eq:chain:rule:expanded}
$$

but 

Even though this result might not make things more clear, to try and clarify 

Then equation $\eqref{eq:chain:rule:expanded}$ becomes

$$
\fder \cdot F = \int_{\mbb{R}}\int_{\mbb{R}} (-1)^k \delta^{(k)}(x-z) \fder_\phi(z)G(\phi,x)\ dz\ dx = 
$$


### Nao e tao lixo assim mas tem de ser transformado 


The solution of this equation can be expressed via the exponential, thus

$$
\rho(t) = e^{\int_{t_0}^t \mbf{B}(\tau)d\tau}\rho_0
$$

where we



Let $\chi(x,y)$ be the transition function from point $x\in\mathbb{R}^n$ to point $y\in\mathbb{R}^n$. for this function to be a proper probability transition function all of its "colums" must sum to one this means that 

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


To better motivate the considered approach we define the following function

$$
g(\sigma,\mu,y) \equiv \frac{\exp(-\tfrac{1}{2}\sigma^{-2}(z-\mu)^2)}{\sqrt{2\pi\sigma^2}}
$$

Given a discrete time interval $\Delta t$ the time varying transition function that tells us how the distribution changes from $t$ to $t+\Delta t$ can be expressed as 

$$
A(x,y,t,\Delta t) = g(\sigma(x,t)\sqrt{\Delta t},\mu(x,t)\Delta t,y)
$$

That is we recover a distribution with standard deviation $\sigma(x,t)\sqrt{\Delta t}$ and expected value $\mu(x,t)\Delta t$. This allows us to carefully define the operator equation

$$
\rho(x,t+\Delta t) = A(t,\Delta t)\rho(x,t)
$$

where $A(t,\Delta t)\in\mathcal{A}^2$ is a time varying rank-2 operator. 

### Transition function for the Gaussian distribution

Consider that $\sigma$ and $\mu$ is constant in space and time, then applying the operator twice we get a distribution with a transformation of $\sigma$ and $\mu$

$$
A(\sigma_1^2,\mu_1)A(\sigma_2^2,\mu_2) = A(\sigma_1^2+\sigma_2^2,\mu_1+\mu_2)
$$

with the operator 

$$
A(\sigma,\mu)\rho =  \int_{\mbb{R}}\frac{\exp(-\tfrac{1}{2}\sigma^{-2}(x-z-\mu)^2)}{\sqrt{2\pi\sigma^2}}\rho(z)\ dz
$$

the transition function for the Gaussian distribution will act as a convolution operation, which is not necessarily true for transition functions of more general probability distributions. If $\rho(x)$ is expressible as a Gaussian distribution with standard deviation $\sigma'$ and expected value $\mu'$ then the resulting distribution after applying $A$ to $\rho$ will have  standard deviation ${\sigma'}^2+\sigma^2$ and expected value $\mu+\mu'$, this can be expressed as follows

$$
 A(\sigma^2,\mu)\rho({\sigma'}^2,\mu')  = \rho({\sigma'}^2+\sigma^2,\mu+\mu')
$$

And in principle for some time step $t $ we would require that 

$$
A(\sigma^2t,\mu t)\rho({\sigma'}^2,\mu')  = \rho( {\sigma'}^2+\sigma^2t,\mu'+\mu t)
$$

this expresses clearly that the time update can be done over discrete time intervals and we would always get the same result. Defining $B(t)=A(\sigma^2t,\mu t)$, then 

$$
B(t_1)B(t_2) \rho({\sigma'}^2,\mu') = B(t_1+t_2) \rho({\sigma'}^2,\mu')\rho( {\sigma'}^2+\sigma^2(t_1+t_2),\mu'+\mu(t_1+t_2) )
$$

thus as long as the standard deviation and the expected values of the transition function do not change with time we will have a consistent way to update the transition function that does not depend upon the discretization of the time-grid.

The discretized transition function reads as the following 

$$\boxed{
A(x,z,\varepsilon) = \frac{\exp(-\tfrac{1}{2}\sigma^{-2}\varepsilon^{-1}(x-z-\mu\varepsilon)^2)}{\sqrt{2\pi\sigma^2\varepsilon}}
}$$

with $\varepsilon$ some tiny time step. The more general expression will involve $x$ and $t$ varying $\sigma$ and $\mu$.


Operators that satisfy $B(t)B(\tau-t)=B(\tau)$ we call them convolutional operators or symmetric evolution operators. The gaussian operator with constant coeficients $(\sigma,\mu)$ is such a good example of a symmetric evolution operator. Note however that in real stochastic processes the operator will not have such nice property, but we can still compute the probability in a straightforward fashion. 
## Probabilistic discrete model 

Let $\mbf{A}\in L(\mbb{R}\times \mbb{R}\times \mbb{N})$ be a rank-$3$ hilbert operator. Let us denote the operator via how it acts on elements of $L(\mbb{R}\times\mbb{N})$ 


$$
\mbf{A}(k)\bsym{\psi}(k) = \int_{\mbb{R}} A(x,y,k)\psi(x,k)\ dx
$$

With this symbology we may express the update of the probability $\bsym{\psi}$ as 

$$ 
\bsym{\psi}(k+1) = \mbf{A}(k)\bsym{\psi}(k)
$$

Assuming some initial probability distribution at time $k_0$, solving the diference equation forward we have the closed form solution 

$$
\bsym{\psi}(k) = \mbf{A}(k)\mbf{A}(k-1)\cdots \mbf{A}(k_0)\bsym{\psi}_{k_0}
$$

with $\bsym{\psi}_{k_0}=\bsym{\psi}(k_0)$ the distribution at time $k_0$. The value of the probabiity at time $k$ depends linearly on the probability at time $k_0$ via the operator 

$$
\mbf{A}(k_0,k)\equiv \mbf{A}(k)\mbf{A}(k-1)\cdots \mbf{A}(k_0)
$$


The flow $\mathcal{F}$ associated with the solution of the differential equation is a function from the input at time $k_0$ to time $k$ and can be written as 

$$
\mathcal{F}(k,\bsym{\psi}(k_0)) = \bsym{\psi}(k) = \mbf{A}(k_0,k)\bsym{\psi}(k_0)
$$

The flow function tell us how a distribution gets transformed after $k-k_0$ time steps. The flow function can be interpreted in finantial terms as follows: consider that we have a stock that evolves in discrete time with probability $\bsym{\psi}(k)$. Time $k_0$ is the present, as time passes by $k_0$ increases, as $k_0$ increases we will have a new estimate value of $\bsym{\psi}(k_0)$ given all observations, in fact we will know deterministically what that value is going to be, that is the value can be expressed as $\psi(x,k_0)=\delta(x-x_{k_0})$, where $x_{k_0}$ is the observed value of the stock. So the flow tells us, given this new knowledge what is the probability in the future.

### The continualization of the discrete equations

Recall that we provided an operator $A$ that updates the distribution function for very small steps. Now if we consider non-constant coefficients the transition function will look as follows

$$\boxed{
A(x,z,t,\varepsilon) = \frac{\exp(-\tfrac{1}{2}\sigma^{-2}(z,t)\varepsilon^{-1}(x-z-\mu(z,t)\varepsilon)^2)}{\sqrt{2\pi\sigma^2(z,t)\varepsilon}}
}$$

we want to determine the equation the governs the continuous time evolution of a distribution function. As such we will rewrite the discrete equation as the following

$$
\bsym{\rho}(t+\varepsilon) = \mbf{A}(t,\varepsilon)\bsym{\rho}(t)
$$

subtracting by $\bsym{\rho}(t)$ and dividing by $\varepsilon$ gives us 

$$
\frac{\bsym{\rho}(t+\varepsilon)-\bsym{\rho}(t)}{\varepsilon} = \frac{\mbf{A}(t,\varepsilon)-\mbf{I}}{\varepsilon}\bsym{\rho}(t)
$$

as $\varepsilon$ goes to $0$ the left hand side is straightforward to evaluate, it gives the the derivative in time $t$, while the right hand side is more complicated. Let us assume that in the limit we get a new operator as such we just define 

$$
\mbf{B}(t) \equiv \lim_{\varepsilon \rightarrow 0} \frac{\mbf{A}(t,\varepsilon)-\mbf{I}}{\varepsilon} \label{eq:Boft:limitA}
$$

and the differential equation may then be written as 

$$
\dot{\bsym{\rho}}(t) = \mbf{B}(t)\bsym{\rho}(t)
$$

Note that as $\varepsilon$ goes to $0$ the transition function will become similar to a dirac delta function centered at $z+\mu(z,t)\varepsilon$, thus in some sense equation $\eqref{eq:Boft:limitA}$ will give us the derivative of the delta function. Another view can be taken by approximating the identity $\mbf{I}$ with the Gaussian transition function, and when $\varepsilon$ goes to zero we get some sort of derivative with respect to the Gaussian function. 

Even though I have not given a closed form expression for the transition function $A(t)$, I want to hightlight that as long as we are able to provide a discrete version of the transition function which nicely depends on the time step $\varepsilon$ then we are able to transform that equation to a continuous time differential equation.  

Now note that we can write a closed form solution for the differential equation if we assume $\mbf{B}$ to be constant in time $t$, then 

$$
\dot{\bsym{\rho}}(t) = \mbf{B}\bsym{\rho}(t)\longrightarrow {\bsym{\rho}}(t) = e^{\mbf{B}t}\bsym{\rho}(t_0)
$$

where the exponential function must be interepreted in a functional sense, that is, expressing the exponential as a power series and understanding that powers of $\mbf{B}$ are self compositions of the operator. For instance the square of the operator $\mbf{B}$ is just the integral


$$
\mbf{B}^2 = \int_{\mathbb{R}^n} B(y,z)B(z,x) \ dz
$$



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


