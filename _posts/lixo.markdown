
## LIXO

# Parameterization of manifolds geometric algebra determinant


with $\partial_i=e_i\cdot\nabla$, the derivatives at all possible directions $e_i$, with $e_i$ spanning the entire vector space $\mbb{R}^m$. Even though $\partial_i\sigma(x)\in\mbb{R}^n$ 

$$
d^my = |D^my| = |\und{\sigma}(D^mz)|=|\und{\sigma}(I_m)|d^mz$$ and $$D^my = \und{\sigma}(D^mz)
$$

where we used $\mbf{d}y^i=\bar{\sigma}(\mbf{d}z^i)$ for $i=1,\dots,m$. 

Note that 

$$
\begin{align}
\bar{\sigma}(\mbf{d}z^i) &= \nabla_x \mbf{d}z^i\cdot\sigma(x) = \nabla_x \sigma_i(x)dz_i\\
\und{\sigma}(\mbf{d}z^i) &= \infd z_i\cdot \nabla_x \sigma(x) = \partial_i \sigma\  dz_i
\end{align}
$$
with $\sigma_i={e}_i\cdot\sigma$. 

$$
\infD^mz=\mbf{d}z^1\wedge \mbf{d}z^2\wedge\cdots\wedge dz^m = (e_1dz_1)\wedge(e_2dz_2)\wedge\cdots\wedge(e_ndz_m) = I_m dz_1\ dz_1\ \cdots\ dz_m = I_m d^mz
$$ 

with $d^mz \equiv dz_1\ dz_1\ \cdots\ dz_m$ the infinitesimal volume and $$dz^i=e_idz_i$$ the infinitesimal vector components in the $$i$$-th direction.



## Basis Functions



Let us now consider the discretization of a linear operator $\mbf{A}$. In order to construct the matrix of the operator $\mbf{A}$ we have to chose some basis functions $\omega_1,\omega_2,\cdots,\omega_N$ and their reciprocals $\omega^1,\omega^2,\cdots,\omega^N$  (that is $\omega_i^\ddagger\omega^j = \delta_{ij}$), then the matrix $\mathrm{A}$ of $\mbf{A}$ can be computed as 

$$
A_{ij} = \omega^{i\ddagger} \mbf{A}\omega_j 
$$

To show that we have an equivariance between both operations we chose $\psi=\psi^i\omega_i$ and compute 

$$
\phi = \mbf{A}\psi = \psi^i \mbf{A}\omega_i = \phi^i\omega_i
$$

The components $\phi_i$ of $\phi$ are computed as 

$$
\phi^j = \omega^{j\dagger} \phi = \psi^i \omega^{j\dagger}\mbf{A}\omega_i 
$$

Thus we have

$$
\phi^i  = \sum_j A_{ij}\psi^j
$$

Another way to compute the matrix of $\mbf{A}$ is 

$$
A_{ij} = \omega_i^{\ddagger} \mbf{A}\omega^j 
$$

### Operator Discretization

To consider a discretization of the operator $\mbf{A}$ let us consider a set of collocation points $x_i$ and set 

$$
\omega_i(x) = \delta(x_i-x)
$$

Note however that this choice does not form a set of orthonormal basis functions, instead they form a set of orthogonal basis functions:  

$$
\omega_i^\dagger\omega_j = \int_{\mbb{R}} \delta(x_j-x)\delta(x_i-x)\ dx = \delta(x_i-x_j) = \delta_{ij} \infty
$$

The reciprocal functions $\omega_i$ are indicator functions, they attain the value of one at a single point, consider $\mathrm{1}(0)=1$ and  $\mathrm{1}(x)=1$, $\forall x\neq 0$. 


Consider $\psi=\psi^i\omega_i$ then compute 

$$
\phi = \mbf{A}\psi = \int_{\mbb{R}} A(x,y)\psi^j\delta(y-x_j)\ dx = A(x,x_j)\psi^j
$$

The components $\phi^i$ of $\phi$ can be computed as 

$$
\phi^i = \omega_i^\dagger \phi = \int_{\mbb{R}} \delta(x-x_i) \phi(x) = \phi(x_i)
$$

$$
\phi = \phi^i\omega_i = \phi_i\omega^i
$$




## Discretizing Linear Operators

Digital computers cannot do analog computations rather they are only capable of doing digital and discretized computations. However, this limitation can be sort of overcome when the discretization is done with a discretization width so small that makes this computations behave almost as analog computations. 

Let us start by considering sampled signal at fixed distances, namely assume a sampling domain $\mcal{S}={x_1,x_2,\dots,x_N\in\mbb{R}^n}$ where each $x_i$ satisfies $\|x_i-x_j\|=\alpha(1-\delta_{ij})$. Sampling can be achieved with the use of a sampling operator $\mbf{S}=\mbf{S}(\mcal{S})$, as we will see the kernel of this sampling operator can be written as 

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

with $\psi_k=\psi(x_k)$. Now that we have provided a discretization of signals/functions we need to determine what are the linear functionals, that transform discretized signals into discretized signals. It turns out that the choice for the operator $\mbf{A}$ involves admiting the kernel 

$$
A(x,y) = \sum_{z,w\in\mcal{S}} A_{ij} \delta(x-z)\rho(y-w)
$$

with $\rho(x)=\rho(\|x\|)$ some arbitrary function (Not a dirac-delta function!). Applying $\mbf{A}$ to $\psi_d$ yields 

$$
\begin{split}
\mbf{A}\psi_d &= \int_{\mbb{R}^n}  \sum_{z,w\in\mcal{S}} A'(z,w)\delta(x-z)\rho(y-w) \sum_{u\in\mcal{S}} \delta(y-u) \psi(y)\ d^ny \\
&= \sum_{z\in\mcal{S}}\delta(x-z)\sum_{u\in\mcal{S}}\left(  \sum_{u\in\mcal{S}} A'(z,w)\rho(u-z) \right)\psi(u)\\
&= \sum_{i} \delta(x-x_i)\sum_{k}\left(  \sum_{j} A_{ij}\rho(x_j-x_k) \right)\psi_k
\end{split}
$$

The expression inside of the parenthesis is just matrix-matrix multiplication while the sum in $i$ and $k$ is matrix vector multiplication. However we are interesred in RGB signals for which $\psi$ will have $3$ components one for each color, this means that instead of each $A_{ij}$ being a scalar, they must be a $3\times 3$ matrix, and the product  $\mbf{A}\psi_d$ in discretized form will become a tensor product. 

### Manifold Signals

We want to study signals/functions/fields that are defined only over some $m$-dimensional manifold $\mathcal{M}\subset \mbb{R}^n$. In particular we want to consider the following 

$$
\psi(x) = 
\begin{cases}
0 & \text{if}\ x\notin \mathcal{M}\\
\psi(x) & \text{otherwise}
\end{cases}
$$

This type of signal can be expressed as the application of a manifold operator $\mbf{A}=\mbf{A}(\mathcal{M})$ to some arbitrary field $\phi$. In particuar we may choose $\mbb{A}$ to be a convolutional operator such that its domain of integration is $\mathcal{M}$. Thus

$$
\mbf{A}{\psi} = \int_{\mathcal{M}} A(x,y) \psi(y) d^my \label{eq:Apsi=int_M}
$$

Now consider that $\mcal{M}$ is given by a parameterization $$\mcal{M}=\{y\in\mbb{R}^n \ \\| \ y=\rho(z), \forall z\in\mbb{R}^m\}$$ with $\rho:\mbb{R}^m\mapsto \mbb{R}^n$, under the change of variables $x=\rho(z)$ the above integral becomes 


$$
\mbf{A}\psi = \int_{\mbb{R}^m} A(x,\rho(z))\psi(\rho(z)) \det(\underline{\rho})\ d^mz \label{eq:Apsi=int_M->Rm}
$$

with $\underline{\rho}$ the jacobian matrix of $\rho$, that is $a\cdot\nabla\rho=\underline{\rho} a$. When the operator $\mbf{A}$ is presented either as an integral over $\mcal{M}$ as in \eqref{eq:Apsi=int_M} or parameterized by some function $\rho$ as in $\eqref{eq:Apsi=int_M->Rm}$ it is not clear what is the kernel of the operator $\mbf{A}$. 


## Reconstruction Manifold 

As a computer vision scientist one of my goals is to find a way to reconstruct signals from raw data. Of particular interest we have the problem of reconstructing signals from samples. A sampled signal is usualy but not allways expresses as a sum of dirac-deltas. Consider that we want to find a linear operator $\mbf{A}$ that best transforms $\phi$ into $\psi$, this problem can be posed in the followin manner 


$$
\underset{\rho\in\mcal{F}_\mcal{\mbb{R}^m\rightarrow \mbb{R}^n}}{\text{minimize}} \ \|\mbf{A}(\rho)\psi-\phi\|^2
$$

where $\mbf{A}(\rho)$ is to indicate explicitily that $\mbf{A}$ is parameterized by the parameter function $\rho$. We want to find the global minimizers of the above cost function, however I think there are multiple local minimum. Nevertheless, let us try to compute its functional derivative, but first let us write the cost function as 

$$
J(\rho) = \psi^\dagger\mbf{A}^\dagger(\rho)\mbf{A}(\rho)\psi + \phi^\dagger\phi - 2\phi^\dagger \mbf{A}(\rho)\psi \label{eq:cost:functional:J(rho)}
$$

In order to compute the derivative of $J$ with respect to $\rho$ we start by determining the directional derivative in the $\theta$ direction, but first let us compute 

$$
\begin{split}
\theta\cdot\fder_\rho \mbf{A}(\rho)\psi &= \theta\cdot\fder_\rho \int_{\mbb{R}^m} A(x,\rho(z))\psi(\rho(z)) \det(\underline{\rho})\ d^mz \\
&= \int_{\mbb{R}^m} \det(\underline{\rho})
\theta(z)\cdot\nabla_\rho\ \left[A(x,\rho(z))\psi(\rho(z))\right] + A(x,\rho(z))\psi(\rho(z))\ \underline{\theta}\cdot\nabla_{\underline{\rho}}\det(\underline{\rho}) \ d^mz
\end{split}
\label{eq:dir:der:theta:J}
$$

while $\underline{\theta}$ is a matrix function of $z$, $\nabla_{\underline{\rho}}$ is the matrix derivative, and it is itself a matrix with components $$[\nabla_{\underline{\rho}}]_{ij}= \frac{\partial}{\partial \underline{\rho}_{ij}}$$, then the inner product can be written as the matrix trace $\matrixtrace(\cdot)$ of $\nabla_{\underline{\rho}}$ with $\underline{\theta}$, that is 

$$
\underline{\theta}\cdot\nabla_{\underline{\rho}} = \matrixtrace\!\left(\nabla_{\underline{\rho}}^\top\underline{\theta}\right)
$$

Let us now define the jacobian operator $\mbf{J}$, this operator acts on vector functions like $\theta$ and reuturns a matrix where each entry is the derivative with respect to each basis element, in general we shall write in equivalent form 


$$
\underline{\theta} = \mbf{J}\theta
$$

$\mbf{J}$ is tensor-valued differential operator, that is a rank-3 tensor of partial derivatives. In indexed notation we may write 

$$
[\mbf{J}\theta]_{ij} = J_{ijk}\theta^k
$$

with $J_{ijk} = \delta_{ik}\frac{\partial}{\partial{x_j}}$. With this definition we can establish $\matrixtrace\\!\left(\nabla_{\underline{\rho}}^\top\underline{\theta}\right)$ as a functional equation of $\theta$, thus

$$
\underline{\theta}\cdot\nabla_{\underline{\rho}} = \matrixtrace\!\left(\nabla_{\underline{\rho}}^\top\mbf{J}\theta\right)
$$

This suggests that we define a functional operator $\mbf{L}=\mbf{L}(\rho)$ as  

$$
\mbf{L}\theta \equiv \underline{\theta}\cdot\nabla_{\underline{\rho}} \det(\underline{\rho}) = \matrixtrace\!\left(\theta^\top\mbf{J}^{\top}\nabla_{\underline{\rho}} \det(\underline{\rho}) \right)\label{eq:Ltheta}
$$

[Jacobi's Formula](https://en.wikipedia.org/wiki/Jacobi's_formula) tells us that 

$$
\nabla_{\und{\rho}}\det(\und{\rho}) = \text{adj}(\und{\rho}) = \und{\rho}^{-1}\det(\rho).
$$

This identity provides us with a simplification of $\eqref{eq:Ltheta}$, then

$$
\mbf{L}\theta = \det(\und{\rho})\ \matrixtrace\!\left(\und{\rho}^{-\top}\mbf{J}\theta\right)
$$

With this result we may write $\eqref{eq:dir:der:theta:J}$ as

$$
\theta\cdot\fder_\rho \mbf{A}(\rho)\psi  = [\det(\und{\rho}) \nabla_{\rho}(A(x,\rho)\psi(\rho))]^\dagger\cdot\theta + [A(x,\rho)\psi(\rho)]^\dagger \mbf{L}\theta
$$

Which taking the functional derivative with respect to the function $\theta$ yields 

$$
\fder_\rho \mbf{A}(\rho)\psi  = [\det(\und{\rho}) \nabla_{\rho}(A(x,\rho)\psi(\rho))]^\dagger + [A(x,\rho)\psi(\rho)]^\dagger \mbf{L}
$$


In index notation and setting $\und{\eta}=\und{\rho}^{-\top}$ we find that 

$$
\left[\und{\rho}^{-\top}\und{\theta}\right]_{ik} = \sum_{j} \und{\eta}_{ij}\frac{\partial\theta_j}{\partial x_k}=\mbf{Q}\theta
$$

with 

$$
\mbf{Q}_{ijk} = \und{\eta}_{ij}\frac{\partial}{\partial x_k}
$$

The trace of $\und{\eta}\und{\theta}$ can be expressed in index notation as 


$$
\matrixtrace\!\left(\und{\eta}\!\ \und{\theta}\right) = \sum_{ii}\left[\und{\eta}\!\ \und{\theta}\right]_{ii} =  \sum_{ij}  \und{\eta}_{ij}\frac{\partial\theta_j}{\partial x_i} = \mbf{f}\cdot\theta 
$$

with $$\mbf{f}_j=\sum_{i} \und{\eta}_{ij}\frac{\partial}{\partial x_i} = \sum_{i} \left[\und{\rho}^{-\top}\right]_{ij}\partial_{x_i} $$ or $\mbf{f}=\und{\rho}^{-\top}\nabla$. 


$$
\und{\theta}\cdot\nabla_{\und{\rho}}f(\rho) = (\nabla_{\und{\rho}}f(\rho)\nabla)\cdot\theta
$$

where $\nabla_{\und{\rho}}f(\rho)\nabla$ is to be understood as a matrix-vector product, set $\mathsf{F}=\nabla_{\und{\rho}}f(\rho)$ then we have $\mathsf{F}(\nabla)\cdot\theta$. Where it is to be understood that $\mathsf{F}$ transforms the nabla operator $\nabla$ linearly.


The directional derivative of the considered functional $\beta(\rho)=\mbf{A}(\rho)\psi $ can be written as the following 

$$
\theta\cdot\fder_\rho \beta(\rho) = \theta^\dagger\cdot \mbf{B} + \theta^\dagger \cdot \mbf{G}(\nabla)\label{eq:dir:der:easy}
$$

with the kernel $B$ of $\mbf{B}$ and the kernel $G$ of $\mbf{G}$ given by

$$
\begin{align}
B(x,z) &= [\det(\und{\rho})](x)\nabla_\rho [A(x,\rho(z))\psi(\rho(z))]\\
G(x,z) &= A(x,\rho(z))\psi(\rho(z)) [\nabla_{\und{\rho}}\det(\und{\rho})](z)
\end{align}
$$

Let us rewrite the cost functional using the identity $\beta(\rho)=\mbf{A}(\rho)\psi $ and $\eqref{eq:dir:der:easy}$, then

$$
J(\rho) = \beta^\dagger(\rho)\beta(\rho) + \phi^\dagger\phi - 2\phi^\dagger\beta(\rho)
$$

Now note that 

$$
\fder_\rho\beta^\dagger(\rho) = \fder_\theta \theta\cdot\fder_\rho\beta(\rho) = \mbf{B} + \mbf{G}\nabla
$$

then using the product rule on the cost functional we find that 

$$
\fder_\rho^\dagger J(\rho)  = 2\beta(\rho)^\dagger [\fder_\rho\beta^\dagger(\rho)]^\dagger - 2\phi^\dagger [\fder_\rho\beta^\dagger(\rho)]^\dagger
$$

which reduces to 

$$
\boxed{
(1/2)\fder_\rho J(\rho) = (\mbf{B} + \mbf{G}(\nabla))(\beta(\rho) - \phi)
}
$$





Let us admit the following simplifications:

- In $\eqref{eq:cost:functional:J(rho)}$ a simplification arises when we consider that $\mbf{A}(\rho)$ is orthogonal under $\mcal{M}$ and thus $\psi^\dagger\mbf{A}^\dagger(\rho)\mbf{A}(\rho)\psi=\psi^\dagger\psi$ an example of this is when $A(x,y)=\delta(x-y)$
- Another simplification also arises when we coinsider that the determinant of $\underline{\rho}$ is constant
- Consider $A(x,y)=\delta(x-y)$

Now note that for any function $f$ we have $\nabla_x f(x-y) = -\nabla_y f(x-y)$ this means that 

$$
\theta\cdot\nabla_\rho \delta(x-\rho(z)) = -\theta\cdot\nabla_x\delta(x-\rho(z))
$$



**Multivariate Euler-Lagrange Equations** 

We want to show that the multivariate Euler-Lagrange equations can be written in the form 

$$
\nabla_\psi \mcal{L}-\nabla_x^\top\nabla_{\und{\psi}} = 0
$$


Let 

$$
J(\psi) = \int_{\mbb{R}^n} \mcal{L}(\psi,\und{\psi}) d^nx
$$

Next compute the directional derivative of $J$ in the $\phi$ direction 

$$
\phi\cdot\fder J(\psi) = \lim_{\varepsilon\rightarrow 0}  \int_{\mbb{R}^n} \frac{\mcal{L}(\psi+\varepsilon\phi,\und{\psi}+\varepsilon\und{\phi})-\mcal{L}(\psi,\und{\psi})}{\varepsilon} d^nx
$$

which simplifies into 

$$
\phi\cdot\fder J(\psi)  =   \int_{\mbb{R}^n} (\phi\cdot\nabla_\psi + \und{\phi}\cdot\nabla_{\und{\psi}}) \mcal{L}(\psi,\und{\psi})\ d^nx
$$

Now use $\eqref{eq:matrix:differential:trace}$ and note that $\nabla_{\und{\psi}}$ is a matrix of derivatives, then write:

$$
\und{\phi}\cdot\nabla_{\und{\psi}} = (\nabla_\und{\psi}^\top\nabla_x)\cdot\phi
$$

Next recall that by the chain-rule derivatives of function become minus the derivatives this means that 

$$
\fder_\phi \phi\cdot\fder J(\psi)  = \nabla_\psi \mcal{L}-\nabla_x^\top\nabla_{\und{\psi}}  \mcal{L} 
$$


Consider the following function $F(\psi) = G(\mbf{A}\psi)$ with $\mbf{A}$ some arbitrary linear operator. Now we want to compute the functional derivative of $F$ with respect to $\psi$, assume that $G$ is a simple functional and that $\mbf{A}\psi$ is a scalar function. The chain-rule tells us that 

$$
\fder_\psi G(\mbf{A}\psi) = \mbf{A}^\ddagger \fder_{\alpha} G(\alpha)|_{\alpha=\mbf{A}\psi}  
$$

Now compute the integral 

$$
\begin{split}
\mbf{A}^\ddagger \fder_{\alpha} G(\alpha) &= \int_{\mbb{R}^n}[\mbf{A}^\top](z,x) \fder_\alpha(z) [G(\alpha)](y)dz = \int_{\mbb{R}^n}[\mbf{A}^\top](z,x)[G'(\alpha)](y)\delta(z-y)dz \\ 
&= [G'(\alpha)](y)\int_{\mbb{R}^n}[\mbf{A}^\top](z,x)\delta(z-y)dz
\end{split}
$$

<!-- now let the kernel $A$ of $\mbf{A}(x,y)$ be $A=\mathsf{A}(x,y)^\top\nabla\delta(x-y)$ then  -->


now let the kernel $A$ of $\mbf{A}(x,y)$ be $A(x,y)=-\nabla_x\delta(x-y)$ then 

$$
\begin{split}
\mbf{A}^\ddagger \fder_{\alpha} G(\alpha)  &= -[G'(\alpha)](y)\int_{\mbb{R}^n} \nabla_z^\top  \delta(z-x) \delta(z-y) dz\\
&= [G'(\alpha)](y) \nabla_x^\top  \int_{\mbb{R}^n}  \delta(z-x) \delta(z-y)dz = [G'(\alpha)](y) \nabla_x^\top \delta(y-x)
\end{split}
$$







## Some other things


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

Note that now the integral does no longer depend on the manifold $\mcal{M}$ that we wish to determine, this suggests that the kernel function of $\mbf{A}$ will be 

$$
A(z,y) = \det(\underline{\rho}) A(\rho(z)-y)
$$

In order to evaluate the functional derivative of $\mbf{A}$ we need to compute the directional derivative of the kernel of $\mbf{A}$ then

$$
\phi\cdot\fder A(z,y) = \underline{\phi}\cdot\nabla_{\underline{\rho}}\det(\underline{\rho})  A(\rho(z)-y) + \det(\underline{\rho}) \phi\cdot\nabla_\rho A(\rho(z)-y) 
$$

Where $\underline{\phi}\cdot\nabla_{\underline{\rho}}\equiv \trace(\underline{\phi}^\dagger\nabla_{\underline{\rho}})$ and it is to be understood that $\underline{\phi}$ and $\nabla_{\underline{\rho}}$ are of *degree*-2. 


### Finding the sub-Optimal manifold

Now we want to discuss one such application of the manifold signal problem. Consider that we have two functions $\phi$ and $\psi$ 






We showed how we can take some manifold express a manifold signal as a parameterization and then express the problem of determining the best manifold by considering the parameterization instead. 


Then we may express the kernel function of 


This is clearly non-linear equation on the function $\rho=\rho(z)$. In order to evaluate the functional derivative of this operator we will make use of the directional derivative. In particular we want to compute $\phi\cdot\fder \mbf{A}$, which using the definition is 


$$\begin{split}
\phi\cdot\fder \mbf{A}&= \lim_{\varepsilon\rightarrow 0} \int_{\mbb{R}^m} d^mz\ \frac{\det(\underline{\rho}+\varepsilon \underline{\phi}) A(\rho(z)+\varepsilon \phi(z)-y) - \det(\underline{\rho}) A(\rho(z)-y)}{\varepsilon}\\
&= \int_{\mbb{R}^m} d^mz\ \underline{\phi}\cdot\nabla_{\underline{\rho}}\det(\underline{\rho})  A(\rho(z)-y) + \det(\underline{\rho}) \phi\cdot\nabla_\rho A(\rho(z)-y) 
\end{split}
$$

where we used the product rule.

###  Another attempt
Let us consider instead some scalar valued function $\alpha(x)$ of $x\in\mbb{R}^n$ that expresses the eigenvalues of the operator $\mbf{B}$ then write the spectral decomposition as 

$$
\mbf{B} = \int_{\mbb{R}^n}\alpha(x) \psi_{x}\psi_{x}^\dagger d^n x
$$

now compute the exponential by considering that $\psi_x^\dagger\psi_{x'}=\delta(x-x')$ thus 

$$
\exp(\mbf{B}) = \int_{\mbb{R}^n} \exp(\alpha(x)) \psi_{x}\psi_{x}^\dagger d^n x
$$

Now assume that $\mbf{A}$ has the decomposition 

$$
\mbf{A} = \int_{\mbb{R}^n} \lambda(x)\psi_{x}\psi_{x}^\dagger d^nx
$$

then setting $\mbf{A} = \exp(\mbf{B})$ yields the equation 

$$
\int_{\mbb{R}^n} \lambda(x)\psi_{x}\psi_{x}^\dagger d^nx = \int_{\mbb{R}^n} \exp(\alpha(x)) \psi_{x}\psi_{x}^\dagger d^n x\ \ \Leftrightarrow\ \  \int_{\mbb{R}^n} [\exp(\alpha(x))-\lambda(x)] \psi_{x}\psi_{x}^\dagger d^n x = 0
$$

which has as unique solution $\lambda(x)=\exp(\alpha(x))$. 




<!-- now consider the change of variables $\beta=\exp(\alpha)$ then  -->
<!---->
<!-- $$ -->
<!-- d\beta = d \exp(\alpha) = \exp(\alpha)d\alpha = \beta d\alpha -->
<!-- $$ -->
<!---->
<!-- By which we reccover the weird result: -->
<!---->
<!-- $$ -->
<!-- \exp(\mbf{B}) = \int_{R}  \psi_{\beta}\psi_{\beta}^\dagger d\beta -->
<!-- $$ -->
<!---->
<!-- with $S=\exp(R)$ and $\psi_\beta=\phi_{\log(\beta)}$, a reparameterization of the orthogonal functions. The exponential of a symmetric operator is a projector!!?? -->
<!---->
<!---->



