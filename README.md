# A Jekyll website for my Blog

The goal of this code is to provide notes that I am continuously taking and that I want to share with anyone. The webpage is hosted on github pages and I am always working on something. 



I need to organize all of this bullcrap!! 

- Functional Spaces:
    - Simple functional spaces (everything are scalar functions with scalar arguments)
        - Everything defined as simple  as possible
        - I care about explaining functions spaces in the most straightforward way possible 
    - Vector function spaces (functional spaces with domain $\mbb{R}^n$, vector functions and matrix functions)
        - Here I will complicate things a bit more, by considering vector functions
        - Kernels of linear operators will transform to matrices
        - Normal derivatives turn into gradients 
    - Function spaces on geometric algebras/hypercomplex numbers (It should be almost the same as vector functions spaces), 
        - Do not need to redefine linear operators notation.

*Operators Algebra:*
    - Products
        - Inner product 
        - Outer product 
        - Matrix-matrix product
    - Skew, symmetric and orthogonal operators
        - Projection operators
    - Trace of operators (and Properties)
    - Eigendecomposition of Operators 
    - Polar decomposition of operators

*Old Calculus but operator's algebra*
- Domain transformation operator $\psi(x)\mapsto\psi(f(x))$ (is a linear operator)
    - (particular case) Shift operators 
    - i-normalization operator (or functions f(x)) (Need to introduce *vector functionals*)
- Relating derivative operators with the kernel functions $[d/dx](x,y)=-\delta'(x-y)$
- Relating kernel functions with derivative operators


*Operator's Calculus:*
- Functional derivatives 
    - Directional derivatives
    - Functional divergence
    - Derivative Properties
    - Derivatives of simple functionals: $F(\psi(x))=[F(\psi)](x)$
    - Trivial derivatives of functionals
    - chain rule
        - Using chain rule to derive Euler-Poisson 
    - Taylor series expansion
    - Adjoint and differential operators 

- Extra concepts(Applications??)
    - Manifold signals and manifold estimation by parameterization
        - Estimating the optimal parameterization
        - Manifold projection is a projection operator
    - Constrained optimization of functionals
    - Linear operator equations of the form $\dot{\psi}(t) = \mbf{A}(t)\psi(t) + \rho(t)$ with for both $\mbf{A}$ time indepdent and dependent operator. 
        - Time discretization and discrete equations
        - Solutions to both discrete and continuous 
    - Functional Neural Networks
        - Equivariant and invariant Neural Networks
            - Functionals that preserve the pseudo-Eucclidean inner product
        - Rotation estimation

- Statistics an Probabilities
    - Evolution of probability densities (conditions on the transition functions)
    - Valuing Options from stochastic models

*Vector Functional Analysis*
- Gradients of functions
- Derivatives and directional derivatives of jacobian matrices 
    - Generalized Euler-Poisson 
    - Normal derivatives become jacobian matrices (need to take derivatives of jacobian matrices)


- Computational Applications
    - Discretization of linear operators 


Should I write things with theorems and proofs??? Or will it just make writing much slower??? 

Should I say gradient when talking about normal derivatives and derivative when talking about functional derivatives???

Replace the section *Functional Operators on fields* by *Vector Functional Analysis/Calculus/Algebra*. Instead of $\mcal{V}$ consider some finite pseudo-Euclidean vector space with induced metric $e_i\cdot e_j=\delta_{ij}\eta_{ii}$
