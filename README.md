1D heat transfer - FEM
============

We will solve the unsteady 1D heat transfer equation on a domain $x \in [0,L]$:

$$\rho c_\mathrm{p}\frac{\partial T}{\partial t} = 
\kappa \frac{\partial^2 T}{\partial x^2}$$

where $T$ is the temperature, $\rho$ is the density, $c_\mathrm{p}$ is the heat
capacity at constant pressure and $\kappa$ is the thermal conductivity. It is
assumed that initially the material is at temperature $T_0$:

$$T(T=0)=T_0$$

Three different
types of boundary conditions are considered, namely 
* Dirichlet: $T = T_\mathrm{w}$, where $T_\mathrm{w}$ is the wall temperature
* Neumann: $q=q_0$, where $q_0$ is the imposed heat flux
* Robin: $q=h(T-T_\infty)$, where $h$ is the heat transfer coefficient, and 
$T_\infty$ is the far-field temperature

NOTE: $q$ is defined as the flux **leaving** the material domain: 
$q=\kappa \partial T/\partial x$ at $x=0$, and $q=-\kappa \partial T/\partial x$ at
$x=L$).

A first order Euler scheme is used for the time derivative, thus:

$$\rho c_\mathrm{p}\frac{T_{n+1}-T_n}{\Delta t} = 
\kappa \frac{\partial^2 T_{n+1}}{\partial x^2}$$

which we rewrite to (dropping the $n+1$ indices):

$$T - \Delta t \alpha \frac{\partial^2 T}{\partial x^2} = T_n$$

where $\alpha = \kappa / (\rho c_\mathrm{p})$ is the thermal diffusivity.