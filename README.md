1D heat transfer - FEM
============

We will solve the unsteady 1D heat transfer equation on a domain $x \in [0,L]$:

$$\rho c_\mathrm{p}\frac{\partial T}{\partial t} = 
\kappa \frac{\partial^2 T}{\partial x^2}$$

where $T$ is the temperature, $\rho$ is the density, $c_\mathrm{p}$ is the
specific heat capacity at constant pressure and $\kappa$ is the thermal 
conductivity. It is assumed that initially the material is at temperature $T_0$:

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

which we rewrite to (dropping the $n+1$ indices for readibility):

$$T - \Delta t \alpha \frac{\partial^2 T}{\partial x^2} = T_n$$

where $\alpha = \kappa / (\rho c_\mathrm{p})$ is the thermal diffusivity.

Next, the integral form of the differential equation is written as: find $T$ such that
$$\int v  T dl - \int \Delta t \alpha v \frac{\partial^2 T}{\partial x^2} dl = 
\int v T_n dl$$
for all $v$.

From integration by parts, we get the following weak form: find $T$ such that
$$\int v  T dl + \int \Delta t \alpha \frac{\partial v}{\partial x} 
\frac{\partial T}{\partial x} dl - \Delta t \alpha \left.v\frac{\partial T}{\partial x}\right|_{x=L} + 
\Delta t \alpha \left.v\frac{\partial T}{\partial x}\right|_{x=0} = \int v T_n dl$$
for all $v$.

Neumann boundary conditions
--------
For the case of Neumann boundary conditions, we can fill in the flux directly in
the boundary terms. Rewriting the defition of the flux as given above: 

* at $x=0$ we have $\alpha \partial T/\partial x=q_0/(\rho c_\mathrm{p})$
* at $x=L$ we have $-\alpha \partial T/\partial x=q_0/(\rho c_\mathrm{p})$

This yields:

$$\int v  T dl + \int \Delta t \alpha \frac{\partial v}{\partial x} 
\frac{\partial T}{\partial x} dl = - \Delta t \bar{q}_0 \left.v\right|_{x=L} 
- \Delta t \bar{q}_0 \left.v\right|_{x=0} + \int v T_n dl$$

where $\bar{q}_0 = q_0/(\rho c_\mathrm{p})$.


Robin boundary conditions
--------
For the case of Robin boundary conditions, we can fill in the flux directly in
the boundary terms. Rewriting the defition of the flux as given above:

* at $x=0$ we have $\alpha \partial T/\partial x=h(T-T_\infty)/(\rho c_\mathrm{p})$
* at $x=L$ we have $-\alpha \partial T/\partial x=h(T-T_\infty)/(\rho c_\mathrm{p})$


$$\int v  T dl + \int \Delta t \alpha \frac{\partial v}{\partial x} 
\frac{\partial T}{\partial x} dl 
+ \left.\Delta t \bar{h} T \right|_{x=L}
+ \left.\Delta t \bar{h} T \right|_{x=0}
= \Delta t \left.v \bar{h} T_\infty \right|_{x=L} 
+ \Delta t \left.v \bar{h} T_\infty\right|_{x=0} + \int v T_n dl$$

where $\bar{h} = h/(\rho c_\mathrm{p})$.
