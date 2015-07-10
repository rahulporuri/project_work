This dir contains the code and results pertaining to the 
evolution of tensor perturbations during power law inflation.

We assume that a potential of the form

V = *V_0* exp(phi/*phi_0*)

drives inflation. We numerically solve the equation governing 
the evolution of the scalar field. Runge-Kutta 4 is used to solve the 
second order ODE. Using the solution for the scalar field, we 
numerically solve the equation governing the evolution of 
tensor perturbations. Again, we used RK4 to solve this second order 
ODE. Using this solution, we arrive at the 
tensor power spectrum in the super-Hubble limit.

Note that solving the aforementioned ODE, which is written in 
conformal time *eta* is numerically not plausible. We therefore 
rewrite the ODE in terms of e-fold N where

N = ln(a(t)/*a_0*)

where a(t) is the scale factor of the universe.
