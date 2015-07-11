For inflation driven by a small field potential model of the form

$$ V(\phi) = V_0\left[1-\left(\frac{\phi}{\mu}\right)^p\right]$$

we solve the equations governing the evolution of the scalar field 
driving inflation. Using these solutions, we solve the evolution of 
tensor perturbations.

The equations governing the evolution of the scalar field and the 
tensor perturbations are

$$\ddot{\phi} + 3H\dot{\phi} + \frac{{\rm d}V}{{\rm d}\phi} = 0$$

$$h_{\bf k}'' +\frac{2a'}{a}h_{\bf k}' +k^2h_{\bf k} = 0$$_

where the overprime refers to differentiation with respect to conformal time $\eta$.
For a more detailed understanding of the theory behind our work, refer to the 
.ipynb file.

Running the .py file generates data files containing the numerical 
solutions to the scalar field \phi, the Hubble parameter H, the parameter 
\eps_1 and the tensor power spectrum in the super-Hubble limit.
Running the python script file plot.py or the gnuplot script file 
small_field.plt will generate the relevant plots as .png or .eps files.
