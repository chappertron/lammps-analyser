# fix accelerate/cos command

## Syntax

``` LAMMPS
fix ID group-ID accelerate value
```

-   ID, group-ID are documented in [fix](fix) command
-   accelerate/cos = style name of this fix command
-   value = amplitude of acceleration (in unit of velocity/time)

## Examples

``` LAMMPS
fix 1 all accelerate/cos 0.02e-5
```

## Description

Give each atom a acceleration in x-direction based on its z coordinate.
The acceleration is a periodic function along the z-direction:

$$a_{x}(z) = A \cos \left(\frac{2 \pi z}{l_{z}}\right)$$

where $A$ is the acceleration amplitude, $l_z$ is the $z$-length of the
simulation box. At steady state, the acceleration generates a velocity
profile:

$$v_{x}(z) = V \cos \left(\frac{2 \pi z}{l_{z}}\right)$$

The generated velocity amplitude $V$ is related to the shear viscosity
$\eta$ by:

$$V = \frac{A \rho}{\eta}\left(\frac{l_{z}}{2 \pi}\right)^{2}$$

and it can be obtained from ensemble average of the velocity profile:

$$V = \frac{\sum\limits_i 2 m_{i} v_{i, x} \cos \left(\frac{2 \pi z_i}{l_{z}}\right)}{\sum\limits_i m_{i}},$$

where $m_i$, $v_{i,x}$, and $z_i$ are the mass, $x$-component velocity,
and $z$-coordinate of a particle, respectively.

The velocity amplitude $V$ can be calculated with [compute
viscosity/cos](compute_viscosity_cos), which enables viscosity
calculation with periodic perturbation method, as described by
[Hess\<Hess2\>]{.title-ref}. Because the applied acceleration drives the
system away from equilibration, the calculated shear viscosity is lower
than the intrinsic viscosity due to the shear-thinning effect.
Extrapolation to zero acceleration should generally be performed to
predict the zero-shear viscosity. As the shear stress decreases, the
signal-noise ratio decreases rapidly, the simulation time must be
extended accordingly to get converged result.

In order to get meaningful result, the group ID of this fix should be
all.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to binary restart files. None
of the fix_modify options are relevant to this fix. No global or
per-atom quantities are stored by this fix for access by various output
commands. No parameter of this fix can be used with the start/stop
keywords of the run command. This fix is not invoked during energy
minimization.

## Restrictions

This command is only available when LAMMPS was built with the MISC
package. Since this fix depends on the $z$-coordinate of atoms, it
cannot be used in 2d simulations.

## Related commands

[compute viscosity/cos](compute_viscosity_cos)

## Default

none

------------------------------------------------------------------------

::: {#Hess2}
**(Hess)** Hess, B. Journal of Chemical Physics 2002, 116 (1),
209\--217.
:::
