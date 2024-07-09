# compute dpd command

## Syntax

``` LAMMPS
compute ID group-ID dpd
```

-   ID, group-ID are documented in [compute](compute) command
-   dpd = style name of this compute command

## Examples

``` LAMMPS
compute 1 all dpd
```

## Description

Define a computation that accumulates the total internal conductive
energy ($U^{\text{cond}}$), the total internal mechanical energy
($U^{\text{mech}}$), the total chemical energy ($U^\text{chem}$) and the
*harmonic* average of the internal temperature ($\theta_\text{avg}$) for
the entire system of particles. See the [compute
dpd/atom](compute_dpd_atom) command if you want per-particle internal
energies and internal temperatures.

The system internal properties are computed according to the following
relations:

$$\begin{aligned}
U^\text{cond} = & \sum_{i=1}^{N} u_{i}^\text{cond} \\
U^\text{mech} = & \sum_{i=1}^{N} u_{i}^\text{mech} \\
U^\text{chem} = & \sum_{i=1}^{N} u_{i}^\text{chem} \\
            U = & \sum_{i=1}^{N} (u_{i}^\text{cond}
                  + u_{i}^\text{mech} + u_{i}^\text{chem}) \\
\theta_{avg} = & \biggl(\frac{1}{N}\sum_{i=1}^{N}
                       \frac{1}{\theta_{i}}\biggr)^{-1} \\
\end{aligned}$$

where $N$ is the number of particles in the system.

------------------------------------------------------------------------

## Output info

This compute calculates a global vector of length 5 ($U^\text{cond}$,
$U^\text{mech}$, $U^\text{chem}$, $\theta_\text{avg}$, $N$), which can
be accessed by indices 1 through 5. See the [Howto output](Howto_output)
page for an overview of LAMMPS output options.

The vector values will be in energy and temperature [units](units).

## Restrictions

This command is part of the DPD-REACT package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

This command also requires use of the [atom_style dpd](atom_style)
command.

## Related commands

[compute dpd/atom](compute_dpd_atom), [thermo_style](thermo_style)

## Default

none

------------------------------------------------------------------------

::: {#Larentzos1}
**(Larentzos)** J.P. Larentzos, J.K. Brennan, J.D. Moore, and W.D.
Mattson, \"LAMMPS Implementation of Constant Energy Dissipative Particle
Dynamics (DPD-E)\", ARL-TR-6863, U.S. Army Research Laboratory, Aberdeen
Proving Ground, MD (2014).
:::
