# fix eos/cv command

## Syntax

    fix ID group-ID eos/cv cv

-   ID, group-ID are documented in [fix](fix) command
-   eos/cv = style name of this fix command
-   cv = constant-volume heat capacity (energy/temperature units)

## Examples

``` LAMMPS
fix 1 all eos/cv 0.01
```

## Description

Fix *eos/cv* applies a mesoparticle equation of state to relate the
particle internal energy ($u_i$) to the particle internal temperature
($\theta_i$). The *eos/cv* mesoparticle equation of state requires the
constant-volume heat capacity, and is defined as follows:

$$u_{i} = u^{mech}_{i} + u^{cond}_{i} = C_{V} \theta_{i}$$

where $C_V$ is the constant-volume heat capacity, $u^{cond}$ is the
internal conductive energy, and $u^{mech}$ is the internal mechanical
energy. Note that alternative definitions of the mesoparticle equation
of state are possible.

------------------------------------------------------------------------

## Restrictions

This command is part of the DPD-REACT package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

This command also requires use of the [atom_style dpd](atom_style)
command.

## Related commands

[fix shardlow](fix_shardlow), [pair dpd/fdt](pair_dpd_fdt)

## Default

none

------------------------------------------------------------------------

::: {#Larentzos4}
**(Larentzos)** J.P. Larentzos, J.K. Brennan, J.D. Moore, and W.D.
Mattson, \"LAMMPS Implementation of Constant Energy Dissipative Particle
Dynamics (DPD-E)\", ARL-TR-6863, U.S. Army Research Laboratory, Aberdeen
Proving Ground, MD (2014).
:::
