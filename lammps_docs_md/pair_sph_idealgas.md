# pair_style sph/idealgas command

## Syntax

``` LAMMPS
pair_style sph/idealgas
```

## Examples

``` LAMMPS
pair_style sph/idealgas
pair_coeff * * 1.0 2.4
```

## Description

The sph/idealgas style computes pressure forces between particles
according to the ideal gas equation of state:

$$p = (\gamma - 1) \rho e$$

where $\gamma = 1.4$ is the heat capacity ratio, $\rho$ is the local
density, and e is the internal energy per unit mass. This pair style
also computes Monaghan\'s artificial viscosity to prevent particles from
interpenetrating [(Monaghan)](ideal-Monoghan).

See [this PDF guide](PDF/SPH_LAMMPS_userguide.pdf)\_ to using SPH in
LAMMPS.

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above.

-   $\nu$ artificial viscosity (no units)
-   h kernel function cutoff (distance units)

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

This style does not support mixing. Thus, coefficients for all I,J pairs
must be specified explicitly.

This style does not support the [pair_modify](pair_modify) shift, table,
and tail options.

This style does not write information to [binary restart
files](restart). Thus, you need to re-specify the pair_style and
pair_coeff commands in an input script that reads a restart file.

This style can only be used via the *pair* keyword of the [run_style
respa](run_style) command. It does not support the *inner*, *middle*,
*outer* keywords.

## Restrictions

This pair style is part of the SPH package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[pair_coeff](pair_coeff), pair_sph/rhosum

## Default

none

------------------------------------------------------------------------

::: {#ideal-Monoghan}
**(Monaghan)** Monaghan and Gingold, Journal of Computational Physics,
52, 374-389 (1983).
:::
