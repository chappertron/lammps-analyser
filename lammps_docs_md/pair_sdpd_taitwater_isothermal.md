# pair_style sdpd/taitwater/isothermal command

## Syntax

``` LAMMPS
pair_style sdpd/taitwater/isothermal temperature viscosity seed
```

-   temperature = temperature of the fluid (temperature units)
-   viscosity = dynamic viscosity of the fluid (mass\*distance/time
    units)
-   seed = random number generator seed (positive integer, optional)

## Examples

``` LAMMPS
pair_style sdpd/taitwater/isothermal 300. 1. 28681
pair_coeff * * 1000.0 1430.0 2.4
```

## Description

The sdpd/taitwater/isothermal style computes forces between mesoscopic
particles according to the Smoothed Dissipative Particle Dynamics model
described in this paper by [(Espanol and Revenga)](Espanol_Revenga)
under the following assumptions:

1.  The temperature is constant and uniform.
2.  The shear viscosity is constant and uniform.
3.  The volume viscosity is negligible before the shear viscosity.
4.  The Boltzmann constant is negligible before the heat capacity of a
    single mesoscopic particle of fluid.

The third assumption is true for water in nearly incompressible flows.
The fourth holds true for water for any reasonable size one can imagine
for a mesoscopic particle.

The pressure forces between particles will be computed according to
Tait\'s equation of state:

$$p = B \left[(\frac{\rho}{\rho_0})^{\gamma} - 1\right]$$

where $\gamma = 7$ and $B = c_0^2 \rho_0 / \gamma$, with $\rho_0$ being
the reference density and $c_0$ the reference speed of sound.

The laminar viscosity and the random forces will be computed according
to formulas described in [(Espanol and Revenga)](Espanol_Revenga).

:::: warning
::: title
Warning
:::

Similar to [brownian](pair_brownian) and [dpd](pair_dpd) styles, the
[newton](newton) setting for pairwise interactions needs to be on when
running LAMMPS in parallel if you want to ensure linear momentum
conservation. Otherwise random forces generated for pairs straddling
processor boundary will not be equal and opposite.
::::

:::: note
::: title
Note
:::

The actual random seed used will be a mix of what you specify and other
parameters like the MPI ranks. This is to ensure that different MPI
tasks have distinct seeds.
::::

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above.

-   $\rho_0$ reference density (mass/volume units)
-   $c_0$ reference soundspeed (distance/time units)
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

This pair style is part of the DPD-SMOOTH package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair coeff](pair_coeff), [pair sph/rhosum](pair_sph_rhosum), [pair
sph/taitwater](pair_sph_taitwater)

## Default

The default seed is 0 (before mixing).

------------------------------------------------------------------------

::: {#Espanol_Revenga}
**(Espanol and Revenga)** Espanol, Revenga, Physical Review E, 67,
026705 (2003).
:::
