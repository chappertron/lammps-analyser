# angle_style mm3 command

## Syntax

``` LAMMPS
angle_style mm3
```

## Examples

``` LAMMPS
angle_style mm3
angle_coeff 1 100.0 107.0
```

## Description

The *mm3* angle style uses the potential that is anharmonic in the angle
as defined in [(Allinger)](mm3-allinger1989)

$$E = K (\theta - \theta_0)^2 \left[ 1 - 0.014(\theta - \theta_0) + 5.6(10)^{-5} (\theta - \theta_0)^2 - 7.0(10)^{-7} (\theta - \theta_0)^3 + 9(10)^{-10} (\theta - \theta_0)^4 \right]$$

where $\theta_0$ is the equilibrium value of the angle, and $K$ is a
prefactor. The anharmonic prefactors have units $\deg^{-n}$, for example
$-0.014 \deg^{-1}$, $5.6
\cdot 10^{-5} \deg^{-2}$, \...

The following coefficients must be defined for each angle type via the
[angle_coeff](angle_coeff) command as in the example above, or in the
data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy)
-   $\theta_0$ (degrees)

$\theta_0$ is specified in degrees, but LAMMPS converts it to radians
internally; hence $K$ is effectively energy per radian\^2.

## Restrictions

This angle style can only be used if LAMMPS was built with the YAFF
package. See the [Build package](Build_package) doc page for more info.

## Related commands

[angle_coeff](angle_coeff)

## Default

none
