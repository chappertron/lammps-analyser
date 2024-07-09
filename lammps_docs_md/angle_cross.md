# angle_style cross command

## Syntax

``` LAMMPS
angle_style cross
```

## Examples

``` LAMMPS
angle_style cross
angle_coeff 1 200.0 100.0 100.0 1.25 1.25 107.0
```

## Description

The *cross* angle style uses a potential that couples the bond stretches
of a bend with the angle stretch of that bend:

$$E = K_{SS} \left(r_{12}-r_{12,0}\right)\left(r_{32}-r_{32,0}\right) + K_{BS0}\left(r_{12}-r_{12,0}\right)\left(\theta-\theta_0\right) + K_{BS1}\left(r_{32}-r_{32,0}\right)\left(\theta-\theta_0\right)$$

where $r_{12,0}$ is the rest value of the bond length between atom 1 and
2, $r_{32,0}$ is the rest value of the bond length between atom 3 and 2,
and $\theta_0$ is the rest value of the angle. $K_{SS}$ is the force
constant of the bond stretch-bond stretch term and $K_{BS0}$ and
$K_{BS1}$ are the force constants of the bond stretch-angle stretch
terms.

The following coefficients must be defined for each angle type via the
[angle_coeff](angle_coeff) command as in the example above, or in the
data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K_{SS}$ (energy/distance\^2)
-   $K_{BS0}$ (energy/distance)
-   $K_{BS1}$ (energy/distance)
-   $r_{12,0}$ (distance)
-   $r_{32,0}$ (distance)
-   $\theta_0$ (degrees)

$\theta_0$ is specified in degrees, but LAMMPS converts it to radians
internally; hence the $K_{BS0}$ and $K_{BS1}$ are effectively
energy/distance per radian.

## Restrictions

This angle style can only be used if LAMMPS was built with the YAFF
package. See the [Build package](Build_package) doc page for more info.

## Related commands

[angle_coeff](angle_coeff)

## Default

none
