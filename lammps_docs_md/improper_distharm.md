# improper_style distharm command

## Syntax

``` LAMMPS
improper_style distharm
```

## Examples

``` LAMMPS
improper_style distharm
improper_coeff 1 25.0 0.5
```

## Description

The *distharm* improper style uses the potential

$$E = K (d - d_0)^2$$

where $d$ is the oriented distance between the central atom and the
plane formed by the other three atoms. If the 4 atoms in an improper
quadruplet (listed in the data file read by the [read_data](read_data)
command) are ordered I,J,K,L then the L-atom is assumed to be the
central atom. Note that this is different from the convention used in
the improper_style distance. The distance $d$ is oriented and can take
on negative values. This may lead to unwanted behavior if $d_0$ is not
equal to zero.

The following coefficients must be defined for each improper type via
the improper_coeff command as in the example above, or in the data file
or restart files read by the read_data or read_restart commands:

-   $K$ (energy/distance\^2)
-   $d_0$ (distance)

------------------------------------------------------------------------

## Restrictions

This improper style can only be used if LAMMPS was built with the YAFF
package. See the [Build package](Build_package) doc page for more info.

## Related commands

[improper_coeff](improper_coeff)

## Default

none
