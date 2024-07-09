# bond_style mm3 command

## Syntax

``` LAMMPS
bond_style mm3
```

## Examples

``` LAMMPS
bond_style mm3
bond_coeff 1 100.0 107.0
```

## Description

The *mm3* bond style uses the potential that is anharmonic in the bond
as defined in [(Allinger)](mm3-allinger1989)

$$E = K (r - r_0)^2 \left[ 1 - 2.55(r-r_0) + \frac{7}{12} 2.55^2(r-r_0)^2 \right]$$

where $r_0$ is the equilibrium value of the bond, and $K$ is a
prefactor. The anharmonic prefactors have units $\AA^{-n}$:
$-2.55 \AA^{-1}$ and $\frac{7}{12} 2.55^2 \AA^{-2}$. The code takes care
of the necessary unit conversion for these factors internally. Note that
the MM3 papers contain an error in Eq (1): $\frac{7}{12} 2.55$ should be
replaced with $\frac{7}{12} 2.55^2$

The following coefficients must be defined for each bond type via the
[bond_coeff](bond_coeff) command as in the example above, or in the data
file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy/distance\^2)
-   $r_0$ (distance)

## Restrictions

This bond style can only be used if LAMMPS was built with the YAFF
package. See the [Build package](Build_package) doc page for more info.

## Related commands

[bond_coeff](bond_coeff)

## Default

none

------------------------------------------------------------------------

::: {#mm3-allinger1989}
**(Allinger)** Allinger, Yuh, Lii, JACS, 111(23), 8551-8566 (1989),
:::
