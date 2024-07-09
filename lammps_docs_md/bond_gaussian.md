# bond_style gaussian command

## Syntax

``` LAMMPS
bond_style gaussian
```

## Examples

``` LAMMPS
bond_style gaussian
bond_coeff 1 300.0 2 0.0128 0.375 3.37 0.0730 0.148 3.63
```

## Description

The *gaussian* bond style uses the potential:

$$E = -k_B T ln\left(\sum_{i=1}^{n} \frac{A_i}{w_i \sqrt{\pi/2}} exp\left( \frac{-2(r-r_{i})^2}{w_i^2}\right)\right)$$

This analytical form is a suitable potential for obtaining mesoscale
effective force fields which can reproduce target atomistic
distributions [(Milano)](Milano0)

The following coefficients must be defined for each bond type via the
[bond_coeff](bond_coeff) command as in the example above, or in the data
file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $T$ temperature at which the potential was derived
-   $n$ (integer \>=1)
-   $A_1$ (\> 0, distance)
-   $w_1$ (\> 0, distance)
-   $r_1$ (\>= 0, distance)
-   \...
-   $A_n$ (\> 0, distance)
-   $w_n$ (\> 0, distance)
-   $r_n$ (\>= 0, distance)

## Restrictions

This bond style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package. See the [Build package](Build_package) doc page
for more info.

## Related commands

[bond_coeff](bond_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Milano0}
**(Milano)** G. Milano, S. Goudeau, F. Mueller-Plathe, J. Polym. Sci. B
Polym. Phys. 43, 871 (2005).
:::
