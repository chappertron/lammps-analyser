# angle_style gaussian command

## Syntax

``` LAMMPS
angle_style gaussian
```

## Examples

``` LAMMPS
angle_style gaussian
angle_coeff 1 300.0 2 0.0128 0.375 80.0 0.0730 0.148 123.0
```

## Description

The *gaussian* angle style uses the potential:

$$E = -k_B T ln\left(\sum_{i=1}^{n} \frac{A_i}{w_i \sqrt{\pi/2}} exp\left( \frac{-2(\theta-\theta_{i})^2}{w_i^2}\right) \right)$$

This analytical form is a suitable potential for obtaining mesoscale
effective force fields which can reproduce target atomistic
distributions [(Milano)](Milano1).

The following coefficients must be defined for each angle type via the
[angle_coeff](angle_coeff) command as in the example above, or in the
data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $T$ temperature at which the potential was derived
-   $n$ (integer \>=1)
-   $A_1$ (\> 0, radians)
-   $w_1$ (\> 0, radians)
-   $\theta_1$ (degrees)
-   \...
-   $A_n$ (\> 0, radians)
-   $w_n$ (\> 0, radians)
-   $\theta_n$ (degrees)

## Restrictions

This angle style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package. See the [Build package](Build_package) doc page
for more info.

## Related commands

[angle_coeff](angle_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Milano1}
**(Milano)** G. Milano, S. Goudeau, F. Mueller-Plathe, J. Polym. Sci. B
Polym. Phys. 43, 871 (2005).
:::
