# bond_style mesocnt command

## Syntax

``` LAMMPS
bond_style mesocnt
```

## Examples

``` LAMMPS
bond_style mesocnt
bond_coeff 1 C 10 10 20.0
bond_coeff 4 custom 800.0 10.0
```

## Description

::: versionadded
15Sep2022
:::

The *mesocnt* bond style is a wrapper for the [harmonic](bond_harmonic)
style, and uses the potential

$$E = K (r - r_0)^2$$

where $r_0$ is the equilibrium bond distance. Note that the usual 1/2
factor is included in $K$. The style implements parameterization presets
of $K$ for mesoscopic simulations of carbon nanotubes based on the
atomistic simulations of [(Srivastava)](Srivastava_1).

Other presets can be readily implemented in the future.

The following coefficients must be defined for each bond type via the
[bond_coeff](bond_coeff) command as in the example above, or in the data
file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   preset = *C* or *custom*
-   additional parameters depending on preset

Preset *C* is for carbon nanotubes, and the additional parameters are:

-   chiral index $n$ (unitless)
-   chiral index $m$ (unitless)
-   $r_0$ (distance)

Preset *custom* is simply a direct wrapper for the
[harmonic](bond_harmonic) style, and the additional parameters are:

-   $K$ (energy/distance\^2)
-   $r_0$ (distance)

## Restrictions

This bond style can only be used if LAMMPS was built with the MOLECULE
and MESONT packages. See the [Build package](Build_package) page for
more info.

## Related commands

[bond_coeff](bond_coeff), [delete_bonds](delete_bonds)

## Default

none

------------------------------------------------------------------------

::: {#Srivastava_1}
**(Srivastava)** Zhigilei, Wei and Srivastava, Phys. Rev. B 71, 165417
(2005).
:::
