# pair_style born/gauss command

## Syntax

``` LAMMPS
pair_style born/gauss cutoff
```

-   born/gauss = name of the pair style
-   cutoff = global cutoff (distance units)

## Examples

``` LAMMPS
pair_style born/gauss 10.0
pair_coeff 1 1 1 1 8.2464e13 12.48 0.042644277 0.44 3.56
```

## Description

::: versionadded
28Mar2023
:::

Pair style *born/gauss* computes pairwise interactions from a
combination of a Born-Mayer repulsive term and a Gaussian attractive
term according to [(Bomont)](Bomont):

$$E = A_0 \exp \left( -\alpha r \right) - A_1 \exp\left[ -\beta \left(r - r_0 \right)^2 \right]
    \qquad r < r_c$$

$r_c$ is the cutoff.

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $A_0$ (energy units)
-   $\alpha$ (1/distance units)
-   $A_1$ (energy units)
-   $\beta$ (1/(distance units)\^2)
-   $r_0$ (distance units)
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global cutoff is
used.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

This pair style does not support mixing. Thus, coefficients for all I,J
pairs must be specified explicitly.

This pair style supports the [pair_modify](pair_modify) shift option for
the energy of the pair interaction.

The [pair_modify](pair_modify) table options are not relevant for this
pair style.

This pair style does not support the [pair_modify](pair_modify) tail
option for adding long-range tail corrections to energy and pressure.

This pair style writes its information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

This pair style is only enabled if LAMMPS was built with the EXTRA-PAIR
package. See the [Build package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff), [pair_style born](pair_born)

## Default

none

------------------------------------------------------------------------

::: {#Bomont}
**(Bomont)** Bomont, Bretonnet, J. Chem. Phys. 124, 054504 (2006)
:::
