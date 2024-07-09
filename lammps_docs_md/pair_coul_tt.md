# pair_style coul/tt command

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *coul/tt*
-   args = list of arguments for a particular style

<!-- -->

    *coul/tt* args = n cutoff
      n = degree of polynomial
      cutoff = global cutoff (distance units)

## Examples

``` LAMMPS
pair_style hybrid/overlay ... coul/tt 4 12.0
pair_coeff 1 2  coul/tt 4.5 1.0
pair_coeff 1 2  coul/tt 4.0 1.0 4 12.0
pair_coeff 1 3* coul/tt 4.5 1.0 4
```

Example input scripts available: examples/PACKAGES/drude

## Description

The *coul/tt* pair style is meant to be used with force fields that
include explicit polarization through Drude dipoles.

The *coul/tt* pair style should be used as a sub-style within in the
[pair_style hybrid/overlay](pair_hybrid) command, in conjunction with a
main pair style including Coulomb interactions and *thole* pair style,
or with *lj/cut/thole/long* pair style that is equivalent to the
combination of preceding two.

The *coul/tt* pair styles compute the charge-dipole Coulomb interaction
damped at short distances by a function

$$f_{n,ij}(r) = 1 - c_{ij} \cdot e^{-b_{ij} r} \sum_{k=0}^n \frac{(b_{ij} r)^k}{k!}$$

This function results from an adaptation to the Coulomb interaction
[(Salanne)](Salanne1) of the damping function originally proposed by
[Tang Toennies](TangToennies1) for van der Waals interactions.

The polynomial takes the degree 4 for damping the Coulomb interaction.
The parameters $b_{ij}$ and $c_{ij}$ could be determined from
first-principle calculations for small, mainly mono-atomic, ions
[(Salanne)](Salanne1), or else treated as empirical for large molecules.

In pair styles with Drude induced dipoles, this damping function is
typically applied to the interactions between a Drude charge (either
$q_{D,i}$ on a Drude particle or $-q_{D,i}$ on the respective Drude
core)) and a charge on a non-polarizable atom, $q_{j}$.

The Tang-Toennies function could also be used to damp electrostatic
interactions between the (non-polarizable part of the) charge of a core,
$q_{i}-q_{D,i}$, and the Drude charge of another, $-q_{D,j}$. The
$b_{ij}$ and $c_{ij}$ are equal to $b_{ji}$ and $c_{ji}$ in the case of
core-core interactions.

For pair_style *coul/tt*, the following coefficients must be defined for
each pair of atoms types via the [pair_coeff](pair_coeff) command as in
the example above.

-   $b_{ij}$
-   $c_{ij}$
-   degree of polynomial (positive integer)
-   cutoff (distance units)

The last two coefficients are optional. If not specified the global
degree of the polynomial or the global cutoff specified in the
pair_style command are used. In order to specify a cutoff (forth
argument), the degree of the polynomial (third argument) must also be
specified.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

The *coul/tt* pair style does not support mixing. Thus, coefficients for
all I,J pairs must be specified explicitly.

## Restrictions

These pair styles are part of the DRUDE package. They are only enabled
if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

This pair_style should currently not be used with the [charmm dihedral
style](dihedral_charmm) if the latter has non-zero 1-4 weighting
factors. This is because the *coul/tt* pair style does not know which
pairs are 1-4 partners of which dihedrals.

## Related commands

[fix drude](fix_drude), [fix langevin/drude](fix_langevin_drude), [fix
drude/transform](fix_drude_transform), [compute
temp/drude](compute_temp_drude), [pair_style thole](pair_thole)

## Default

none

------------------------------------------------------------------------

::: {#Thole3}
**(Thole)** Chem Phys, 59, 341 (1981).
:::

::: {#Salanne1}
**(Salanne)** Salanne, Rotenberg, Jahn, Vuilleumier, Simon, Christian
and Madden, Theor Chem Acc, 131, 1143 (2012).
:::

::: {#TangToennies1}
**(Tang and Toennies)** J Chem Phys, 80, 3726 (1984).
:::
