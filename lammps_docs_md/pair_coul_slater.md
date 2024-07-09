# pair_style coul/slater command

# pair_style coul/slater/cut command

# pair_style coul/slater/long command

## Syntax

``` LAMMPS
pair_style coul/slater/cut lambda cutoff
pair_style coul/slater/long lambda cutoff
```

lambda = decay length of the charge (distance units) cutoff = cutoff
(distance units)

## Examples

``` LAMMPS
pair_style coul/slater/cut 1.0 3.5
pair_coeff * *
pair_coeff 2 2 2.5

pair_style coul/slater/long 1.0 12.0
pair_coeff * *
pair_coeff 1 1 5.0
```

## Description

Styles *coul/slater* compute electrostatic interactions in mesoscopic
models which employ potentials without explicit excluded-volume
interactions. The goal is to prevent artificial ionic pair formation by
including a charge distribution in the Coulomb potential, following the
formulation of [(Melchor)](Melchor):

$$E  =  \frac{Cq_iq_j}{\epsilon r} \left( 1- \left( 1 + \frac{r_{ij}}{\lambda} exp\left( -2r_{ij}/\lambda \right) \right) \right)                       \qquad r < r_c$$

where $r_c$ is the cutoff distance and $\lambda$ is the decay length of
the charge. C is the same Coulomb conversion factor as in the
pair_styles coul/cut and coul/long. In this way the Coulomb interaction
between ions is corrected at small distances r. For the
*coul/slater/cut* style, the potential energy for distances larger than
the cutoff is zero, while for the *coul/slater/long*, the long-range
interactions are computed either by the Ewald or the PPPM technique.

Phenomena that can be captured at a mesoscopic level using this type of
electrostatic interactions include the formation of
polyelectrolyte-surfactant aggregates, charge stabilization of colloidal
suspensions, and the formation of complexes driven by charged species in
biological systems. [(Vaiwala)](Vaiwala).

The cutoff distance is optional. If it is not used, the default global
value specified in the pair_style command is used. For each pair of atom
types, a specific cutoff distance can be defined via the
[pair_coeff](pair_coeff) command as in the example above, or in the data
file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $r_c$ (distance units)

The global decay length of the charge ($\lambda$) specified in the
pair_style command is used for all pairs.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

For atom type pairs I,J and I != J, the cutoff distance for the
*coul/slater* styles can be mixed. The default mix value is *geometric*.
See the \"pair_modify\" command for details.

The [pair_modify](pair_modify) shift and table options are not relevant
for these pair styles.

These pair styles do not support the [pair_modify](pair_modify) tail
option for adding long-range tail corrections to energy and pressure.

These pair styles write their information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

## Restrictions

The *coul/slater/long* style requires the long-range solvers included in
the KSPACE package.

These styles are part of the EXTRA-PAIR package. They are only enabled
if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff), [pair_style, hybrid/overlay](pair_hybrid),
[kspace_style](kspace_style)

## Default

none

------------------------------------------------------------------------

::: {#Melchor}
**(Melchor)** Gonzalez-Melchor, Mayoral, Velazquez, and Alejandre, J
Chem Phys, 125, 224107 (2006).
:::

::: {#Vaiwala}
**(Vaiwala)** Vaiwala, Jadhav, and Thaokar, J Chem Phys, 146, 124904
(2017).
:::
