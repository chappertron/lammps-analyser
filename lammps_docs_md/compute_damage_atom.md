# compute damage/atom command

## Syntax

``` LAMMPS
compute ID group-ID damage/atom
```

-   ID, group-ID are documented in [compute](compute) command
-   damage/atom = style name of this compute command

## Examples

``` LAMMPS
compute 1 all damage/atom
```

## Description

Define a computation that calculates the per-atom damage for each atom
in a group. This is a quantity relevant for [Peridynamics
models](pair_peri). See [this document](PDF/PDLammps_overview.pdf)\_ for
an overview of LAMMPS commands for Peridynamics modeling.

The \"damage\" of a Peridynamics particles is based on the bond breakage
between the particle and its neighbors. If all the bonds are broken the
particle is considered to be fully damaged.

See the [Peridynamics Howto](Howto_peri) for a formal definition of
\"damage\" and more details about Peridynamics as it is implemented in
LAMMPS.

This command can be used with all the Peridynamic pair styles.

The damage value will be 0.0 for atoms not in the specified compute
group.

## Output info

This compute calculates a per-atom vector, which can be accessed by any
command that uses per-atom values from a compute as input. See the
[Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-atom vector values are unitless numbers (damage) $\ge 0.0$.

## Restrictions

This compute is part of the PERI package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[compute dilatation/atom](compute_dilatation_atom), [compute
plasticity/atom](compute_plasticity_atom)

## Default

none
