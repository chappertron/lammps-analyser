# compute ke/atom command

## Syntax

``` LAMMPS
compute ID group-ID ke/atom
```

-   ID, group-ID are documented in [compute](compute) command
-   ke/atom = style name of this compute command

## Examples

``` LAMMPS
compute 1 all ke/atom
```

## Description

Define a computation that calculates the per-atom translational kinetic
energy for each atom in a group.

The kinetic energy is simply $\frac12 m v^2$, where $m$ is the mass and
$v$ is the velocity of each atom.

The value of the kinetic energy will be 0.0 for atoms not in the
specified compute group.

## Output info

This compute calculates a per-atom vector, which can be accessed by any
command that uses per-atom values from a compute as input. See the
[Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-atom vector values will be in energy [units](units).

## Restrictions

> none

## Related commands

[dump custom](dump)

## Default

none
