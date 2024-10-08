# compute nbond/atom command

## Syntax

``` LAMMPS
compute ID group-ID nbond/atom
```

-   ID, group-ID are documented in [compute](compute) command
-   nbond/atom = style name of this compute command

## Examples

``` LAMMPS
compute 1 all nbond/atom
```

## Description

::: versionadded
4May2022
:::

Define a computation that computes the number of bonds each atom is part
of. Bonds which are broken are not counted in the tally. See the [Howto
broken bonds](Howto_bpm) page for more information. The number of bonds
will be zero for atoms not in the specified compute group. This compute
does not depend on Newton bond settings.

## Output info

This compute calculates a per-atom vector, which can be accessed by any
command that uses per-atom values from a compute as input. See the
[Howto output](Howto_output) doc page for an overview of LAMMPS output
options.

## Restrictions

This compute is part of the BPM package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

## Default

none
