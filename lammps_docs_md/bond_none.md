# bond_style none command

## Syntax

``` LAMMPS
bond_style none
```

## Examples

``` LAMMPS
bond_style none
```

## Description

Using a bond style of none means bond forces and energies are not
computed, even if pairs of bonded atoms were listed in the data file
read by the [read_data](read_data) command.

See the [bond_style zero](bond_zero) command for a way to calculate bond
statistics, but compute no bond interactions.

## Restrictions

> none

## Related commands

none

[bond_style zero](bond_zero)

## Default

none
