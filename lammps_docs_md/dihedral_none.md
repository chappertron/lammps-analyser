# dihedral_style none command

## Syntax

``` LAMMPS
dihedral_style none
```

## Examples

``` LAMMPS
dihedral_style none
```

## Description

Using a dihedral style of none means dihedral forces and energies are
not computed, even if quadruplets of dihedral atoms were listed in the
data file read by the [read_data](read_data) command.

See the [dihedral_style zero](dihedral_zero) command for a way to
calculate dihedral statistics, but compute no dihedral interactions.

## Restrictions

> none

## Related commands

[dihedral_style zero](dihedral_zero)

## Default

none
