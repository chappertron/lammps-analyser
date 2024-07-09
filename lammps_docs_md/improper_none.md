# improper_style none command

## Syntax

``` LAMMPS
improper_style none
```

## Examples

``` LAMMPS
improper_style none
```

## Description

Using an improper style of none means improper forces and energies are
not computed, even if quadruplets of improper atoms were listed in the
data file read by the [read_data](read_data) command.

See the [improper_style zero](improper_zero) command for a way to
calculate improper statistics, but compute no improper interactions.

## Restrictions

> none

## Related commands

[improper_style zero](improper_zero)

## Default

none
