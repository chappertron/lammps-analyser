# dielectric command

## Syntax

``` LAMMPS
dielectric value
```

-   value = dielectric constant

## Examples

``` LAMMPS
dielectric 2.0
```

## Description

Set the dielectric constant for Coulombic interactions (pairwise and
long-range) to this value. The constant is unitless, since it is used to
reduce the strength of the interactions. The value is used in the
denominator of the formulas for Coulombic interactions (e.g., a value of
4.0 reduces the Coulombic interactions to 25% of their default
strength). See the [pair_style](pair_style) command for more details.

## Restrictions

> none

## Related commands

[pair_style](pair_style)

## Default

``` LAMMPS
dielectric 1.0
```
