# fix_modify AtC computes command

## Syntax

    fix_modify <AtC fixID> computes <add|delete> <per-atom compute-ID> <volume|number>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   computes = name of the AtC sub-command
-   *add* or *delete* = add or delete the calculation of an equivalent
    continuum field for the specified per-atom compute as volume or
    number density quantity
-   per-atom compute-ID = ID of a per-atom compute; fields can be
    calculated for all per-atom computes available in LAMMPS
-   *volume* or *number* = select whether the created field is a
    per-unit-volume quantity or a per-atom quantity as weighted by
    kernel functions

## Examples

``` LAMMPS
compute virial all stress/atom
fix_modify AtC computes add virial volume
fix_modify AtC computes delete virial

compute centrosymmetry all centro/atom
fix_modify AtC computes add centrosymmetry number
```

## Description

Calculates continuum fields corresponding to specified per-atom
[computes](compute) created by LAMMPS.

## Restrictions

Must be used with [fix atc hardy](fix_atc). The per-atom compute must be
specified before the corresponding continuum field can be requested.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC fields](atc_hardy_fields)
-   [compute](compute)

## Default

None.
