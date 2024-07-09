# fix_modify AtC pair_interactions command

# fix_modify AtC bond_interactions command

## Syntax

    fix_modify <AtC fixID> pair_interactions <on|off>
    fix_modify <AtC fixID> bond_interactions <on|off>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   *pair_interactions* or *bond_interactions* = name of the AtC
    sub-command
-   *on* or *off* = activate or deactivate

## Examples

``` LAMMPS
fix_modify AtC pair_interactions off
fix_modify AtC bond_interactions on
```

## Description

Include bonds and/or pairs in stress and heat flux computations.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

*pair_interactions*: on, *bond_interactions*: off
