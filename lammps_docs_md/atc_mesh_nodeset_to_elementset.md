# fix_modify AtC mesh nodeset_to_elementset command

## Syntax

    fix_modify <AtC fixID> mesh nodeset_to_elementset <nodeset_id> <elementset_id> <max/min>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   mesh nodeset_to_elementset = name of the AtC sub-command
-   nodeset_id = id of desired nodeset from which to create the
    elementset
-   elementset_id = id to assign to the collection of FE elements
-   \<max/min\> = flag to choose either the maximal or minimal
    elementset

## Examples

``` LAMMPS
fix_modify AtC mesh nodeset_to_elementset myNodeset myElementset min
```

## Description

Command to create an elementset from an existing nodeset. Either the
minimal element set of elements with all nodes in the set, or maximal
element set with all elements with at least one node in the set, can be
created.

## Restrictions

None.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)
-   [fix_modify AtC mesh create_elementset](atc_mesh_create_elementset)
-   [fix_modify AtC mesh delete_elements](atc_mesh_delete_elements)

## Default

Unless specified, the maximal element set is created.
