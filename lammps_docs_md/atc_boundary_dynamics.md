# fix_modify AtC boundary_dynamics command

## Syntax

    fix_modify <AtC fixID> boundary_dynamics <on|damped_harmonic|prescribed|coupled|none>

-   AtC fixID = ID of [fix atc](fix_atc) instance
-   boundary_dynamics = name of the AtC sub-command
-   *on* or *damped_harmonic* *prescribed* *coupled* *none*

## Description

Sets different schemes for controlling boundary atoms. *on* will
integrate the boundary atoms using the velocity-verlet algorithm.
*damped_harmonic* uses a mass/spring/dashpot for the boundary atoms with
added arguments of the damping and spring constants followed by the
ratio of the boundary type mass to the desired mass. *prescribed* forces
the boundary atoms to follow the finite element displacement. *coupled*
does the same.

## Restrictions

Boundary atoms must be specified. When using swaps between internal and
boundary atoms, the initial configuration must have already correctly
partitioned the two.

## Related AtC commands

-   [fix_modify AtC command overview](atc_fix_modify)

## Default

*prescribed on*
