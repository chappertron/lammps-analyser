# compute ke/rigid command

## Syntax

``` LAMMPS
compute ID group-ID ke/rigid fix-ID
```

-   ID, group-ID are documented in [compute](compute) command
-   ke = style name of this compute command
-   fix-ID = ID of rigid body fix

## Examples

``` LAMMPS
compute 1 all ke/rigid myRigid
```

## Description

Define a computation that calculates the translational kinetic energy of
a collection of rigid bodies, as defined by one of the [fix
rigid](fix_rigid) command variants.

The kinetic energy of each rigid body is computed as
$\frac12 M V_\text{cm}^2$, where $M$ is the total mass of the rigid
body, and $V_\text{cm}$ is its center-of-mass velocity.

The *fix-ID* should be the ID of one of the [fix rigid](fix_rigid)
commands which defines the rigid bodies. The group specified in the
compute command is ignored. The kinetic energy of all the rigid bodies
defined by the fix rigid command in included in the calculation.

## Output info

This compute calculates a global scalar (the summed KE of all the rigid
bodies). This value can be used by any command that uses a global scalar
value from a compute as input. See the [Howto output](Howto_output) page
for an overview of LAMMPS output options.

The scalar value calculated by this compute is \"extensive\". The scalar
value will be in energy [units](units).

## Restrictions

This compute is part of the RIGID package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[compute erotate/rigid](compute_erotate_rigid)

## Default

none
