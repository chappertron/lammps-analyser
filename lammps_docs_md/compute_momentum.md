# compute momentum command

## Syntax

``` LAMMPS
compute ID group-ID momentum
```

-   ID, group-ID are documented in [compute](compute) command
-   momentum = style name of this compute command

## Examples

``` LAMMPS
compute 1 all momentum
```

## Description

Define a computation that calculates the translational momentum *p* of a
group of particles. It is computed as the sum
$\vec{p} = \sum_i m_i \cdot \vec{v}_i$ over all particles in the compute
group, where *m* and *v* are the mass and velocity vector of the
particle, respectively.

## Output info

This compute calculates a global vector (the summed momentum) of length
3. This value can be used by any command that uses a global vector value
from a compute as input. See the [Howto output](Howto_output) page for
an overview of LAMMPS output options.

The vector value calculated by this compute is \"extensive\". The vector
value will be in mass\*velocity [units](units).

## Restrictions

This compute is part of the EXTRA-COMPUTE package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

## Default

none
