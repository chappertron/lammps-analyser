# compute temp/body command

## Syntax

``` LAMMPS
compute ID group-ID temp/body keyword value ...
```

-   ID, group-ID are documented in [compute](compute) command

-   temp/body = style name of this compute command

-   zero or more keyword/value pairs may be appended

-   keyword = *bias* or *dof*

        *bias* value = bias-ID
          bias-ID = ID of a temperature compute that removes a velocity bias
        *dof* value = *all* or *rotate*
          all = compute temperature of translational and rotational degrees of freedom
          rotate = compute temperature of just rotational degrees of freedom

## Examples

``` LAMMPS
compute 1 all temp/body
compute myTemp mobile temp/body bias tempCOM
compute myTemp mobile temp/body dof rotate
```

## Description

Define a computation that calculates the temperature of a group of body
particles, including a contribution from both their translational and
rotational kinetic energy. This differs from the usual [compute
temp](compute_temp) command, which assumes point particles with only
translational kinetic energy.

Only body particles can be included in the group. For 3d particles, each
has 6 degrees of freedom (3 translational, 3 rotational). For 2d body
particles, each has 3 degrees of freedom (2 translational, 1
rotational).

:::: note
::: title
Note
:::

This choice for degrees of freedom (DOF) assumes that all body particles
in your model will freely rotate, sampling all their rotational DOF. It
is possible to use a combination of interaction potentials and fixes
that induce no torque or otherwise constrain some of all of your
particles so that this is not the case. Then there are less DOF and you
should use the [compute_modify extra](compute_modify) command to adjust
the DOF accordingly.
::::

The translational kinetic energy is computed the same as is described by
the [compute temp](compute_temp) command. The rotational kinetic energy
is computed as $\frac12 I \omega^2$, where $I$ is the moment of inertia
tensor for the aspherical particle and $\omega$ is its angular velocity,
which is computed from its angular momentum.

A kinetic energy tensor, stored as a 6-element vector, is also
calculated by this compute. The formula for the components of the tensor
is the same as the above formula, except that $v^2$ and $\omega^2$ are
replaced by $v_x v_y$ and $\omega_x \omega_y$ for the math:\`xy\`
component, and the appropriate elements of the inertia tensor are used.
The six components of the vector are ordered $xx$, $yy$, $zz$, $xy$,
$xz$, $yz$.

The number of atoms contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
[compute_modify](compute_modify) command if this is not the case.

This compute subtracts out translational degrees-of-freedom due to fixes
that constrain molecular motion, such as [fix shake](fix_shake) and [fix
rigid](fix_rigid). This means the temperature of groups of atoms that
include these constraints will be computed correctly. If needed, the
subtracted degrees-of-freedom can be altered using the *extra* option of
the [compute_modify](compute_modify) command.

See the [Howto thermostat](Howto_thermostat) page for a discussion of
different ways to compute temperature and perform thermostatting.

------------------------------------------------------------------------

The keyword/value option pairs are used in the following ways.

For the *bias* keyword, *bias-ID* refers to the ID of a temperature
compute that removes a \"bias\" velocity from each atom. This allows
compute temp/sphere to compute its thermal temperature after the
translational kinetic energy components have been altered in a
prescribed way (e.g., to remove a flow velocity profile). Thermostats
that use this compute will work with this bias term. See the doc pages
for individual computes that calculate a temperature and the doc pages
for fixes that perform thermostatting for more details.

For the *dof* keyword, a setting of *all* calculates a temperature that
includes both translational and rotational degrees of freedom. A setting
of *rotate* calculates a temperature that includes only rotational
degrees of freedom.

------------------------------------------------------------------------

## Output info

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1\--6.
These values can be used by any command that uses global scalar or
vector values from a compute as input. See the [Howto
output](Howto_output) page for an overview of LAMMPS output options.

The scalar value calculated by this compute is \"intensive\". The vector
values are \"extensive\".

The scalar value will be in temperature [units](units). The vector
values will be in energy [units](units).

## Restrictions

This compute is part of the BODY package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

This compute requires that atoms store angular momentum and a quaternion
as defined by the [atom_style body](atom_style) command.

## Related commands

[compute temp](compute_temp)

## Default

none
