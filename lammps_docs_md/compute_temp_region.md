# compute temp/region command

## Syntax

``` LAMMPS
compute ID group-ID temp/region region-ID
```

-   ID, group-ID are documented in [compute](compute) command
-   temp/region = style name of this compute command
-   region-ID = ID of region to use for choosing atoms

## Examples

``` LAMMPS
compute mine flow temp/region boundary
```

## Description

Define a computation that calculates the temperature of a group of atoms
in a geometric region. This can be useful for thermostatting one portion
of the simulation box. For example, a McDLT simulation where one side is
cooled, and the other side is heated. A compute of this style can be
used by any command that computes a temperature (e.g.,
[thermo_modify](thermo_modify), [fix temp/rescale](fix_temp_rescale)).

Note that a *region*-style temperature can be used to thermostat with
[fix temp/rescale](fix_temp_rescale) or [fix langevin](fix_langevin),
but should probably not be used with Nose\--Hoover style fixes ([fix
nvt](fix_nh), [fix npt](fix_nh), or [fix nph](fix_nh)) if the degrees of
freedom included in the computed temperature vary with time.

The temperature is calculated by the formula

$$\text{KE} = \frac{\text{dim}}{2} N k_B T,$$

where KE = is the total kinetic energy of the group of atoms (sum of
$\frac12 m v^2$), dim = 2 or 3 is the dimensionality of the simulation,
$N$ is the number of atoms in both the group and region, $k_B$ is the
Boltzmann constant, and $T$ temperature.

A kinetic energy tensor, stored as a six-element vector, is also
calculated by this compute for use in the computation of a pressure
tensor. The formula for the components of the tensor is the same as the
above formula, except that $v^2$ is replaced by $v_x v_y$ for the $xy$
component, and so on. The six components of the vector are ordered $xx$,
$yy$, $zz$, $xy$, $xz$, $yz$.

The number of atoms contributing to the temperature is calculated each
time the temperature is evaluated since it is assumed atoms can
enter/leave the region. Thus there is no need to use the *dynamic*
option of the [compute_modify](compute_modify) command for this compute
style.

The removal of atoms outside the region by this fix is essentially
computing the temperature after a \"bias\" has been removed, which in
this case is the velocity of any atoms outside the region. If this
compute is used with a fix command that performs thermostatting then
this bias will be subtracted from each atom, thermostatting of the
remaining thermal velocity will be performed, and the bias will be added
back in. Thermostatting fixes that work in this way include [fix
nvt](fix_nh), [fix temp/rescale](fix_temp_rescale), [fix
temp/berendsen](fix_temp_berendsen), and [fix langevin](fix_langevin).
This means that when this compute is used to calculate the temperature
for any of the thermostatting fixes via the [fix modify
temp](fix_modify) command, the thermostat will operate only on atoms
that are currently in the geometric region.

Unlike other compute styles that calculate temperature, this compute
does not subtract out degrees-of-freedom due to fixes that constrain
motion, such as [fix shake](fix_shake) and [fix rigid](fix_rigid). This
is because those degrees of freedom (e.g., a constrained bond) could
apply to sets of atoms that straddle the region boundary, and hence the
concept is somewhat ill-defined. If needed the number of subtracted
degrees of freedom can be set explicitly using the *extra* option of the
[compute_modify](compute_modify) command.

See the [Howto thermostat](Howto_thermostat) page for a discussion of
different ways to compute temperature and perform thermostatting.

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

> none

## Related commands

[compute temp](compute_temp), [compute pressure](compute_pressure)

## Default

none
