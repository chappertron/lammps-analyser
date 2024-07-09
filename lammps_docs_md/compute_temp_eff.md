# compute temp/eff command

## Syntax

``` LAMMPS
compute ID group-ID temp/eff
```

-   ID, group-ID are documented in [compute](compute) command
-   temp/eff = style name of this compute command

## Examples

``` LAMMPS
compute 1 all temp/eff
compute myTemp mobile temp/eff
```

## Description

Define a computation that calculates the temperature of a group of
nuclei and electrons in the [electron force field](pair_eff) model. A
compute of this style can be used by commands that compute a temperature
(e.g., [thermo_modify](thermo_modify), [fix npt/eff](fix_nh_eff)).

The temperature is calculated by the formula

$$\text{KE} = \frac{\text{dim}}{2} N k_B T,$$

where KE is the total kinetic energy of the group of atoms (sum of
$\frac12 m v^2$ for nuclei and sum of $\frac12 (m v^2 + \frac34 m s^2$)
for electrons, where $s$ includes the radial electron velocity
contributions), dim = 2 or 3 is the dimensionality of the simulation,
$N$ is the number of atoms (only total number of nuclei in the eFF (see
the [pair_eff](pair_style) command) in the group, $k_B$ is the Boltzmann
constant, and $T$ is the absolute temperature. This expression is summed
over all nuclear and electronic degrees of freedom, essentially by
setting the kinetic contribution to the heat capacity to $\frac32 k$
(where only nuclei contribute). This subtlety is valid for temperatures
well below the Fermi temperature, which for densities two to five times
the density of liquid hydrogen ranges from 86,000 to 170,000 K.

:::: note
::: title
Note
:::

For eFF models, in order to override the default temperature reported by
LAMMPS in the thermodynamic quantities reported via the [thermo](thermo)
command, the user should apply a [thermo_modify](thermo_modify) command,
as shown in the following example:
::::

``` LAMMPS
compute         effTemp all temp/eff
thermo_style    custom step etotal pe ke temp press
thermo_modify   temp effTemp
```

A six-component kinetic energy tensor is also calculated by this compute
for use in the computation of a pressure tensor. The formula for the
components of the tensor is the same as the above formula, except that
$v^2$ is replaced by $v_x v_y$ for the $xy$ component, etc. For the eFF,
again, the radial electronic velocities are also considered.

The number of atoms contributing to the temperature is assumed to be
constant for the duration of the run; use the *dynamic* option of the
[compute_modify](compute_modify) command if this is not the case.

This compute subtracts out degrees-of-freedom due to fixes that
constrain molecular motion, such as [fix shake](fix_shake) and [fix
rigid](fix_rigid). This means the temperature of groups of atoms that
include these constraints will be computed correctly. If needed, the
subtracted degrees-of-freedom can be altered using the *extra* option of
the [compute_modify](compute_modify) command.

See the [Howto thermostat](Howto_thermostat) page for a discussion of
different ways to compute temperature and perform thermostatting.

## Output info

The scalar value calculated by this compute is \"intensive\", meaning it
is independent of the number of atoms in the simulation. The vector
values are \"extensive\", meaning they scale with the number of atoms in
the simulation.

## Restrictions

This compute is part of the EFF package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[compute temp/partial](compute_temp_partial), [compute
temp/region](compute_temp_region), [compute pressure](compute_pressure)

## Default

none
