# compute ke/eff command

## Syntax

``` LAMMPS
compute ID group-ID ke/eff
```

-   ID, group-ID are documented in [compute](compute) command
-   ke/eff = style name of this compute command

## Examples

``` LAMMPS
compute 1 all ke/eff
```

## Description

Define a computation that calculates the kinetic energy of motion of a
group of eFF particles (nuclei and electrons), as modeled with the
[electronic force field](pair_eff).

The kinetic energy for each nucleus is computed as $\frac{1}{2} m
v^2$ and the kinetic energy for each electron is computed as
$\frac{1}{2}(m_e v^2 + \frac{3}{4} m_e s^2)$, where $m$ corresponds to
the nuclear mass, $m_e$ to the electron mass, $v$ to the translational
velocity of each particle, and $s$ to the radial velocity of the
electron, respectively.

There is a subtle difference between the quantity calculated by this
compute and the kinetic energy calculated by the *ke* or *etotal*
keyword used in thermodynamic output, as specified by the
[thermo_style](thermo_style) command. For this compute, kinetic energy
is \"translational\" and \"radial\" (only for electrons) kinetic energy,
calculated by the simple formula above. For thermodynamic output, the
*ke* keyword infers kinetic energy from the temperature of the system
with $\frac{1}{2} k_B T$ of energy for each degree of freedom. For the
eFF temperature computation via the [compute temp_eff](compute_temp_eff)
command, these are the same. But different computes that calculate
temperature can subtract out different non-thermal components of
velocity and/or include other degrees of freedom.

:::: warning
::: title
Warning
:::

The temperature in eFF models should be monitored via the [compute
temp/eff](compute_temp_eff) command, which can be printed with
thermodynamic output by using the [thermo_modify](thermo_modify)
command, as shown in the following example:
::::

``` LAMMPS
compute         effTemp all temp/eff
thermo_style    custom step etotal pe ke temp press
thermo_modify   temp effTemp
```

See [compute temp/eff](compute_temp_eff).

## Output info

This compute calculates a global scalar (the KE). This value can be used
by any command that uses a global scalar value from a compute as input.
See the [Howto output](Howto_output) page for an overview of LAMMPS
output options.

The scalar value calculated by this compute is \"extensive\". The scalar
value will be in energy [units](units).

## Restrictions

This compute is part of the EFF package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

none

## Default

none
