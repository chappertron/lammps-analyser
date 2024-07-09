# Barostats

Barostatting means controlling the pressure in an MD simulation.
[Thermostatting](Howto_thermostat) means controlling the temperature of
the particles. Since the pressure includes a kinetic component due to
particle velocities, both these operations require calculation of the
temperature. Typically a target temperature (T) and/or pressure (P) is
specified by the user, and the thermostat or barostat attempts to
equilibrate the system to the requested T and/or P.

Barostatting in LAMMPS is performed by [fixes](fix). Two barostatting
methods are currently available: Nose-Hoover (npt and nph) and
Berendsen:

-   [fix npt](fix_nh)
-   [fix npt/sphere](fix_npt_sphere)
-   [fix npt/asphere](fix_npt_asphere)
-   [fix nph](fix_nh)
-   [fix press/berendsen](fix_press_berendsen)

The [fix npt](fix_nh) commands include a Nose-Hoover thermostat and
barostat. [Fix nph](fix_nh) is just a Nose/Hoover barostat; it does no
thermostatting. Both [fix nph](fix_nh) and [fix
press/berendsen](fix_press_berendsen) can be used in conjunction with
any of the thermostatting fixes.

As with the [thermostats](Howto_thermostat), [fix npt](fix_nh) and [fix
nph](fix_nh) only use translational motion of the particles in computing
T and P and performing thermo/barostatting. [Fix
npt/sphere](fix_npt_sphere) and [fix npt/asphere](fix_npt_asphere)
thermo/barostat using not only translation velocities but also
rotational velocities for spherical and aspherical particles.

All of the barostatting fixes use the [compute
pressure](compute_pressure) compute to calculate a current pressure. By
default, this compute is created with a simple [compute
temp](compute_temp) (see the last argument of the [compute
pressure](compute_pressure) command), which is used to calculated the
kinetic component of the pressure. The barostatting fixes can also use
temperature computes that remove bias for the purpose of computing the
kinetic component which contributes to the current pressure. See the doc
pages for the individual fixes and for the [fix_modify](fix_modify)
command for instructions on how to assign a temperature or pressure
compute to a barostatting fix.

:::: note
::: title
Note
:::

As with the thermostats, the Nose/Hoover methods ([fix npt](fix_nh) and
[fix nph](fix_nh)) perform time integration. [Fix
press/berendsen](fix_press_berendsen) does NOT, so it should be used
with one of the constant NVE fixes or with one of the NVT fixes.
::::

Thermodynamic output, which can be setup via the
[thermo_style](thermo_style) command, often includes pressure values. As
explained on the page for the [thermo_style](thermo_style) command, the
default pressure is setup by the thermo command itself. It is NOT the
pressure associated with any barostatting fix you have defined or with
any compute you have defined that calculates a pressure. The doc pages
for the barostatting fixes explain the ID of the pressure compute they
create. Thus if you want to view these pressures, you need to specify
them explicitly via the [thermo_style custom](thermo_style) command. Or
you can use the [thermo_modify](thermo_modify) command to re-define what
pressure compute is used for default thermodynamic output.
