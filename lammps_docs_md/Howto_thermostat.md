# Thermostats

Thermostatting means controlling the temperature of particles in an MD
simulation. [Barostatting](Howto_barostat) means controlling the
pressure. Since the pressure includes a kinetic component due to
particle velocities, both these operations require calculation of the
temperature. Typically a target temperature (T) and/or pressure (P) is
specified by the user, and the thermostat or barostat attempts to
equilibrate the system to the requested T and/or P.

Thermostatting in LAMMPS is performed by [fixes](fix), or in one case by
a pair style. Several thermostatting fixes are available: Nose-Hoover
(nvt), Berendsen, CSVR, Langevin, and direct rescaling (temp/rescale).
Dissipative particle dynamics (DPD) thermostatting can be invoked via
the *dpd/tstat* pair style:

-   [fix nvt](fix_nh)
-   [fix nvt/sphere](fix_nvt_sphere)
-   [fix nvt/asphere](fix_nvt_asphere)
-   [fix nvt/sllod](fix_nvt_sllod)
-   [fix temp/berendsen](fix_temp_berendsen)
-   [fix temp/csvr](fix_temp_csvr)
-   [fix langevin](fix_langevin)
-   [fix temp/rescale](fix_temp_rescale)
-   [pair_style dpd/tstat](pair_dpd)

[Fix nvt](fix_nh) only thermostats the translational velocity of
particles. [Fix nvt/sllod](fix_nvt_sllod) also does this, except that it
subtracts out a velocity bias due to a deforming box and integrates the
SLLOD equations of motion. See the [Howto nemd](Howto_nemd) page for
further details. [Fix nvt/sphere](fix_nvt_sphere) and [fix
nvt/asphere](fix_nvt_asphere) thermostat not only translation velocities
but also rotational velocities for spherical and aspherical particles.

:::: note
::: title
Note
:::

A recent (2017) book by [(Daivis and Todd)](Daivis-thermostat) discusses
use of the SLLOD method and non-equilibrium MD (NEMD) thermostatting
generally, for both simple and complex fluids, e.g. molecular systems.
The latter can be tricky to do correctly.
::::

DPD thermostatting alters pairwise interactions in a manner analogous to
the per-particle thermostatting of [fix langevin](fix_langevin).

Any of the thermostatting fixes can be instructed to use custom
temperature computes that remove bias which has two effects: first, the
current calculated temperature, which is compared to the requested
target temperature, is calculated with the velocity bias removed;
second, the thermostat adjusts only the thermal temperature component of
the particle\'s velocities, which are the velocities with the bias
removed. The removed bias is then added back to the adjusted velocities.
See the doc pages for the individual fixes and for the
[fix_modify](fix_modify) command for instructions on how to assign a
temperature compute to a thermostatting fix.

For example, you can apply a thermostat only to atoms in a spatial
region by using it in conjunction with [compute
temp/region](compute_temp_region). Or you can apply a thermostat to only
the x and z components of velocity by using it with [compute
temp/partial](compute_temp_partial). Of you could thermostat only the
thermal temperature of a streaming flow of particles without affecting
the streaming velocity, by using [compute
temp/profile](compute_temp_profile).

Below is a list of custom temperature computes that can be used like
that:

-   [compute_temp_asphere]{.title-ref}
-   [compute_temp_body]{.title-ref}
-   [compute_temp_chunk]{.title-ref}
-   [compute_temp_com]{.title-ref}
-   [compute_temp_deform]{.title-ref}
-   [compute_temp_partial]{.title-ref}
-   [compute_temp_profile]{.title-ref}
-   [compute_temp_ramp]{.title-ref}
-   [compute_temp_region]{.title-ref}
-   [compute_temp_rotate]{.title-ref}
-   [compute_temp_sphere]{.title-ref}

:::: note
::: title
Note
:::

Only the nvt fixes perform time integration, meaning they update the
velocities and positions of particles due to forces and velocities
respectively. The other thermostat fixes only adjust velocities; they do
NOT perform time integration updates. Thus they should be used in
conjunction with a constant NVE integration fix such as these:
::::

-   [fix nve](fix_nve)
-   [fix nve/sphere](fix_nve_sphere)
-   [fix nve/asphere](fix_nve_asphere)

Thermodynamic output, which can be setup via the
[thermo_style](thermo_style) command, often includes temperature values.
As explained on the page for the [thermo_style](thermo_style) command,
the default temperature is setup by the thermo command itself. It is NOT
the temperature associated with any thermostatting fix you have defined
or with any compute you have defined that calculates a temperature. The
doc pages for the thermostatting fixes explain the ID of the temperature
compute they create. Thus if you want to view these temperatures, you
need to specify them explicitly via the [thermo_style
custom](thermo_style) command. Or you can use the
[thermo_modify](thermo_modify) command to re-define what temperature
compute is used for default thermodynamic output.

------------------------------------------------------------------------

::: {#Daivis-thermostat}
**(Daivis and Todd)** Daivis and Todd, Nonequilibrium Molecular Dynamics
(book), Cambridge University Press,
<https://doi.org/10.1017/9781139017848>, (2017).
:::
