# Calculate temperature

Temperature is computed as kinetic energy divided by some number of
degrees of freedom (and the Boltzmann constant). Since kinetic energy is
a function of particle velocity, there is often a need to distinguish
between a particle\'s advection velocity (due to some aggregate motion
of particles) and its thermal velocity. The sum of the two is the
particle\'s total velocity, but the latter is often what is wanted to
compute a temperature.

LAMMPS has several options for computing temperatures, any of which can
be used in [thermostatting](Howto_thermostat) and
[barostatting](Howto_barostat). These [compute commands](compute)
calculate temperature:

-   [compute temp](compute_temp)
-   [compute temp/sphere](compute_temp_sphere)
-   [compute temp/asphere](compute_temp_asphere)
-   [compute temp/com](compute_temp_com)
-   [compute temp/deform](compute_temp_deform)
-   [compute temp/partial](compute_temp_partial)
-   [compute temp/profile](compute_temp_profile)
-   [compute temp/ramp](compute_temp_ramp)
-   [compute temp/region](compute_temp_region)

All but the first 3 calculate velocity biases directly (e.g. advection
velocities) that are removed when computing the thermal temperature.
[Compute temp/sphere](compute_temp_sphere) and [compute
temp/asphere](compute_temp_asphere) compute kinetic energy for
finite-size particles that includes rotational degrees of freedom. They
both allow for velocity biases indirectly, via an optional extra
argument which is another temperature compute that subtracts a velocity
bias. This allows the translational velocity of spherical or aspherical
particles to be adjusted in prescribed ways.
