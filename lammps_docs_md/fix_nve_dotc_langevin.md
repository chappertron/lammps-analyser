# fix nve/dotc/langevin command

## Syntax

    fix ID group-ID nve/dotc/langevin Tstart Tstop damp seed keyword value

-   ID, group-ID are documented in [fix](fix) command

-   nve/dotc/langevin = style name of this fix command

-   Tstart,Tstop = desired temperature at start/end of run (temperature
    units)

-   damp = damping parameter (time units)

-   seed = random number seed to use for white noise (positive integer)

-   keyword = *angmom*

        *angmom* value = factor
          factor = do thermostat rotational degrees of freedom via the angular momentum and apply numeric scale factor as discussed below

## Examples

``` LAMMPS
fix 1 all nve/dotc/langevin 1.0 1.0 0.03 457145 angmom 10
fix 1 all nve/dotc/langevin 0.1 0.1 78.9375 457145 angmom 10
```

## Description

Apply a rigid-body Langevin-type integrator of the kind \"Langevin C\"
as described in [(Davidchack)](Davidchack5) to a group of atoms, which
models an interaction with an implicit background solvent. This command
performs Brownian dynamics (BD) via a technique that splits the
integration into a deterministic Hamiltonian part and the
Ornstein-Uhlenbeck process for noise and damping. The quaternion degrees
of freedom are updated though an evolution operator which performs a
rotation in quaternion space, preserves the quaternion norm and is akin
to [(Miller)](Miller5).

In terms of syntax this command has been closely modelled on the [fix
langevin](fix_langevin) and its *angmom* option. But it combines the
[fix nve](fix_nve) and the [fix langevin](fix_langevin) in one single
command. The main feature is improved stability over the standard
integrator, permitting slightly larger timestep sizes.

:::: note
::: title
Note
:::

Unlike the [fix langevin](fix_langevin) this command performs also time
integration of the translational and quaternion degrees of freedom.
::::

The total force on each atom will have the form:

$$\begin{aligned}
F =   & F_c + F_f + F_r \\
F_f = & - \frac{m}{\mathrm{damp}} v \\
F_r \propto & \sqrt{\frac{k_B T m}{dt~\mathrm{damp}}}
\end{aligned}$$

$F_c$ is the conservative force computed via the usual inter-particle
interactions ([pair_style](pair_style), [bond_style](bond_style), etc).
The $F_f$ and $F_r$ terms are implicitly taken into account by this fix
on a per-particle basis.

$F_f$ is a frictional drag or viscous damping term proportional to the
particle\'s velocity. The proportionality constant for each atom is
computed as $\frac{m}{\mathrm{damp}}$, where *m* is the mass of the
particle and damp is the damping factor specified by the user.

$F_r$ is a force due to solvent atoms at a temperature $T$ randomly
bumping into the particle. As derived from the fluctuation/dissipation
theorem, its magnitude as shown above is proportional to
$\sqrt{\frac{k_B T m}{dt~\mathrm{damp}}}$, where $k_B$ is the Boltzmann
constant, $T$ is the desired temperature, *m* is the mass of the
particle, *dt* is the timestep size, and damp is the damping factor.
Random numbers are used to randomize the direction and magnitude of this
force as described in [(Dunweg)](Dunweg5), where a uniform random number
is used (instead of a Gaussian random number) for speed.

------------------------------------------------------------------------

*Tstart* and *Tstop* have to be constant values, i.e. they cannot be
variables. If used together with the oxDNA force field for
coarse-grained simulation of DNA please note that T = 0.1 in oxDNA units
corresponds to T = 300 K.

The *damp* parameter is specified in time units and determines how
rapidly the temperature is relaxed. For example, a value of 0.03 means
to relax the temperature in a timespan of (roughly) 0.03 time units
$\tau$ (see the [units](units) command). The damp factor can be thought
of as inversely related to the viscosity of the solvent, i.e. a small
relaxation time implies a high-viscosity solvent and vice versa. See the
discussion about gamma and viscosity in the documentation for the [fix
viscous](fix_viscous) command for more details. Note that the value
78.9375 in the second example above corresponds to a diffusion constant,
which is about an order of magnitude larger than realistic ones. This
has been used to sample configurations faster in Brownian dynamics
simulations.

The random \# *seed* must be a positive integer. A Marsaglia random
number generator is used. Each processor uses the input seed to generate
its own unique seed and its own stream of random numbers. Thus the
dynamics of the system will not be identical on two runs on different
numbers of processors.

The keyword/value option has to be used in the following way:

This fix has to be used together with the *angmom* keyword. The
particles are always considered to have a finite size. The keyword
*angmom* enables thermostatting of the rotational degrees of freedom in
addition to the usual translational degrees of freedom.

The scale factor after the *angmom* keyword gives the ratio of the
rotational to the translational friction coefficient.

An example input file can be found in
examples/PACKAGES/cgdna/examples/duplex2/. Further details of the
implementation and stability of the integrators are contained in
[(Henrich)](Henrich5). The preprint version of the article can be found
[here](PDF/CG-DNA.pdf)\_.

------------------------------------------------------------------------

## Restrictions

These pair styles can only be used if LAMMPS was built with the
[CG-DNA](PKG-CG-DNA) package and the MOLECULE and ASPHERE package. See
the [Build package](Build_package) page for more info.

## Related commands

[fix nve](fix_nve), [fix langevin](fix_langevin), [fix
nve/dot](fix_nve_dot), [bond_style oxdna/fene](bond_oxdna), [bond_style
oxdna2/fene](bond_oxdna), [pair_style oxdna/excv](pair_oxdna),
[pair_style oxdna2/excv](pair_oxdna2)

## Default

none

------------------------------------------------------------------------

::: {#Davidchack5}
**(Davidchack)** R.L Davidchack, T.E. Ouldridge, M.V. Tretyakov. J.
Chem. Phys. 142, 144114 (2015).
:::

::: {#Miller5}
**(Miller)** T. F. Miller III, M. Eleftheriou, P. Pattnaik, A. Ndirango,
G. J. Martyna, J. Chem. Phys., 116, 8649-8659 (2002).
:::

::: {#Dunweg5}
**(Dunweg)** B. Dunweg, W. Paul, Int. J. Mod. Phys. C, 2, 817-27 (1991).
:::

::: {#Henrich5}
**(Henrich)** O. Henrich, Y. A. Gutierrez-Fosado, T. Curk, T. E.
Ouldridge, Eur. Phys. J. E 41, 57 (2018).
:::
