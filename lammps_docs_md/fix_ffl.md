# fix ffl command

## Syntax

    fix ID id-group ffl tau Tstart Tstop seed [flip-type]

-   ID, group-ID are documented in [fix](fix) command

-   ffl = style name of this fix command

-   tau = thermostat parameter (positive real)

-   Tstart, Tstop = temperature ramp during the run

-   seed = random number seed to use for generating noise (positive
    integer)

-   one more value may be appended

        flip-type  = determines the flipping type, can be chosen between rescale - no_flip - hard - soft, if no flip type is given, rescale will be chosen by default

## Examples

``` LAMMPS
fix 3 boundary ffl 10 300 300 31415
fix 1 all ffl 100 500 500 9265 soft
```

## Description

Apply a Fast-Forward Langevin Equation (FFL) thermostat as described in
[(Hijazi)](Hijazi). Contrary to [fix langevin](fix_langevin), this fix
performs both thermostatting and evolution of the Hamiltonian equations
of motion, so it should not be used together with [fix nve](fix_nve) \--
at least not on the same atom groups.

The time-evolution of a single particle undergoing Langevin dynamics is
described by the equations

$$\frac {dq}{dt} = \frac{p}{m},$$

$$\frac {dp}{dt} = -\gamma p + W + F,$$

where $F$ is the physical force, $\gamma$ is the friction coefficient,
and $W$ is a Gaussian random force.

The friction coefficient is the inverse of the thermostat parameter :
$\gamma = 1/\tau$, with $\tau$ the thermostat parameter *tau*. The
thermostat parameter is given in the time units, $\gamma$ is in inverse
time units.

Equilibrium sampling a temperature T is obtained by specifying the
target value as the *Tstart* and *Tstop* arguments, so that the internal
constants depending on the temperature are computed automatically.

The random number *seed* must be a positive integer. A Marsaglia random
number generator is used. Each processor uses the input seed to generate
its own unique seed and its own stream of random numbers. Thus the
dynamics of the system will not be identical on two runs on different
numbers of processors.

The flipping type *flip-type* can be chosen between 4 types described in
[(Hijazi)](Hijazi). The flipping operation occurs during the
thermostatting step and it flips the momenta of the atoms. If no_flip is
chosen, no flip will be executed and the integration will be the same as
a standard Langevin thermostat [(Bussi)](Bussi3). The other flipping
types are : rescale - hard - soft.

## Restart, fix_modify, output, run start/stop, minimize info

The instantaneous values of the extended variables are written to
[binary restart files](restart). Because the state of the random number
generator is not saved in restart files, this means you cannot do
\"exact\" restarts with this fix, where the simulation continues on the
same as if no restart had taken place. However, in a statistical sense,
a restarted simulation should produce the same behavior. Note however
that you should use a different seed each time you restart, otherwise
the same sequence of random numbers will be used each time, which might
lead to stochastic synchronization and subtle artifacts in the sampling.

The cumulative energy change in the system imposed by this fix is
included in the [thermodynamic output](thermo_style) keywords *ecouple*
and *econserve*. See the [thermo_style](thermo_style) doc page for
details.

This fix computes a global scalar which can be accessed by various
[output commands](Howto_output). The scalar is the same cumulative
energy change due to this fix described in the previous paragraph. The
scalar value calculated by this fix is \"extensive\".

This fix can ramp its target temperature over multiple runs, using the
*start* and *stop* keywords of the [run](run) command. See the
[run](run) command for details of how to do this.

This fix is not invoked during [energy minimization](minimize).

## Restrictions

In order to perform constant-pressure simulations please use [fix
press/berendsen](fix_press_berendsen), rather than [fix npt](fix_nh), to
avoid duplicate integration of the equations of motion.

This fix is part of the EXTRA-FIX package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

## Related commands

[fix nvt](fix_nh), [fix temp/rescale](fix_temp_rescale), [fix
viscous](fix_viscous), [fix nvt](fix_nh), [pair_style
dpd/tstat](pair_dpd), [fix gld](fix_gld), [fix gle](fix_gle)

------------------------------------------------------------------------

:::: {#Hijazi}
::: {#Bussi3}
**(Hijazi)** M. Hijazi, D. M. Wilkins, M. Ceriotti, J. Chem. Phys. 148,
184109 (2018)
:::
::::

**(Bussi)** G. Bussi, M. Parrinello, Phs. Rev. E 75, 056707 (2007)
