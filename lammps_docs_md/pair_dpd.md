# pair_style dpd command

Accelerator Variants: *dpd/gpu*, *dpd/intel*, *dpd/kk*, *dpd/omp*

# pair_style dpd/tstat command

Accelerator Variants: *dpd/tstat/gpu*, *dpd/tstat/kk*, *dpd/tstat/omp*

## Syntax

``` LAMMPS
pair_style dpd T cutoff seed
pair_style dpd/tstat Tstart Tstop cutoff seed
```

-   T = temperature (temperature units) (dpd only)
-   Tstart,Tstop = desired temperature at start/end of run (temperature
    units) (dpd/tstat only)
-   cutoff = global cutoff for DPD interactions (distance units)
-   seed = random \# seed (positive integer)

## Examples

``` LAMMPS
pair_style dpd 1.0 2.5 34387
pair_coeff * * 3.0 1.0
pair_coeff 1 1 3.0 1.0 1.0

pair_style hybrid/overlay lj/cut 2.5 dpd/tstat 1.0 1.0 2.5 34387
pair_coeff * * lj/cut 1.0 1.0
pair_coeff * * dpd/tstat 1.0
```

## Description

Style *dpd* computes a force field for dissipative particle dynamics
(DPD) following the exposition in [(Groot)](Groot1).

Style *dpd/tstat* invokes a DPD thermostat on pairwise interactions,
which is equivalent to the non-conservative portion of the DPD force
field. This pairwise thermostat can be used in conjunction with any
[pair style](pair_style), and instead of per-particle thermostats like
[fix langevin](fix_langevin) or ensemble thermostats like Nose Hoover as
implemented by [fix nvt](fix_nh). To use *dpd/tstat* as a thermostat for
another pair style, use the [pair_style hybrid/overlay](pair_hybrid)
command to compute both the desired pair interaction and the thermostat
for each pair of particles.

For style *dpd*, the force on atom I due to atom J is given as a sum of
3 terms

$$\begin{aligned}
\vec{f}  = & (F^C + F^D + F^R) \hat{r_{ij}} \qquad \qquad r < r_c \\
F^C      = & A w(r) \\
F^D      = & - \gamma w^2(r) (\hat{r_{ij}} \bullet \vec{v}_{ij}) \\
F^R      = & \sigma w(r) \alpha (\Delta t)^{-1/2} \\
w(r)     = & 1 - \frac{r}{r_c}
\end{aligned}$$

where $F^C$ is a conservative force, $F^D$ is a dissipative force, and
$F^R$ is a random force. $\hat{r_{ij}}$ is a unit vector in the
direction $r_i - r_j$, $\vec{v}_{ij}$ is the vector difference in
velocities of the two atoms $\vec{v}_i -
\vec{v}_j$, $\alpha$ is a Gaussian random number with zero mean and unit
variance, *dt* is the timestep size, and $w(r)$ is a weighting factor
that varies between 0 and 1. $r_c$ is the pairwise cutoff. $\sigma$ is
set equal to $\sqrt{2 k_B T
\gamma}$, where $k_B$ is the Boltzmann constant and *T* is the
temperature parameter in the pair_style command.

For style *dpd/tstat*, the force on atom I due to atom J is the same as
the above equation, except that the conservative $F^C$ term is dropped.
Also, during the run, *T* is set each timestep to a ramped value from
*Tstart* to *Tstop*.

For style *dpd*, the pairwise energy associated with style *dpd* is only
due to the conservative force term $F^C$, and is shifted to be zero at
the cutoff distance $r_c$. The pairwise virial is calculated using all 3
terms. For style *dpd/tstat* there is no pairwise energy, but the last
two terms of the formula make a contribution to the virial.

For style *dpd*, the following coefficients must be defined for each
pair of atoms types via the [pair_coeff](pair_coeff) command as in the
examples above, or in the data file or restart files read by the
[read_data](read_data) or [read_restart](read_restart) commands:

-   A (force units)
-   $\gamma$ (force/velocity units)
-   cutoff (distance units)

The cutoff coefficient is optional. If not specified, the global DPD
cutoff is used. Note that sigma is set equal to sqrt(2 T gamma), where T
is the temperature set by the [pair_style](pair_style) command so it
does not need to be specified.

For style *dpd/tstat*, the coefficients defined for each pair of atoms
types via the [pair_coeff](pair_coeff) command are:

-   $\gamma$ (force/velocity units)
-   cutoff (distance units)

The cutoff coefficient is optional.

The GPU-accelerated versions of these styles are implemented based on
the work of [(Afshar)](Afshar) and [(Phillips)](Phillips).

:::: note
::: title
Note
:::

If you are modeling DPD polymer chains, you may want to use the
[pair_style srp](pair_srp) command in conjunction with these pair
styles. It is a soft segmental repulsive potential (SRP) that can
prevent DPD polymer chains from crossing each other.
::::

:::: note
::: title
Note
:::

The virial calculation for pressure when using these pair styles
includes all the components of force listed above, including the random
force. Since the random force depends on random numbers, everything that
changes the order of atoms in the neighbor list (e.g. different number
of MPI ranks or a different neighbor list skin distance) will also
change the sequence in which the random numbers are applied and thus the
individual forces and therefore also the virial/pressure.
::::

:::: note
::: title
Note
:::

For more consistent time integration and force computation you may
consider using [fix mvv/dpd](fix_mvv_dpd) instead of [fix nve](fix_nve).
::::

------------------------------------------------------------------------

Styles with a *gpu*, *intel*, *kk*, *omp*, or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the [Accelerator packages](Speed_packages)
page. The accelerated styles take the same arguments and should produce
the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, INTEL, KOKKOS, OPENMP, and
OPT packages, respectively. They are only enabled if LAMMPS was built
with those packages. See the [Build package](Build_package) page for
more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the [-suffix command-line
switch](Run_options) when you invoke LAMMPS, or you can use the
[suffix](suffix) command in your input script.

See the [Accelerator packages](Speed_packages) page for more
instructions on how to use the accelerated styles effectively.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

These pair styles do not support mixing. Thus, coefficients for all I,J
pairs must be specified explicitly.

These pair styles do not support the [pair_modify](pair_modify) shift
option for the energy of the pair interaction. Note that as discussed
above, the energy due to the conservative $F^C$ term is already shifted
to be 0.0 at the cutoff distance $r_c$.

The [pair_modify](pair_modify) table option is not relevant for these
pair styles.

These pair styles do not support the [pair_modify](pair_modify) tail
option for adding long-range tail corrections to energy and pressure.

These pair styles write their information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file. Note that the
user-specified random number seed is stored in the restart file, so when
a simulation is restarted, each processor will re-initialize its random
number generator the same way it did initially. This means the random
forces will be random, but will not be the same as they would have been
if the original simulation had continued past the restart time.

These pair styles can only be used via the *pair* keyword of the
[run_style respa](run_style) command. They do not support the *inner*,
*middle*, *outer* keywords.

The *dpd/tstat* style can ramp its target temperature over multiple
runs, using the *start* and *stop* keywords of the [run](run) command.
See the [run](run) command for details of how to do this.

------------------------------------------------------------------------

## Restrictions

These styles are part of the DPD-BASIC package. They are only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

The default frequency for rebuilding neighbor lists is every 10 steps
(see the [neigh_modify](neigh_modify) command). This may be too
infrequent for style *dpd* simulations since particles move rapidly and
can overlap by large amounts. If this setting yields a non-zero number
of \"dangerous\" reneighborings (printed at the end of a simulation),
you should experiment with forcing reneighboring more often and see if
system energies/trajectories change.

These pair styles requires you to use the [comm_modify vel
yes](comm_modify) command so that velocities are stored by ghost atoms.

These pair styles will not restart exactly when using the
[read_restart](read_restart) command, though they should provide
statistically similar results. This is because the forces they compute
depend on atom velocities. See the [read_restart](read_restart) command
for more details.

## Related commands

[pair_style dpd/ext](pair_dpd_ext), [pair_coeff](pair_coeff), [fix
nvt](fix_nh), [fix langevin](fix_langevin), [pair_style srp](pair_srp),
[fix mvv/dpd](fix_mvv_dpd).

## Default

none

------------------------------------------------------------------------

::: {#Groot1}
**(Groot)** Groot and Warren, J Chem Phys, 107, 4423-35 (1997).
:::

::: {#Afshar}
**(Afshar)** Afshar, F. Schmid, A. Pishevar, S. Worley, Comput Phys
Comm, 184, 1119-1128 (2013).
:::

::: {#Phillips}
**(Phillips)** C. L. Phillips, J. A. Anderson, S. C. Glotzer, Comput
Phys Comm, 230, 7191-7201 (2011).
:::
