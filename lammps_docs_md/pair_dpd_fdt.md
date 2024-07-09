# pair_style dpd/fdt command

# pair_style dpd/fdt/energy command

Accelerator Variants: *dpd/fdt/energy/kk*

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *dpd/fdt* or *dpd/fdt/energy*
-   args = list of arguments for a particular style

<!-- -->

    *dpd/fdt* args = T cutoff seed
      T = temperature (temperature units)
      cutoff = global cutoff for DPD interactions (distance units)
      seed = random # seed (positive integer)
    *dpd/fdt/energy* args = cutoff seed
      cutoff = global cutoff for DPD interactions (distance units)
      seed = random # seed (positive integer)

## Examples

``` LAMMPS
pair_style dpd/fdt 300.0 2.5 34387
pair_coeff * * 3.0 1.0 2.5

pair_style dpd/fdt/energy 2.5 34387
pair_coeff * * 3.0 1.0 0.1 2.5
```

## Description

Styles *dpd/fdt* and *dpd/fdt/energy* compute the force for dissipative
particle dynamics (DPD) simulations. The *dpd/fdt* style is used to
perform DPD simulations under isothermal and isobaric conditions, while
the *dpd/fdt/energy* style is used to perform DPD simulations under
isoenergetic and isoenthalpic conditions (see [(Lisal)](Lisal3)). For
DPD simulations in general, the force on atom I due to atom J is given
as a sum of 3 terms

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
velocities of the two atoms, $\vec{v}_i -
\vec{v}_j$, $\alpha$ is a Gaussian random number with zero mean and unit
variance, *dt* is the timestep size, and $w(r)$ is a weighting factor
that varies between 0 and 1, $r_c$ is the pairwise cutoff. Note that
alternative definitions of the weighting function exist, but would have
to be implemented as a separate pair style command.

For style *dpd/fdt*, the fluctuation-dissipation theorem defines
$\gamma$ to be set equal to $\sigma^2/(2 T)$, where *T* is the set point
temperature specified as a pair style parameter in the above examples.
The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   A (force units)
-   $\sigma$ (force\*time\^(1/2) units)
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global DPD
cutoff is used.

Style *dpd/fdt/energy* is used to perform DPD simulations under
isoenergetic and isoenthalpic conditions. The fluctuation-dissipation
theorem defines $\gamma$ to be set equal to $\sigma^2/(2
\theta)$, where $\theta$ is the average internal temperature for the
pair. The particle internal temperature is related to the particle
internal energy through a mesoparticle equation of state (see [fix
eos](fix)). The differential internal conductive and mechanical energies
are computed within style *dpd/fdt/energy* as:

$$\begin{aligned}
du_{i}^{cond}  = & \kappa_{ij}(\frac{1}{\theta_{i}}-\frac{1}{\theta_{j}})\omega_{ij}^{2} + \alpha_{ij}\omega_{ij}\zeta_{ij}^{q}(\Delta{t})^{-1/2} \\
du_{i}^{mech}  = & -\frac{1}{2}\gamma_{ij}\omega_{ij}^{2}(\frac{\vec{r}_{ij}}{r_{ij}}\bullet\vec{v}_{ij})^{2} -
\frac{\sigma^{2}_{ij}}{4}(\frac{1}{m_{i}}+\frac{1}{m_{j}})\omega_{ij}^{2} -
\frac{1}{2}\sigma_{ij}\omega_{ij}(\frac{\vec{r}_{ij}}{r_{ij}}\bullet\vec{v}_{ij})\zeta_{ij}(\Delta{t})^{-1/2}
\end{aligned}$$

where

$$\begin{aligned}
\alpha_{ij}^{2}  = & 2k_{B}\kappa_{ij} \\
\sigma^{2}_{ij}  = & 2\gamma_{ij}k_{B}\Theta_{ij} \\
\Theta_{ij}^{-1}  = & \frac{1}{2}(\frac{1}{\theta_{i}}+\frac{1}{\theta_{j}})
\end{aligned}$$

$\zeta_ij^q$ is a second Gaussian random number with zero mean and unit
variance that is used to compute the internal conductive energy. The
fluctuation-dissipation theorem defines $alpha^2$ to be set equal to
$2k_B\kappa$, where $\kappa$ is the mesoparticle thermal conductivity
parameter. The following coefficients must be defined for each pair of
atoms types via the [pair_coeff](pair_coeff) command as in the examples
above, or in the data file or restart files read by the
[read_data](read_data) or [read_restart](read_restart) commands:

-   A (force units)
-   $\sigma$ (force\*time\^(1/2) units)
-   $\kappa$ (energy\*temperature/time units)
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global DPD
cutoff is used.

The pairwise energy associated with styles *dpd/fdt* and
*dpd/fdt/energy* is only due to the conservative force term $F^C$, and
is shifted to be zero at the cutoff distance $r_c$. The pairwise virial
is calculated using only the conservative term.

The forces computed through the *dpd/fdt* and *dpd/fdt/energy* styles
can be integrated with the velocity-Verlet integration scheme or the
Shardlow splitting integration scheme described by [(Lisal)](Lisal3). In
the cases when these pair styles are combined with the [fix
shardlow](fix_shardlow), these pair styles differ from the other dpd
styles in that the dissipative and random forces are split from the
force calculation and are not computed within the pair style. Thus, only
the conservative force is computed by the pair style, while the
stochastic integration of the dissipative and random forces are handled
through the Shardlow splitting algorithm approach. The Shardlow
splitting algorithm is advantageous, especially when performing DPD
under isoenergetic conditions, as it allows significantly larger
timesteps to be taken.

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

## Restrictions

These commands are part of the DPD-REACT package. They are only enabled
if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

Pair styles *dpd/fdt* and *dpd/fdt/energy* require use of the
[comm_modify vel yes](comm_modify) option so that velocities are stored
by ghost atoms.

Pair style *dpd/fdt/energy* requires [atom_style dpd](atom_style) to be
used in order to properly account for the particle internal energies and
temperatures.

## Related commands

[pair_coeff](pair_coeff), [fix shardlow](fix_shardlow)

## Default

none

------------------------------------------------------------------------

::: {#Lisal3}
**(Lisal)** M. Lisal, J.K. Brennan, J. Bonet Avalos, \"Dissipative
particle dynamics at isothermal, isobaric, isoenergetic, and
isoenthalpic conditions using Shardlow-like splitting algorithms.\", J.
Chem. Phys., 135, 204105 (2011).
:::
