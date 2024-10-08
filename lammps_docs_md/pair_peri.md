# pair_style peri/pmb command

Accelerator Variants: *peri/pmb/omp*

# pair_style peri/lps command

Accelerator Variants: *peri/lps/omp*

# pair_style peri/ves command

# pair_style peri/eps command

## Syntax

``` LAMMPS
pair_style style
```

-   style = *peri/pmb* or *peri/lps* or *peri/ves* or *peri/eps*

## Examples

``` LAMMPS
pair_style peri/pmb
pair_coeff * * 1.6863e22 0.0015001 0.0005 0.25

pair_style peri/lps
pair_coeff * * 14.9e9 14.9e9 0.0015001 0.0005 0.25

pair_style peri/ves
pair_coeff * * 14.9e9 14.9e9 0.0015001 0.0005 0.25 0.5 0.001

pair_style peri/eps
pair_coeff * * 14.9e9 14.9e9 0.0015001 0.0005 0.25 118.43
```

## Description

The peridynamic pair styles implement material models that can be used
at the mesoscopic and macroscopic scales. See [this
document](PDF/PDLammps_overview.pdf)\_ for an overview of LAMMPS
commands for Peridynamics modeling.

Style *peri/pmb* implements the Peridynamic bond-based prototype
microelastic brittle (PMB) model.

Style *peri/lps* implements the Peridynamic state-based linear
peridynamic solid (LPS) model.

Style *peri/ves* implements the Peridynamic state-based linear
peridynamic viscoelastic solid (VES) model.

Style *peri/eps* implements the Peridynamic state-based elastic-plastic
solid (EPS) model.

The canonical papers on Peridynamics are [(Silling 2000)](Silling2000)
and [(Silling 2007)](Silling2007). The implementation of Peridynamics in
LAMMPS is described in [(Parks)](Parks). Also see the [Peridynamics
Howto](Howto_peri) for more details about its implementation.

The peridynamic VES and EPS models in PDLAMMPS were implemented by R.
Rahman and J. T. Foster at University of Texas at San Antonio. The
original VES formulation is described in \"(Mitchell2011)\" and the
original EPS formulation is in \"(Mitchell2011a)\". Additional PDF docs
that describe the VES and EPS implementations are include in the LAMMPS
distribution in [doc/PDF/PDLammps_VES.pdf](PDF/PDLammps_VES.pdf)\_ and
[doc/PDF/PDLammps_EPS.pdf](PDF/PDLammps_EPS.pdf)\_. For questions
regarding the VES and EPS models in LAMMPS you can contact R. Rahman
(rezwanur.rahman at utsa.edu).

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below.

For the *peri/pmb* style:

-   c (energy/distance/volume\^2 units)
-   horizon (distance units)
-   s00 (unitless)
-   $\alpha$ (unitless)

C is the effectively a spring constant for Peridynamic bonds, the
horizon is a cutoff distance for truncating interactions, and s00 and
$\alpha$ are used as a bond breaking criteria. The units of c are such
that c/distance = stiffness/volume\^2, where stiffness is
energy/distance\^2 and volume is distance\^3. See the users guide for
more details.

For the *peri/lps* style:

-   K (force/area units)
-   G (force/area units)
-   horizon (distance units)
-   s00 (unitless)
-   $\alpha$ (unitless)

K is the bulk modulus and G is the shear modulus. The horizon is a
cutoff distance for truncating interactions, and s00 and $\alpha$ are
used as a bond breaking criteria. See the users guide for more details.

For the *peri/ves* style:

-   K (force/area units)
-   G (force/area units)
-   horizon (distance units)
-   s00 (unitless)
-   $\alpha$ (unitless)
-   m_lambdai (unitless)
-   m_taubi (unitless)

K is the bulk modulus and G is the shear modulus. The horizon is a
cutoff distance for truncating interactions, and s00 and $\alpha$ are
used as a bond breaking criteria. m_lambdai and m_taubi are the
viscoelastic relaxation parameter and time constant, respectively.
m_lambdai varies within zero to one. For very small values of m_lambdai
the viscoelastic model responds very similar to a linear elastic model.
For details please see the description in \"(Mitchell2011)\".

For the *peri/eps* style:

-   K (force/area units)
-   G (force/area units)
-   horizon (distance units)
-   s00 (unitless)
-   $\alpha$ (unitless)
-   m_yield_stress (force/area units)

K is the bulk modulus and G is the shear modulus. The horizon is a
cutoff distance and s00 and $\alpha$ are used as a bond breaking
criteria. m_yield_stress is the yield stress of the material. For
details please see the description in \"(Mitchell2011a)\".

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
option.

The [pair_modify](pair_modify) table and tail options are not relevant
for these pair styles.

These pair styles write their information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
[run_style respa](run_style) command. They do not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

All of these styles are part of the PERI package. They are only enabled
if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff)

## Default

none

------------------------------------------------------------------------

::: {#Parks}
**(Parks)** Parks, Lehoucq, Plimpton, Silling, Comp Phys Comm, 179(11),
777-783 (2008).
:::

::: {#Silling2000}
**(Silling 2000)** Silling, J Mech Phys Solids, 48, 175-209 (2000).
:::

::: {#Silling2007}
**(Silling 2007)** Silling, Epton, Weckner, Xu, Askari, J Elasticity,
88, 151-184 (2007).
:::

::: {#Mitchell2011}
**(Mitchell2011)** Mitchell. A non-local, ordinary-state-based
viscoelasticity model for peridynamics. Sandia National Lab Report,
8064:1-28 (2011).
:::

::: {#Mitchell2011a}
**(Mitchell2011a)** Mitchell. A Nonlocal, Ordinary, State-Based
Plasticity Model for Peridynamics. Sandia National Lab Report, 3166:1-34
(2011).
:::
