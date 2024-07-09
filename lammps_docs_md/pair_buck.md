# pair_style buck command

Accelerator Variants: *buck/gpu*, *buck/intel*, *buck/kk*, *buck/omp*

# pair_style buck/coul/cut command

Accelerator Variants: *buck/coul/cut/gpu*, *buck/coul/cut/intel*,
*buck/coul/cut/kk*, *buck/coul/cut/omp*

# pair_style buck/coul/long command

Accelerator Variants: *buck/coul/long/gpu*, *buck/coul/long/intel*,
*buck/coul/long/kk*, *buck/coul/long/omp*

# pair_style buck/coul/msm command

Accelerator Variants: *buck/coul/msm/omp*

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *buck* or *buck/coul/cut* or *buck/coul/long* or
    *buck/coul/msm*
-   args = list of arguments for a particular style

<!-- -->

    *buck* args = cutoff
      cutoff = global cutoff for Buckingham interactions (distance units)
    *buck/coul/cut* args = cutoff (cutoff2)
      cutoff = global cutoff for Buckingham (and Coulombic if only 1 arg) (distance units)
      cutoff2 = global cutoff for Coulombic (optional) (distance units)
    *buck/coul/long* args = cutoff (cutoff2)
      cutoff = global cutoff for Buckingham (and Coulombic if only 1 arg) (distance units)
      cutoff2 = global cutoff for Coulombic (optional) (distance units)
    *buck/coul/msm* args = cutoff (cutoff2)
      cutoff = global cutoff for Buckingham (and Coulombic if only 1 arg) (distance units)
      cutoff2 = global cutoff for Coulombic (optional) (distance units)

## Examples

``` LAMMPS
pair_style buck 2.5
pair_coeff * * 100.0 1.5 200.0
pair_coeff * * 100.0 1.5 200.0 3.0

pair_style buck/coul/cut 10.0
pair_style buck/coul/cut 10.0 8.0
pair_coeff * * 100.0 1.5 200.0
pair_coeff 1 1 100.0 1.5 200.0 9.0
pair_coeff 1 1 100.0 1.5 200.0 9.0 8.0

pair_style buck/coul/long 10.0
pair_style buck/coul/long 10.0 8.0
pair_coeff * * 100.0 1.5 200.0
pair_coeff 1 1 100.0 1.5 200.0 9.0

pair_style buck/coul/msm 10.0
pair_style buck/coul/msm 10.0 8.0
pair_coeff * * 100.0 1.5 200.0
pair_coeff 1 1 100.0 1.5 200.0 9.0
```

## Description

The *buck* style computes a Buckingham potential (exp/6 instead of
Lennard-Jones 12/6) given by

$$E = A e^{-r / \rho} - \frac{C}{r^6} \qquad r < r_c$$

where $\rho$ is an ionic-pair dependent length parameter, and $r_c$ is
the cutoff on both terms.

The styles with *coul/cut* or *coul/long* or *coul/msm* add a Coulombic
term as described for the [lj/cut](pair_lj) pair styles. For
*buck/coul/long* and *buc/coul/msm*, an additional damping factor is
applied to the Coulombic term so it can be used in conjunction with the
[kspace_style](kspace_style) command and its *ewald* or *pppm* or *msm*
option. The Coulombic cutoff specified for this style means that
pairwise interactions within this distance are computed directly;
interactions outside that distance are computed in reciprocal space.

If one cutoff is specified for the *born/coul/cut* and *born/coul/long*
and *born/coul/msm* styles, it is used for both the A,C and Coulombic
terms. If two cutoffs are specified, the first is used as the cutoff for
the A,C terms, and the second is the cutoff for the Coulombic term.

Note that these potentials are related to the [Born-Mayer-Huggins
potential](pair_born).

:::: note
::: title
Note
:::

For all these pair styles, the terms with A and C are always cutoff. The
additional Coulombic term can be cutoff or long-range (no cutoff)
depending on whether the style name includes coul/cut or coul/long or
coul/msm. If you wish the C/r\^6 term to be long-range (no cutoff), then
see the [pair_style buck/long/coul/long](pair_buck_long) command.
::::

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   A (energy units)
-   $\rho$ (distance units)
-   C (energy-distance\^6 units)
-   cutoff (distance units)
-   cutoff2 (distance units)

The second coefficient, $\rho$, must be greater than zero. The
coefficients A, $\rho$, and C can be written as analytical expressions
of $\epsilon$ and $\sigma$, in analogy to the Lennard-Jones potential
[(Khrapak)](Khrapak).

The latter 2 coefficients are optional. If not specified, the global A,C
and Coulombic cutoffs are used. If only one cutoff is specified, it is
used as the cutoff for both A,C and Coulombic interactions for this type
pair. If both coefficients are specified, they are used as the A,C and
Coulombic cutoffs for this type pair. You cannot specify 2 cutoffs for
style *buck*, since it has no Coulombic terms. For *buck/coul/long* only
the LJ cutoff can be specified since a Coulombic cutoff cannot be
specified for an individual I,J type pair. All type pairs use the same
global Coulombic cutoff specified in the pair_style command.

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

These styles support the [pair_modify](pair_modify) shift option for the
energy of the exp() and 1/r\^6 portion of the pair interaction.

The *buck/coul/long* pair style supports the [pair_modify](pair_modify)
table option to tabulate the short-range portion of the long-range
Coulombic interaction.

These styles support the pair_modify tail option for adding long-range
tail corrections to energy and pressure for the A,C terms in the pair
interaction.

These styles write their information to [binary restart files](restart),
so pair_style and pair_coeff commands do not need to be specified in an
input script that reads a restart file.

These styles can only be used via the *pair* keyword of the [run_style
respa](run_style) command. They do not support the *inner*, *middle*,
*outer* keywords.

## Restrictions

The *buck/coul/long* style is part of the KSPACE package. They are only
enabled if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff), [pair_style born](pair_born)

## Default

none

::: {#Khrapak}
**(Khrapak)** Khrapak, Chaudhuri, and Morfill, J Chem Phys, 134, 054120
(2011).
:::
