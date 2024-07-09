# pair_style lj/cut command

Accelerator Variants: *lj/cut/gpu*, *lj/cut/intel*, *lj/cut/kk*,
*lj/cut/opt*, *lj/cut/omp*

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *lj/cut*
-   args = list of arguments for a particular style

<!-- -->

    *lj/cut* args = cutoff
      cutoff = global cutoff for Lennard Jones interactions (distance units)

## Examples

``` LAMMPS
pair_style lj/cut 2.5
pair_coeff * * 1 1
pair_coeff 1 1 1 1.1 2.8
```

## Description

The *lj/cut* styles compute the standard 12/6 Lennard-Jones potential,
given by

$$E = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
    \left(\frac{\sigma}{r}\right)^6 \right]   \qquad r < r_c$$

$r_c$ is the cutoff.

See the [lj/cut/coul](pair_lj_cut_coul) styles to add a Coulombic
pairwise interaction and the [lj/cut/tip4p](pair_lj_cut_tip4p) styles to
add the TIP4P water model.

## Coefficients

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $\epsilon$ (energy units)
-   $\sigma$ (distance units)
-   LJ cutoff (distance units)

The last coefficient is optional. If not specified, the global LJ cutoff
specified in the pair_style command is used.

Note that $\sigma$ is defined in the LJ formula as the zero-crossing
distance for the potential, *not* as the energy minimum at
$r_0 = 2^{\frac{1}{6}} \sigma$. The \_[same]() potential function
becomes:

$$E = \epsilon \left[ \left(\frac{r_0}{r}\right)^{12} -
     2 \left(\frac{r_0}{r}\right)^6 \right]  \qquad r < r_c$$

When using the minimum as reference width. In the literature both
formulations are used, but the describe the same potential, only the
$\sigma$ value must be computed by $\sigma = r_0 /
2^{\frac{1}{6}}$ for use with LAMMPS, if this latter formulation is
used.

------------------------------------------------------------------------

A version of these styles with a soft core, *lj/cut/soft*, suitable for
use in free energy calculations, is part of the FEP package and is
documented with the [pair_style \*/soft](pair_fep_soft) styles.

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

For atom type pairs I,J and I != J, the epsilon and sigma coefficients
and cutoff distance for all of the lj/cut pair styles can be mixed. The
default mix value is *geometric*. See the \"pair_modify\" command for
details.

All of the *lj/cut* pair styles support the [pair_modify](pair_modify)
shift option for the energy of the Lennard-Jones portion of the pair
interaction.

All of the *lj/cut* pair styles support the [pair_modify](pair_modify)
tail option for adding a long-range tail correction to the energy and
pressure for the Lennard-Jones portion of the pair interaction.

All of the *lj/cut* pair styles write their information to [binary
restart files](restart), so pair_style and pair_coeff commands do not
need to be specified in an input script that reads a restart file.

The *lj/cut* pair styles support the use of the *inner*, *middle*, and
*outer* keywords of the [run_style respa](run_style) command, meaning
the pairwise forces can be partitioned by distance at different levels
of the rRESPA hierarchy. The other styles only support the *pair*
keyword of run_style respa. See the [run_style](run_style) command for
details.

------------------------------------------------------------------------

## Related commands

-   [pair_coeff](pair_coeff)
-   [pair_style lj/cut/coul/cut](pair_lj_cut_coul)
-   [pair_style lj/cut/coul/debye](pair_lj_cut_coul)
-   [pair_style lj/cut/coul/dsf](pair_lj_cut_coul)
-   [pair_style lj/cut/coul/long](pair_lj_cut_coul)
-   [pair_style lj/cut/coul/msm](pair_lj_cut_coul)
-   [pair_style lj/cut/coul/wolf](pair_lj_cut_coul)
-   [pair_style lj/cut/tip4p/cut](pair_lj_cut_tip4p)
-   [pair_style lj/cut/tip4p/long](pair_lj_cut_tip4p)

## Default

none
