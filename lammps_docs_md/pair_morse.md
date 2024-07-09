# pair_style morse command

Accelerator Variants: *morse/gpu*, *morse/omp*, *morse/opt*, *morse/kk*

# pair_style morse/smooth/linear command

Accelerator Variants: *morse/smooth/linear/omp*

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *morse* or *morse/smooth/linear* or *morse/soft*
-   args = list of arguments for a particular style

<!-- -->

    *morse* args = cutoff
      cutoff = global cutoff for Morse interactions (distance units)
    *morse/smooth/linear* args = cutoff
      cutoff = global cutoff for Morse interactions (distance units)

## Examples

``` LAMMPS
pair_style morse 2.5
pair_style morse/smooth/linear 2.5
pair_coeff * * 100.0 2.0 1.5
pair_coeff 1 1 100.0 2.0 1.5 3.0
```

## Description

Style *morse* computes pairwise interactions with the formula

$$E = D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right]
    \qquad r < r_c$$

Rc is the cutoff.

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $D_0$ (energy units)
-   $\alpha$ (1/distance units)
-   $r_0$ (distance units)
-   cutoff (distance units)

The last coefficient is optional. If not specified, the global morse
cutoff is used.

------------------------------------------------------------------------

The *morse/smooth/linear* variant is similar to the lj/smooth/linear
variant in that it adds to the potential a shift and a linear term so
that both, potential energy and force, go to zero at the cut-off:

$$\begin{aligned}
\phi\left(r\right) & =  D_0 \left[ e^{- 2 \alpha (r - r_0)} - 2 e^{- \alpha (r - r_0)} \right] \qquad r < r_c \\
E\left(r\right) & =  \phi\left(r\right)  - \phi\left(R_c\right) - \left(r - R_c\right) \left.\frac{d\phi}{d r} \right|_{r=R_c}       \qquad r < R_c
\end{aligned}$$

The syntax of the pair_style and pair_coeff commands are the same for
the *morse* and *morse/smooth/linear* styles.

------------------------------------------------------------------------

A version of the *morse* style with a soft core, *morse/soft*, suitable
for use in free energy calculations, is part of the FEP package and is
documented with the [pair_style \*/soft](pair_fep_soft) styles. The
version with soft core is only available if LAMMPS was built with that
package. See the [Build package](Build_package) page for more info.

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

None of these pair styles support mixing. Thus, coefficients for all I,J
pairs must be specified explicitly.

All of these pair styles support the [pair_modify](pair_modify) shift
option for the energy of the pair interaction.

The [pair_modify](pair_modify) table options are not relevant for the
Morse pair styles.

None of these pair styles support the [pair_modify](pair_modify) tail
option for adding long-range tail corrections to energy and pressure.

All of these pair styles write their information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
[run_style respa](run_style) command. They do not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

The *morse/smooth/linear* pair style is only enabled if LAMMPS was built
with the EXTRA-PAIR package. See the [Build package](Build_package) page
for more info.

## Related commands

[pair_coeff](pair_coeff), [pair_style \*/soft](pair_fep_soft)

## Default

none
