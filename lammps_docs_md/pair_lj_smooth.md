# pair_style lj/smooth command

Accelerator Variants: *lj/smooth/gpu*, *lj/smooth/omp*

## Syntax

``` LAMMPS
pair_style lj/smooth Rin Rc
```

-   Rin = inner cutoff beyond which force smoothing will be applied
    (distance units)
-   Rc = outer cutoff for lj/smooth interactions (distance units)

## Examples

``` LAMMPS
pair_style lj/smooth 8.0 10.0
pair_coeff * * 10.0 1.5
pair_coeff 1 1 20.0 1.3 7.0 9.0
```

## Description

Style *lj/smooth* computes a LJ interaction with a force smoothing
applied between the inner and outer cutoff.

$$\begin{aligned}
E & =  4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} -
                      \left(\frac{\sigma}{r}\right)^6 \right]
                      \qquad r < r_{in} \\
F & =  C_1 + C_2 (r - r_{in}) + C_3 (r - r_{in})^2 + C_4 (r - r_{in})^3
                    \qquad r_{in} < r < r_c
\end{aligned}$$

The polynomial coefficients C1, C2, C3, C4 are computed by LAMMPS to
cause the force to vary smoothly from the inner cutoff $r_{in}$ to the
outer cutoff $r_c$.

At the inner cutoff the force and its first derivative will match the
non-smoothed LJ formula. At the outer cutoff the force and its first
derivative will be 0.0. The inner cutoff cannot be 0.0.

:::: note
::: title
Note
:::

this force smoothing causes the energy to be discontinuous both in its
values and first derivative. This can lead to poor energy conservation
and may require the use of a thermostat. Plot the energy and force
resulting from this formula via the [pair_write](pair_write) command to
see the effect.
::::

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   $\epsilon$ (energy units)
-   $\sigma$ (distance units)
-   $r_{in}$ (distance units)
-   $r_c$ (distance units)

The last 2 coefficients are optional inner and outer cutoffs. If not
specified, the global values for $r_{in}$ and $r_c$ are used.

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

For atom type pairs I,J and I != J, the epsilon, sigma, Rin coefficients
and the cutoff distance for this pair style can be mixed. Rin is a
cutoff value and is mixed like the cutoff. The other coefficients are
mixed according to the pair_modify mix option. The default mix value is
*geometric*. See the \"pair_modify\" command for details.

This pair style supports the [pair_modify](pair_modify) shift option for
the energy of the pair interaction.

The [pair_modify](pair_modify) table option is not relevant for this
pair style.

This pair style does not support the [pair_modify](pair_modify) tail
option for adding long-range tail corrections to energy and pressure,
since the energy of the pair interaction is smoothed to 0.0 at the
cutoff.

This pair style writes its information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

------------------------------------------------------------------------

## Restrictions

This pair style is part of the EXTRA-PAIR package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff), [pair lj/smooth/linear](pair_lj_smooth_linear)

## Default

none
