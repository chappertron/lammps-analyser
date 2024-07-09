# pair_style coul/cut/dielectric command

# pair_style coul/long/dielectric command

# pair_style lj/cut/coul/cut/dielectric command

Accelerator Variants: *lj/cut/coul/cut/dielectric/omp*

# pair_style lj/cut/coul/debye/dielectric command

Accelerator Variants: *lj/cut/coul/debye/dielectric/omp*

# pair_style lj/cut/coul/long/dielectric command

Accelerator Variants: *lj/cut/coul/long/dielectric/omp*

# pair_style lj/cut/coul/msm/dielectric command

# pair_style lj/long/coul/long/dielectric command

## Syntax

``` LAMMPS
pair_style style args
```

-   style = *lj/cut/coul/cut/dielectric* or
    *lj/cut/coul/long/dielectric* or *lj/cut/coul/msm/dielectric* or
    *lj/long/coul/msm/dielectric*
-   args = list of arguments for a particular style

## Examples

``` LAMMPS
pair_style coul/cut/dielectric 10.0
pair_coeff * *
pair_coeff 1 1 9.0

pair_style lj/cut/coul/cut/dielectric 10.0
pair_style lj/cut/coul/cut/dielectric 10.0 8.0
pair_coeff * * 100.0 3.0
pair_coeff 1 1 100.0 3.5 9.0

pair_style lj/cut/coul/long/dielectric 10.0
pair_style lj/cut/coul/long/dielectric 10.0 8.0
pair_coeff * * 100.0 3.0
pair_coeff 1 1 100.0 3.5 9.0
```

Used in input scripts:

> examples/PACKAGES/dielectric/in.confined
>     examples/PACKAGES/dielectric/in.nopbc

## Description

All these pair styles are derived from the corresponding pair styles
without the *dielectric* suffix. In addition to computing atom forces
and energies, these pair styles compute the electric field vector at
each atom, which are intended to be used by the [fix
polarize](fix_polarize) commands to compute induced charges at
interfaces between two regions of different dielectric constant.

These pair styles should be used with [atom_style
dielectric](atom_style).

The styles lj/cut/coul/long/dielectric, lj/cut/coul/msm/dielectric, and
lj/long/coul/long/dielectric should be used with their kspace style
counterparts, namely, pppm/dielectric, pppm/disp/dielectric, and
msm/dielectric, respectively.

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
and cutoff distances for this pair style can be mixed. The default mix
algorithm is *geometric*. See the [pair_modify](pair_modify)\" command
for details.

The [pair_modify](pair_modify) table option is not relevant for this
pair style.

These pair styles write its information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

These pair styles can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

## Restrictions

These styles are part of the DIELECTRIC package. They are only enabled
if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_coeff](pair_coeff), [fix polarize](fix_polarize),
[read_data](read_data)

## Default

none
