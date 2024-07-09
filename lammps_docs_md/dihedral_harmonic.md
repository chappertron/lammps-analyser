# dihedral_style harmonic command

Accelerator Variants: *harmonic/intel*, *harmonic/kk*, *harmonic/omp*

## Syntax

``` LAMMPS
dihedral_style harmonic
```

## Examples

``` LAMMPS
dihedral_style harmonic
dihedral_coeff 1 80.0 1 2
```

## Description

The *harmonic* dihedral style uses the potential

$$E = K [ 1 + d  \cos (n \phi) ]$$

The following coefficients must be defined for each dihedral type via
the [dihedral_coeff](dihedral_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy)
-   $d$ (+1 or -1)
-   $n$ (integer \>= 0)

:::: note
::: title
Note
:::

Here are important points to take note of when defining LAMMPS dihedral
coefficients for the harmonic style, so that they are compatible with
how harmonic dihedrals are defined by other force fields:
::::

-   The LAMMPS convention is that the trans position = 180 degrees,
    while in some force fields trans = 0 degrees.
-   Some force fields reverse the sign convention on $d$.
-   Some force fields let $n$ be positive or negative which corresponds
    to $d = 1$ or $d = -1$ for the harmonic style.

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

This dihedral style can only be used if LAMMPS was built with the
MOLECULE package. See the [Build package](Build_package) doc page for
more info.

## Related commands

[dihedral_coeff](dihedral_coeff)

## Default

none
