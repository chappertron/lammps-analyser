# improper_style cvff command

Accelerator Variants: *cvff/intel*, *cvff/omp*

## Syntax

``` LAMMPS
improper_style cvff
```

## Examples

``` LAMMPS
improper_style cvff
improper_coeff 1 80.0 -1 4
```

## Description

The *cvff* improper style uses the potential

$$E = K [1 + d  \cos (n \phi) ]$$

where phi is the improper dihedral angle.

If the 4 atoms in an improper quadruplet (listed in the data file read
by the [read_data](read_data) command) are ordered I,J,K,L then the
improper dihedral angle is between the plane of I,J,K and the plane of
J,K,L. Note that because this is effectively a dihedral angle, the
formula for this improper style is the same as for [dihedral_style
harmonic](dihedral_harmonic).

Note that defining 4 atoms to interact in this way, does not mean that
bonds necessarily exist between I-J, J-K, or K-L, as they would in a
linear dihedral. Normally, the bonds I-J, I-K, I-L would exist for an
improper to be defined between the 4 atoms.

The following coefficients must be defined for each improper type via
the [improper_coeff](improper_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands:

-   $K$ (energy)
-   $d$ (+1 or -1)
-   $n$ (0,1,2,3,4,6)

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

This improper style can only be used if LAMMPS was built with the
MOLECULE package. See the [Build package](Build_package) doc page for
more info.

## Related commands

[improper_coeff](improper_coeff)

## Default

none
