# fix freeze command

Accelerator Variants: *freeze/kk*

## Syntax

    fix ID group-ID freeze

-   ID, group-ID are documented in [fix](fix) command
-   freeze = style name of this fix command

## Examples

``` LAMMPS
fix 2 bottom freeze
```

## Description

Zero out the force and torque on a granular particle. This is useful for
preventing certain particles from moving in a simulation. The [granular
pair styles](pair_gran) also detect if this fix has been defined and
compute interactions between frozen and non-frozen particles
appropriately, as if the frozen particle has infinite mass. A similar
functionality for normal (point) particles can be obtained using [fix
setforce](fix_setforce).

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

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix.

This fix computes a global 3-vector of forces, which can be accessed by
various [output commands](Howto_output). This is the total force on the
group of atoms before the forces on individual atoms are changed by the
fix. The vector values calculated by this fix are \"extensive\".

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

This fix is part of the GRANULAR package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

There can only be a single freeze fix defined. This is because other the
[granular pair styles](pair_gran) treat frozen particles differently and
need to be able to reference a single group to which this fix is
applied.

## Related commands

[atom_style sphere](atom_style), [fix setforce](fix_setforce)

## Default

none
