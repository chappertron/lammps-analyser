# fix reaxff/bonds command

Accelerator Variants: *reaxff/bonds/kk*

## Syntax

    fix ID group-ID reaxff/bonds Nevery filename

-   ID, group-ID are documented in [fix](fix) command
-   reax/bonds = style name of this fix command
-   Nevery = output interval in timesteps
-   filename = name of output file

## Examples

``` LAMMPS
fix 1 all reaxff/bonds 100 bonds.reaxff
```

## Description

Write out the bond information computed by the ReaxFF potential
specified by [pair_style reaxff](pair_reaxff) in the exact same format
as the original stand-alone ReaxFF code of Adri van Duin. The bond
information is written to *filename* on timesteps that are multiples of
*Nevery*, including timestep 0. For time-averaged chemical species
analysis, please see the [fix reaxff/species](fix_reaxff_species)
command.

The specified group-ID is ignored by this fix.

The format of the output file should be reasonably self-explanatory. The
meaning of the column header abbreviations is as follows:

-   id = atom id
-   type = atom type
-   nb = number of bonds
-   id_1 = atom id of first bond
-   id_nb = atom id of Nth bond
-   mol = molecule id
-   bo_1 = bond order of first bond
-   bo_nb = bond order of Nth bond
-   abo = atom bond order (sum of all bonds)
-   nlp = number of lone pairs
-   q = atomic charge

If the filename ends with \".gz\", the output file is written in gzipped
format. A gzipped dump file will be about 3x smaller than the text
version, but will also take longer to write.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is not invoked during [energy
minimization](minimize).

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

The fix reaxff/bonds command requires that the [pair_style
reaxff](pair_reaxff) is invoked. This fix is part of the REAXFF package.
It is only enabled if LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

To write gzipped bond files, you must compile LAMMPS with the
-DLAMMPS_GZIP option.

## Related commands

[pair_style reaxff](pair_reaxff), [fix
reaxff/species](fix_reaxff_species)

## Default

none
