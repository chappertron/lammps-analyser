# fix spring/chunk command

## Syntax

    fix ID group-ID spring/chunk K chunkID comID

-   ID, group-ID are documented in [fix](fix) command
-   spring/chunk = style name of this fix command
-   K = spring constant for each chunk (force/distance units)
-   chunkID = ID of [compute chunk/atom](compute_chunk_atom) command
-   comID = ID of [compute com/chunk](compute_com_chunk) command

## Examples

``` LAMMPS
fix restrain all spring/chunk 100 chunkID comID
```

## Description

Apply a spring force to the center-of-mass (COM) of chunks of atoms as
defined by the [compute chunk/atom](compute_chunk_atom) command. Chunks
can be molecules or spatial bins or other groupings of atoms. This is a
way of tethering each chunk to its initial COM coordinates.

The *chunkID* is the ID of a compute chunk/atom command defined in the
input script. It is used to define the chunks. The *comID* is the ID of
a compute com/chunk command defined in the input script. It is used to
compute the COMs of each chunk.

At the beginning of the first [run](run) or [minimize](minimize) command
after this fix is defined, the initial COM of each chunk is calculated
and stored as R0m, where M is the chunk number. Thereafter, at every
timestep (or minimization iteration), the current COM of each chunk is
calculated as Rm. A restoring force of magnitude K (Rm - R0m) Mi / Mm is
applied to each atom in each chunk where *K* is the specified spring
constant, Mi is the mass of the atom, and Mm is the total mass of all
atoms in the chunk. Note that *K* thus represents the spring constant
for the total force on each chunk of atoms, not for a spring applied to
each atom.

## Restart, fix_modify, output, run start/stop, minimize info

This fix writes the locations of the initial per-chunk center of mass
coordinates to [binary restart files](restart). See the
[read_restart](read_restart) command for info on how to re-specify a fix
in an input script that reads a restart file, so that the fix continues
in an uninterrupted fashion. Since this fix depends on an instance of
[compute chunk/atom](compute_chunk_atom) it will check when reading the
restart if the chunk still exists and will define the same number of
chunks. The restart data is only applied when the number of chunks
matches. Otherwise the center of mass coordinates are recomputed.

The [fix_modify](fix_modify) *energy* option is supported by this fix to
add the energy stored in all the springs to the global potential energy
of the system as part of [thermodynamic output](thermo_style). The
default setting for this fix is [fix_modify energy no](fix_modify).

The [fix_modify](fix_modify) *respa* option is supported by this fix.
This allows to set at which level of the [r-RESPA](run_style) integrator
the fix is adding its forces. Default is the outermost level.

This fix computes a global scalar which can be accessed by various
[output commands](Howto_output). The scalar is the energy of all the
springs, i.e. 0.5 \* K \* r\^2 per-spring.

The scalar value calculated by this fix is \"extensive\".

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command.

The forces due to this fix are imposed during an energy minimization,
invoked by the [minimize](minimize) command.

:::: note
::: title
Note
:::

If you want the spring energies to be included in the total potential
energy of the system (the quantity being minimized), you MUST enable the
[fix_modify](fix_modify) *energy* option for this fix.
::::

## Restrictions

> none

## Related commands

[fix spring](fix_spring), [fix spring/self](fix_spring_self), [fix
spring/rg](fix_spring_rg)

## Default

none
