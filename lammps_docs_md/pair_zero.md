# pair_style zero command

## Syntax

``` LAMMPS
pair_style zero cutoff [nocoeff] [full]
```

-   zero = style name of this pair style
-   cutoff = global cutoff (distance units)
-   nocoeff = ignore all pair_coeff parameters (optional)
-   full = build full neighbor list (optional)

## Examples

``` LAMMPS
pair_style zero 10.0
pair_style zero 5.0 nocoeff
pair_coeff * *
pair_coeff 1 2*4 3.0
```

## Description

Define a global or per-type cutoff length for the purpose of building a
neighbor list and acquiring ghost atoms, but do not compute any pairwise
forces or energies.

This can be useful for fixes or computes which require a neighbor list
to enumerate pairs of atoms within some cutoff distance, but when
pairwise forces are not otherwise needed. Examples are the [fix
bond/create](fix_bond_create), [compute rdf](compute_rdf), [compute
voronoi/atom](compute_voronoi_atom) commands.

Note that the [comm_modify cutoff](comm_modify) command can be used to
ensure communication of ghost atoms even when a pair style is not
defined, but it will not trigger neighbor list generation.

The optional *nocoeff* flag allows to read data files with a PairCoeff
section for any pair style. Similarly, any pair_coeff commands will only
be checked for the atom type numbers and the rest ignored. In this case,
only the global cutoff will be used.

::: versionadded
3Nov2022
:::

The optional *full* flag builds a full neighbor list instead of the
default half neighbor list.

The following coefficients must be defined for each pair of atoms types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands, or by mixing as described below:

-   cutoff (distance units)

This coefficient is optional. If not specified, the global cutoff
specified in the pair_style command is used. If the pair_style has been
specified with the optional *nocoeff* flag, then a cutoff pair
coefficient is ignored.

------------------------------------------------------------------------

## Mixing, shift, table, tail correction, restart, rRESPA info

The cutoff distance for this pair style can be mixed. The default mix
value is *geometric*. See the \"pair_modify\" command for details.

This pair style does not support the [pair_modify](pair_modify) shift,
table, and tail options.

This pair style writes its information to [binary restart
files](restart), so pair_style and pair_coeff commands do not need to be
specified in an input script that reads a restart file.

This pair style supports the use of the *inner*, *middle*, and *outer*
keywords of the [run_style respa](run_style) command.

------------------------------------------------------------------------

## Restrictions

none

## Related commands

[pair_style none](pair_none)

## Default

none
