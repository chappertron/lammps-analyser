# mass command

## Syntax

    mass I value

-   I = atom type (see asterisk form below), or type label
-   value = mass

## Examples

``` LAMMPS
mass 1 1.0
mass * 62.5
mass 2* 62.5

labelmap atom 1 C
mass C 12.01
```

## Description

Set the mass for all atoms of one or more atom types. Per-type mass
values can also be set in the [read_data](read_data) data file using the
\"Masses\" keyword. See the [units](units) command for what mass units
to use.

The I index can be specified in one of several ways. An explicit numeric
value can be used, as in the first example above. Or I can be a type
label, which is an alphanumeric string defined by the
[labelmap](labelmap) command or in a section of a data file read by the
[read_data](read_data) command, and which converts internally to a
numeric type. Or a wild-card asterisk can be used to set the mass for
multiple atom types. This takes the form \"\*\" or \"\*n\" or \"n\*\" or
\"m\*n\", where m and n are numbers. If N = the number of atom types,
then an asterisk with no numeric values means all types from 1 to N. A
leading asterisk means all types from 1 to n (inclusive). A trailing
asterisk means all types from n to N (inclusive). A middle asterisk
means all types from m to n (inclusive).

A line in a [data file](read_data) that follows the \"Masses\" keyword
specifies mass using the same format as the arguments of the mass
command in an input script, except that no wild-card asterisk can be
used. For example, under the \"Masses\" section of a data file, the line
that corresponds to the first example above would be listed as

    1 1.0

Note that the mass command can only be used if the [atom
style](atom_style) requires per-type atom mass to be set. Currently, all
but the *sphere* and *ellipsoid* and *peri* styles do. They require mass
to be set for individual particles, not types. Per-atom masses are
defined in the data file read by the [read_data](read_data) command, or
set to default values by the [create_atoms](create_atoms) command.
Per-atom masses can also be set to new values by the [set mass](set) or
[set density](set) commands.

Also note that [pair_style eam](pair_eam) and [pair_style bop](pair_bop)
commands define the masses of atom types in their respective potential
files, in which case the mass command is normally not used.

If you define a [hybrid atom style](atom_style) which includes one (or
more) sub-styles which require per-type mass and one (or more)
sub-styles which require per-atom mass, then you must define both.
However, in this case the per-type mass will be ignored; only the
per-atom mass will be used by LAMMPS.

## Restrictions

This command must come after the simulation box is defined by a
[read_data](read_data), [read_restart](read_restart), or
[create_box](create_box) command.

All masses must be defined before a simulation is run. They must also
all be defined before a [velocity](velocity) or [fix shake](fix_shake)
command is used.

The mass assigned to any type or atom must be \> 0.0.

## Related commands

none

## Default

none
