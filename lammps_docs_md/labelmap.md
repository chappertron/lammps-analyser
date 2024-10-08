# labelmap command

## Syntax

``` LAMMPS
labelmap option args
```

-   *option* = *atom* or *bond* or *angle* or *dihedral* or *improper*
    or *clear* or *write*

        *clear* = no args
        *write* arg = filename
        *atom* or *bond* or *angle* or *dihedral* or *improper*
          args = list of one or more numeric-type/type-label pairs

## Examples

``` LAMMPS
labelmap atom 3 carbon 4 'c3"' 5 "c1'" 6 "c#"
labelmap atom $(label2type(atom,carbon)) C  # change type label from 'carbon' to 'C'
labelmap clear
labelmap write mymap.include
labelmap bond 1 carbonyl 2 nitrile 3 """ c1'-c2" """
```

## Description

::: versionadded
15Sep2022
:::

Define alphanumeric type labels to associate with one or more numeric
atom, bond, angle, dihedral or improper types. A collection of type
labels for all atom types, bond types, etc. is stored as a label map.

The label map can also be defined by the [read_data](read_data) command
when it reads these sections in a data file: Atom Type Labels, Bond Type
Labels, etc. See the [Howto type labels](Howto_type_labels) doc page for
a general discussion of how type labels can be used.

Valid type labels can contain any alphanumeric character, but must not
start with a number, a \'#\', or a \'\*\' character. They can contain
other standard ASCII characters such as angular or square brackets
\'\<\' and \'\>\' or \'\[\' and \'\]\', parenthesis \'(\' and \')\',
dash \'-\', underscore \'\_\', plus \'+\' and equals \'=\' signs and
more. They must not contain blanks or any other whitespace. Note that
type labels must be put in single or double quotation marks if they
contain the \'#\' character or if they contain a double (\") or single
quotation mark (\'). If the label contains both a single and a double
quotation mark, then triple quotation (\"\"\") must be used. When
enclosing a type label with quotation marks, the LAMMPS input parser may
require adding leading or trailing blanks around the type label so it
can identify the enclosing quotation marks. Those blanks will be removed
when defining the label.

A *labelmap* command can only modify the label map for one type-kind
(atom types, bond types, etc). Any number of numeric-type/type-label
pairs may follow. If a type label already exists for the same numeric
type, it will be overwritten. Type labels must be unique; assigning the
same type label to multiple numeric types within the same type-kind is
not allowed. When reading and writing data files, it is required that
there is a label defined for *every* numeric type within a given
type-kind in order to write out the type label section for that
type-kind.

The *clear* option resets the label map and thus discards all previous
settings.

The *write* option takes a filename as argument and writes the current
label mappings to a file as a sequence of *labelmap* commands, so the
file can be copied into a new LAMMPS input file or read in using the
[include](include) command.

------------------------------------------------------------------------

## Restrictions

This command must come after the simulation box is defined by a
[read_data](read_data), [read_restart](read_restart), or
[create_box](create_box) command.

Label maps are currently not supported when using the KOKKOS package.

## Related commands

[read_data](read_data), [write_data](write_data), [molecule](molecule),
[fix bond/react](fix_bond_react)

## Default

none
