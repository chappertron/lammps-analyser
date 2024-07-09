# write_coeff command

## Syntax

``` LAMMPS
write_coeff file
```

-   file = name of data file to write out

## Examples

``` LAMMPS
write_coeff polymer.coeff
```

## Description

Write a text format file with the currently defined force field
coefficients in a way, that it can be read by LAMMPS with the
[include](include) command. In combination with the nocoeff option of
[write_data](write_data) this can be used to move the Coeffs sections
from a data file into a separate file.

:::: note
::: title
Note
:::

The write_coeff command is not yet fully implemented as some pair styles
do not output their coefficient information. This means you will need to
add/copy this information manually.
::::

------------------------------------------------------------------------

## Restrictions

none

## Related commands

[read_data](read_data), [write_restart](write_restart),
[write_data](write_data)
