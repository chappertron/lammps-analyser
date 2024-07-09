# bond_style harmonic/restrain command

## Syntax

``` LAMMPS
bond_style harmonic/restrain
```

## Examples

``` LAMMPS
bond_style harmonic
bond_coeff 5 80.0
```

## Description

::: versionadded
28Mar2023
:::

The *harmonic/restrain* bond style uses the potential

$$E = K (r - r_{t=0})^2$$

where $r_{t=0}$ is the distance between the bonded atoms at the
beginning of the first [run](run) or [minimize](minimize) command after
the bond style has been defined (*t=0*). Note that the usual 1/2 factor
is included in $K$. This will effectively restrain bonds to their
initial length, whatever that is. This is where this bond style differs
from [bond style harmonic](bond_harmonic) where the bond length is set
through the per bond type coefficients.

The following coefficient must be defined for each bond type via the
[bond_coeff](bond_coeff) command as in the example above, or in the data
file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands

-   $K$ (energy/distance\^2)

This bond style differs from other options to add harmonic restraints
like [fix restrain](fix_restrain) or [pair style list](pair_list) or
[fix colvars](fix_colvars) in that it requires a bond topology, and thus
the defined bonds will trigger exclusion of special neighbors from the
neighbor list according to the [special_bonds](special_bonds) settings.

## Restart info

This bond style supports the [write_restart](write_restart) and
[read_restart](read_restart) commands. The state of the initial bond
lengths is stored with restart files and read back.

## Restrictions

This bond style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package. See the [Build package](Build_package) page for
more info.

This bond style maintains internal data to determine the original bond
lengths $r_{t=0}$. This information will be written to [binary restart
files](write_restart) but **not** to [data files](write_data). Thus,
continuing a simulation is *only* possible with
[read_restart](read_restart). When using the [read_data
command](read_data), the reference bond lengths $r_{t=0}$ will be
re-initialized from the current geometry.

This bond style cannot be used with [fix shake or fix
rattle](fix_shake), with [fix filter/corotate](fix_filter_corotate), or
any [tip4p pair style](pair_lj_cut_tip4p) since there is no specific
equilibrium distance for a given bond type.

## Related commands

[bond_coeff](bond_coeff), [bond_harmonic](bond_harmonic), [fix
restrain](fix_restrain), [pair style list](pair_list)

## Default

none
