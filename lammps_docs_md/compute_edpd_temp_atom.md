# compute edpd/temp/atom command

## Syntax

``` LAMMPS
compute ID group-ID edpd/temp/atom
```

-   ID, group-ID are documented in [compute](compute) command
-   edpd/temp/atom = style name of this compute command

## Examples

``` LAMMPS
compute 1 all edpd/temp/atom
```

## Description

Define a computation that calculates the per-atom temperature for each
eDPD particle in a group.

The temperature is a local temperature derived from the internal energy
of each eDPD particle based on the local equilibrium hypothesis. For
more details please see [(Espanol1997)](Espanol1997) and
[(Li2014)](Li2014a).

## Output info

This compute calculates a per-atom vector, which can be accessed by any
command that uses per-atom values from a compute as input. See the
[Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-atom vector values will be in temperature [units](units).

## Restrictions

This compute is part of the DPD-MESO package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[pair_style edpd](pair_mesodpd)

## Default

none

------------------------------------------------------------------------

::: {#Espanol1997}
**(Espanol1997)** Espanol, Europhys Lett, 40(6): 631-636 (1997). DOI:
10.1209/epl/i1997-00515-8
:::

::: {#Li2014a}
**(Li2014)** Li, Tang, Lei, Caswell, Karniadakis, J Comput Phys, 265:
113-127 (2014). DOI: 10.1016/j.jcp.2014.02.003.
:::
