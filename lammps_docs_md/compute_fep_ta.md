# compute fep/ta command

## Syntax

``` LAMMPS
compute ID group-ID fep/ta temp plane scale_factor keyword value ...
```

-   ID, group-ID are documented in the [compute](compute) command

-   fep/ta = name of this compute command

-   temp = external temperature (as specified for constant-temperature
    run)

-   plane = *xy* or *xz* or *yz*

-   scale_factor = multiplicative factor for change in plane area

-   zero or more keyword/value pairs may be appended

-   keyword = *tail*

        *tail* value = *no* or *yes*
          *no* = ignore tail correction to pair energies (usually small in fep)
          *yes* = include tail correction to pair energies

## Examples

``` LAMMPS
compute 1 all fep/ta 298 xy 1.0005
```

## Description

::: versionadded
4May2022
:::

Define a computation that calculates the change in the free energy due
to a test-area (TA) perturbation [(Gloor)](Gloor). The test-area
approach can be used to determine the interfacial tension of the system
in a single simulation:

$$\gamma = \lim_{\Delta \mathcal{A} \to 0} \left( \frac{\Delta A_{0 \to 1 }}{\Delta \mathcal{A}}\right)_{N,V,T}
= - \frac{k_B T}{\Delta \mathcal{A}} \ln \left\langle \exp\left(\frac{-(U_1 - U_0)}{k_B T}\right) \right\rangle_0$$

During the perturbation, both axes of *plane* are scaled by multiplying
$\sqrt{\mathrm{scale\_factor}}$, while the other axis divided by
$\mathrm{scale\_factor}$ such that the overall volume of the system is
maintained.

The *tail* keyword controls the calculation of the tail correction to
\"van der Waals\" pair energies beyond the cutoff, if this has been
activated via the [pair_modify](pair_modify) command. If the
perturbation is small, the tail contribution to the energy difference
between the reference and perturbed systems should be negligible.

------------------------------------------------------------------------

## Output info

This compute calculates a global vector of length 3 which contains the
energy difference $(U_1-U_0)$ as c_ID\[1\], the Boltzmann factor
$\exp\bigl(-(U_1-U_0)/k_B T\bigr)$, as c_ID\[2\] and the change in the
*plane* area $\Delta \mathcal{A}$ as c_ID\[3\]. $U_1$ is the potential
energy of the perturbed state and $U_0$ is the potential energy of the
reference state. The energies include kspace terms if these are used in
the simulation.

These output results can be used by any command that uses a global
scalar or vector from a compute as input. See the [Howto
output](Howto_output) page for an overview of LAMMPS output options. For
example, the computed values can be averaged using [fix
ave/time](fix_ave_time).

## Restrictions

Constraints, like fix shake, may lead to incorrect values for energy
difference.

This compute is distributed as the FEP package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

## Related commands

[compute fep](compute_fep)

## Default

The option defaults are *tail* = *no*.

------------------------------------------------------------------------

::: {#Gloor}
**(Gloor)** Gloor, J Chem Phys, 123, 134703 (2005)
:::
