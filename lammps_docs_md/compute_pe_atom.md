# compute pe/atom command

## Syntax

``` LAMMPS
compute ID group-ID pe/atom keyword ...
```

-   ID, group-ID are documented in [compute](compute) command
-   pe/atom = style name of this compute command
-   zero or more keywords may be appended
-   keyword = *pair* or *bond* or *angle* or *dihedral* or *improper* or
    *kspace* or *fix*

## Examples

``` LAMMPS
compute 1 all pe/atom
compute 1 all pe/atom pair
compute 1 all pe/atom pair bond
```

## Description

Define a computation that computes the per-atom potential energy for
each atom in a group. See the [compute pe](compute_pe) command if you
want the potential energy of the entire system.

The per-atom energy is calculated by the various pair, bond, etc
potentials defined for the simulation. If no extra keywords are listed,
then the potential energy is the sum of pair, bond, angle, dihedral,
improper, $k$-space (long-range), and fix energy (i.e., it is as though
all the keywords were listed). If any extra keywords are listed, then
only those components are summed to compute the potential energy.

Note that the energy of each atom is due to its interaction with all
other atoms in the simulation, not just with other atoms in the group.

For an energy contribution produced by a small set of atoms (e.g., 4
atoms in a dihedral or 3 atoms in a Tersoff 3-body interaction), that
energy is assigned in equal portions to each atom in the set (e.g., 1/4
of the dihedral energy to each of the four atoms).

The [dihedral_style charmm](dihedral_charmm) style calculates pairwise
interactions between 1\--4 atoms. The energy contribution of these terms
is included in the pair energy, not the dihedral energy.

The KSpace contribution is calculated using the method in
[(Heyes)](Heyes1) for the Ewald method and a related method for PPPM, as
specified by the [kspace_style pppm](kspace_style) command. For PPPM,
the calculation requires 1 extra FFT each timestep that per-atom energy
is calculated. This [document](PDF/kspace.pdf)\_ describes how the
long-range per-atom energy calculation is performed.

Various fixes can contribute to the per-atom potential energy of the
system if the *fix* contribution is included. See the doc pages for
[individual fixes](fix) for details of which ones compute a per-atom
potential energy.

:::: note
::: title
Note
:::

The [fix_modify energy yes](fix_modify) command must also be specified
if a fix is to contribute per-atom potential energy to this command.
::::

As an example of per-atom potential energy compared to total potential
energy, these lines in an input script should yield the same result in
the last 2 columns of thermo output:

``` LAMMPS
compute        peratom all pe/atom
compute        pe all reduce sum c_peratom
thermo_style   custom step temp etotal press pe c_pe
```

:::: note
::: title
Note
:::

The per-atom energy does not include any Lennard-Jones tail corrections
to the energy added by the [pair_modify tail yes](pair_modify) command,
since those are contributions to the global system energy.
::::

## Output info

This compute calculates a per-atom vector, which can be accessed by any
command that uses per-atom values from a compute as input. See the
[Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-atom vector values will be in energy [units](units).

## Restrictions

## Related commands

[compute pe](compute_pe), [compute stress/atom](compute_stress_atom)

## Default

none

------------------------------------------------------------------------

::: {#Heyes1}
**(Heyes)** Heyes, Phys Rev B 49, 755 (1994),
:::
