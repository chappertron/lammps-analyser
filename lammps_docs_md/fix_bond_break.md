# fix bond/break command

## Syntax

``` LAMMPS
fix ID group-ID bond/break Nevery bondtype Rmax keyword values ...
```

-   ID, group-ID are documented in [fix](fix) command

-   bond/break = style name of this fix command

-   Nevery = attempt bond breaking every this many steps

-   bondtype = type of bonds to break

-   Rmax = bond longer than Rmax can break (distance units)

-   zero or more keyword/value pairs may be appended

-   keyword = *prob*

        *prob* values = fraction seed
          fraction = break a bond with this probability if otherwise eligible
          seed = random number seed (positive integer)

## Examples

``` LAMMPS
fix 5 all bond/break 10 2 1.2
fix 5 polymer bond/break 1 1 2.0 prob 0.5 49829
```

## Description

Break bonds between pairs of atoms as a simulation runs according to
specified criteria. This can be used to model the dissolution of a
polymer network due to stretching of the simulation box or other
deformations. In this context, a bond means an interaction between a
pair of atoms computed by the [bond_style](bond_style) command. Once the
bond is broken it will be permanently deleted, as will all angle,
dihedral, and improper interactions that bond is part of.

This is different than a [pair-wise](pair_style) bond-order potential
such as Tersoff or AIREBO which infers bonds and many-body interactions
based on the current geometry of a small cluster of atoms and
effectively creates and destroys bonds and higher-order many-body
interactions from timestep to timestep as atoms move.

A check for possible bond breakage is performed every *Nevery*
timesteps. If two bonded atoms $i$ and $j$ are farther than the distance
*Rmax* from each other, the bond is of type *bondtype*, and both $i$ and
$j$ are in the specified fix group, then the bond between $i$ and $j$ is
labeled as a \"possible\" bond to break.

If several bonds involving an atom are stretched, it may have multiple
possible bonds to break. Every atom checks its list of possible bonds to
break and labels the longest such bond as its \"sole\" bond to break.
After this is done, if atom $i$ is bonded to atom $j$ in its sole bond,
and atom $j$ is bonded to atom $j$ in its sole bond, then the bond
between $i$ and $j$ is \"eligible\" to be broken.

Note that these rules mean an atom will only be part of at most one
broken bond on a given time step. It also means that if atom $i$ chooses
atom $j$ as its sole partner, but atom $j$ chooses atom $k$ as its sole
partner (because $R_{jk} > R_{ij}$), then this means atom $i$ will not
be part of a broken bond on this time step, even if it has other
possible bond partners.

The *prob* keyword can effect whether an eligible bond is actually
broken. The *fraction* setting must be a value between 0.0 and 1.0. A
uniform random number between 0.0 and 1.0 is generated and the eligible
bond is only broken if the random number is less than *fraction*.

When a bond is broken, data structures within LAMMPS that store bond
topologies are updated to reflect the breakage. Likewise, if the bond is
part of a 3-body (angle) or 4-body (dihedral, improper) interaction,
that interaction is removed as well. These changes typically affect
pair-wise interactions between atoms that used to be part of bonds,
angles, etc.

:::: note
::: title
Note
:::

One data structure that is not updated when a bond breaks are the
molecule IDs stored by each atom. Even though one molecule becomes two
molecules due to the broken bond, all atoms in both new molecules retain
their original molecule IDs.
::::

Computationally, each time step this fix is invoked, it loops over all
the bonds in the system and computes distances between pairs of bonded
atoms. It also communicates between neighboring processors to coordinate
which bonds are broken. Moreover, if any bonds are broken, neighbor
lists must be immediately updated on the same time step. This is to
ensure that any pair-wise interactions that should be turned \"on\" due
to a bond breaking, because they are no longer excluded by the presence
of the bond and the settings of the [special_bonds](special_bonds)
command, will be immediately recognized. All of these operations
increase the cost of a time step. Thus, you should be cautious about
invoking this fix too frequently.

You can dump out snapshots of the current bond topology via the [dump
local](dump) command.

:::: note
::: title
Note
:::

Breaking a bond typically alters the energy of a system. You should be
careful not to choose bond breaking criteria that induce a dramatic
change in energy. For example, if you define a very stiff harmonic bond
and break it when two atoms are separated by a distance far from the
equilibrium bond length, then the two atoms will be dramatically
released when the bond is broken. More generally, you may need to
thermostat your system to compensate for energy changes resulting from
broken bonds (as well as angles, dihedrals, and impropers).
::::

See the [Howto](Howto_broken_bonds) page on broken bonds for more
information on related features in LAMMPS.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix.

This fix computes two statistics, which it stores in a global vector of
length 2. This vector can be accessed by various [output
commands](Howto_output). The vector values calculated by this fix are
\"intensive\".

The two quantities in the global vector are

> (1) number of bonds broken on the most recent breakage time step
> (2) cumulative number of bonds broken

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

This fix is part of the MC package. It is only enabled if LAMMPS was
built with that package. See the [Build package](Build_package) doc page
for more info.

## Related commands

[fix bond/create](fix_bond_create), [fix bond/react](fix_bond_react),
[fix bond/swap](fix_bond_swap), [dump local](dump),
[special_bonds](special_bonds)

## Default

The option defaults are prob = 1.0.
