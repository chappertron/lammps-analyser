# fix store/force command

## Syntax

    fix ID group-ID store/force

-   ID, group-ID are documented in [fix](fix) command
-   store/force = style name of this fix command

## Examples

``` LAMMPS
fix 1 all store/force
```

## Description

Store the forces on atoms in the group at the point during each timestep
when the fix is invoked, as described below. This is useful for storing
forces before constraints or other boundary conditions are computed
which modify the forces, so that unmodified forces can be [written to a
dump file](dump) or accessed by other [output commands](Howto_output)
that use per-atom quantities.

This fix is invoked at the point in the velocity-Verlet timestepping
immediately after [pair](pair_style), [bond](bond_style),
[angle](angle_style), [dihedral](dihedral_style),
[improper](improper_style), and [long-range](kspace_style) forces have
been calculated. It is the point in the timestep when various fixes that
compute constraint forces are calculated and potentially modify the
force on each atom. Examples of such fixes are [fix shake](fix_shake),
[fix wall](fix_wall), and [fix indent](fix_indent).

:::: note
::: title
Note
:::

The order in which various fixes are applied which operate at the same
point during the timestep, is the same as the order they are specified
in the input script. Thus normally, if you want to store per-atom forces
due to force field interactions, before constraints are applied, you
should list this fix first within that set of fixes, i.e. before other
fixes that apply constraints. However, if you wish to include certain
constraints (e.g. fix shake) in the stored force, then it could be
specified after some fixes and before others.
::::

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix.

This fix produces a per-atom array which can be accessed by various
[output commands](Howto_output). The number of columns for each atom is
3, and the columns store the x,y,z forces on each atom. The per-atom
values be accessed on any timestep.

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command. This fix is not invoked during [energy
minimization](minimize).

## Restrictions

> none

## Related commands

[fix store_state](fix_store_state)

## Default

none
