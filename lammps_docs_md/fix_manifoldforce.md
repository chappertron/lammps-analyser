# fix manifoldforce command

## Syntax

    fix ID group-ID manifoldforce manifold manifold-args ...

-   ID, group-ID are documented in [fix](fix) command
-   manifold = name of the manifold
-   manifold-args = parameters for the manifold

## Examples

``` LAMMPS
fix constrain all manifoldforce sphere 5.0
```

## Description

This fix subtracts each time step from the force the component along the
normal of the specified [manifold](Howto_manifold). This can be used in
combination with [minimize](minimize) to remove overlap between
particles while keeping them (roughly) constrained to the given
manifold, e.g. to set up a run with [fix
nve/manifold/rattle](fix_nve_manifold_rattle). I have found that only
*hftn* and *quickmin* with a very small time step perform adequately
though.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command. This fix is invoked during [energy
minimization](minimize).

------------------------------------------------------------------------

## Restrictions

This fix is part of the MANIFOLD package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

Only use this with *min_style hftn* or *min_style quickmin*. If not, the
constraints will not be satisfied very well at all. A warning is
generated if the *min_style* is incompatible but no error.

------------------------------------------------------------------------

## Related commands

[fix nve/manifold/rattle](fix_nve_manifold_rattle), [fix
nvt/manifold/rattle](fix_nvt_manifold_rattle)
