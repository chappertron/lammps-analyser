# fix viscous/sphere command

## Syntax

    fix ID group-ID viscous/sphere gamma keyword values ...

-   ID, group-ID are documented in [fix](fix) command

-   viscous/sphere = style name of this fix command

-   gamma = damping coefficient (torque/angular velocity units)

-   zero or more keyword/value pairs may be appended

        keyword = *scale*
          *scale* values = *type ratio* or *v_name*
            type = atom type (1-N)
            ratio = factor to scale the damping coefficients by
            v_name = reference to atom style variable *name*

## Examples

``` LAMMPS
fix 1 flow viscous/sphere 0.1
fix 1 damp viscous/sphere 0.5 scale 3 2.5
fix 1 damp viscous/sphere 0.5 scale v_radscale
```

## Description

Add a viscous damping torque to finite-size spherical particles in the
group that is proportional to the angular velocity of the atom. In
granular simulations this can be useful for draining the rotational
kinetic energy from the system in a controlled fashion. If used without
additional thermostatting (to add kinetic energy to the system), it has
the effect of slowly (or rapidly) freezing the system; hence it can also
be used as a simple energy minimization technique.

The damping torque $T_i$ is given by $T_i = - \gamma \omega_i$. The
larger the coefficient, the faster the rotational kinetic energy is
reduced.

If the optional keyword *scale* is used, $\gamma$ can be scaled up or
down by the specified factor for atoms. This factor can be set for
different atom types and thus the *scale* keyword used multiple times
followed by the atom type and the associated scale factor. Alternately
the scaling factor can be computed for each atom (e.g. based on its
radius) by using an [atom-style variable](variable).

:::: note
::: title
Note
:::

You should specify gamma in torque/angular velocity units. This is not
the same as mass/time units, at least for some of the LAMMPS
[units](units) options like \"real\" or \"metal\" that are not
self-consistent.
::::

In the current implementation, rather than have the user specify a
viscosity, $\gamma$ is specified directly in torque/angular velocity
units. If needed, $\gamma$ can be adjusted for atoms of different sizes
(i.e. $\sigma$) by using the *scale* keyword.

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart). None of the [fix_modify](fix_modify) options are
relevant to this fix. No global or per-atom quantities are stored by
this fix for access by various [output commands](Howto_output). No
parameter of this fix can be used with the *start/stop* keywords of the
[run](run) command.

The [fix_modify](fix_modify) *respa* option is supported by this fix.
This allows to set at which level of the [r-RESPA](run_style) integrator
the fix is modifying torques. Default is the outermost level.

The torques due to this fix are imposed during an energy minimization,
invoked by the [minimize](minimize) command. This fix should only be
used with damped dynamics minimizers that allow for non-conservative
forces. See the [min_style](min_style) command for details.

## Restrictions

This fix is part of the EXTRA-FIX package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) page
for more info.

This fix requires that atoms store torque and angular velocity (omega)
and a radius as defined by the [atom_style sphere](atom_style) command.

All particles in the group must be finite-size spheres. They cannot be
point particles.

## Related commands

[fix viscous](fix_viscous), [fix damping/cundall](fix_damping_cundall)

## Default

none
