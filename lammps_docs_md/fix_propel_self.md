# fix propel/self command

## Syntax

    fix ID group-ID propel/self mode magnitude keyword values

-   ID, group-ID are documented in [fix](fix) command

-   propel/self = style name of this fix command

-   mode = *dipole* or *velocity* or *quat*

-   magnitude = magnitude of self-propulsion force

-   zero or one keyword/value pairs may be appended

-   keyword = *qvector*

        *qvector* value = direction of force in ellipsoid frame
         *sx*, *sy*, *sz* = components of *qvector*

## Examples

``` LAMMPS
fix active all propel/self dipole 40.0
fix active all propel/self velocity 10.0
fix active all propel/self quat 15.7 qvector 1.0 0.0 0.0
```

## Description

Add a force to each atom in the group due to a self-propulsion force.
The force is given by

$$F_i = f_P e_i$$

where *i* is the particle the force is being applied to, $f_P$ is the
magnitude of the force, and $e_i$ is the vector direction of the force.
The specification of $e_i$ is based on which of the three keywords
(*dipole* or *velocity* or *quat*) one selects.

For mode *dipole*, $e_i$ is just equal to the dipole vectors of the
atoms in the group. Therefore, if the dipoles are not unit vectors, the
$e_i$ will not be unit vectors.

:::: note
::: title
Note
:::

If another command changes the magnitude of the dipole, this force will
change accordingly (since $|e_i|$ will change, which is physically
equivalent to re-scaling $f_P$ while keeping $|e_i|$ constant), and no
warning will be provided by LAMMPS. This is almost never what you want,
so ensure you are not changing dipole magnitudes with another LAMMPS fix
or pair style. Furthermore, self-propulsion forces (almost) always set
$e_i$ to be a unit vector for all times, so it\'s best to set all the
dipole magnitudes to 1.0 unless you have a good reason not to (see the
[set](set) command on how to do this).
::::

For mode *velocity*, $e_i$ points in the direction of the current
velocity (a unit-vector). This can be interpreted as a
velocity-dependent friction, as proposed by e.g. [(Erdmann)](Erdmann1).

For mode *quat*, $e_i$ points in the direction of a unit vector,
oriented in the coordinate frame of the ellipsoidal particles, which
defaults to point along the x-direction. This default behavior can be
changed by via the *quatvec* keyword.

The optional *quatvec* keyword specifies the direction of
self-propulsion via a unit vector (sx,sy,sz). The arguments *sx*, *sy*,
and *sz*, are defined within the coordinate frame of the atom\'s
ellipsoid. For instance, for an ellipsoid with long axis along its
x-direction, if one wanted the self-propulsion force to also be along
this axis, set *sx* equal to 1 and *sy*, *sz* both equal to zero. This
keyword may only be specified for mode *quat*.

:::: note
::: title
Note
:::

In using keyword *quatvec*, the three arguments *sx*, *sy*, and *sz*
will be automatically normalized to components of a unit vector
internally to avoid users having to explicitly do so themselves.
Therefore, in mode *quat*, the vectors $e_i$ will always be of unit
length.
::::

Along with adding a force contribution, this fix can also contribute to
the virial (pressure) of the system, defined as
$f_P \sum_i <e_i . r_i>/(d V)$, where $r_i$ is the *unwrapped*
coordinate of particle i in the case of periodic boundary conditions.
See [(Winkler)](Winkler1) for a discussion of this active pressure
contribution.

For modes *dipole* and *quat*, this fix is by default included in
pressure computations.

For mode *velocity*, this fix is by default not included in pressure
computations.

:::: note
::: title
Note
:::

In contrast to equilibrium systems, pressure of active systems in
general depends on the geometry of the container. The active pressure
contribution as calculated in this fix is only valid for certain
boundary conditions (spherical walls, rectangular walls, or periodic
boundary conditions). For other geometries, the pressure must be
measured via explicit calculation of the force per unit area on a wall,
and so one must not calculate it using this fix. (Use
[fix_modify](fix_modify) as described below to turn off the virial
contribution of this fix). Again, see [(Winkler)](Winkler1) for
discussion of why this is the case.

Furthermore, when dealing with active systems, the temperature is no
longer well defined. Therefore, one should ensure that the *virial* flag
is used in the [compute pressure](compute_pressure) command (turning off
temperature contributions).
::::

------------------------------------------------------------------------

## Restart, fix_modify, output, run start/stop, minimize info

No information about this fix is written to [binary restart
files](restart).

The [fix_modify](fix_modify) *virial* option is supported by this fix to
add the contribution due to the added forces on atoms to the system\'s
virial as part of [thermodynamic output](thermo_style). The default is
*virial yes* for keywords *dipole* and *quat*. The default is *virial
no* for keyword *velocity*.

No parameter of this fix can be used with the *start/stop* keywords of
the [run](run) command.

## Restrictions

With keyword *dipole*, this fix only works when the DIPOLE package is
enabled. See the [Build package](Build_package) page for more info.

This fix is part of the BROWNIAN package. It is only enabled if LAMMPS
was built with that package. See the [Build package](Build_package) doc
page for more info.

## Related commands

[fix efield](fix_efield) , [fix setforce](fix_setforce), [fix
addforce](fix_addforce)

## Default

none

------------------------------------------------------------------------

::: {#Erdmann1}
**(Erdmann)** U. Erdmann , W. Ebeling, L. Schimansky-Geier, and F.
Schweitzer, Eur. Phys. J. B 15, 105-113, 2000.
:::

::: {#Winkler1}
**(Winkler)** Winkler, Wysocki, and Gompper, Soft Matter, 11, 6680
(2015).
:::
