# Granular models

Granular system are composed of spherical particles with a diameter, as
opposed to point particles. This means they have an angular velocity and
torque can be imparted to them to cause them to rotate.

To run a simulation of a granular model, you will want to use the
following commands:

-   [atom_style sphere](atom_style)
-   [fix nve/sphere](fix_nve_sphere)
-   [fix gravity](fix_gravity)

This compute

-   [compute erotate/sphere](compute_erotate_sphere)

calculates rotational kinetic energy which can be [output with
thermodynamic info](Howto_output). The compute

-   [compute fabric](compute_fabric)

calculates various versions of the fabric tensor for granular and
non-granular pair styles.

Use one of these 4 pair potentials, which compute forces and torques
between interacting pairs of particles:

-   [pair_style gran/history](pair_gran)
-   [pair_style gran/no_history](pair_gran)
-   [pair_style gran/hertzian](pair_gran)
-   [pair_style granular](pair_granular)

These commands implement fix options specific to granular systems:

-   [fix freeze](fix_freeze)
-   [fix pour](fix_pour)
-   [fix viscous](fix_viscous)
-   [fix wall/gran](fix_wall_gran)
-   [fix wall/gran/region](fix_wall_gran_region)

The fix style *freeze* zeroes both the force and torque of frozen atoms,
and should be used for granular system instead of the fix style
*setforce*.

To model heat conduction, one must add the temperature and heatflow atom
variables with: \* [fix property/atom](fix_property_atom) a temperature
integration fix \* [fix heat/flow](fix_heat_flow) and a heat conduction
option defined in both \* [pair_style granular](pair_granular) \* [fix
wall/gran](fix_wall_gran)

For computational efficiency, you can eliminate needless pairwise
computations between frozen atoms by using this command:

-   [neigh_modify](neigh_modify) exclude

:::: note
::: title
Note
:::

By default, for 2d systems, granular particles are still modeled as 3d
spheres, not 2d discs (circles), meaning their moment of inertia will be
the same as in 3d. If you wish to model granular particles in 2d as 2d
discs, see the note on this topic on the [Howto 2d](Howto_2d) doc page,
where 2d simulations are discussed.
::::

To add custom granular contact models, see the [modifying granular
sub-models page](Modify_gran_sub_mod).
