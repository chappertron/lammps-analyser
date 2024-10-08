# create_box command

## Syntax

``` LAMMPS
create_box N region-ID keyword value ...
```

-   N = \# of atom types to use in this simulation

-   region-ID = ID of region to use as simulation domain

-   zero or more keyword/value pairs may be appended

-   keyword = *bond/types* or *angle/types* or *dihedral/types* or
    *improper/types* or *extra/bond/per/atom* or *extra/angle/per/atom*
    or *extra/dihedral/per/atom* or *extra/improper/per/atom* or
    *extra/special/per/atom*

        *bond/types* value = # of bond types
        *angle/types* value = # of angle types
        *dihedral/types* value = # of dihedral types
        *improper/types* value = # of improper types
        *extra/bond/per/atom* value = # of bonds per atom
        *extra/angle/per/atom* value = # of angles per atom
        *extra/dihedral/per/atom* value = # of dihedrals per atom
        *extra/improper/per/atom* value = # of impropers per atom
        *extra/special/per/atom* value = # of special neighbors per atom

## Examples

``` LAMMPS
create_box 2 mybox
create_box 2 mybox bond/types 2 extra/bond/per/atom 1
```

## Description

This command creates a simulation box based on the specified region.
Thus a [region](region) command must first be used to define a geometric
domain. It also partitions the simulation box into a regular 3d grid of
rectangular bricks, one per processor, based on the number of processors
being used and the settings of the [processors](processors) command. The
partitioning can later be changed by the [balance](balance) or [fix
balance](fix_balance) commands.

The argument N is the number of atom types that will be used in the
simulation.

If the region is not of style *prism*, then LAMMPS encloses the region
(block, sphere, etc.) with an axis-aligned orthogonal bounding box which
becomes the simulation domain.

If the region is of style *prism*, LAMMPS creates a non-orthogonal
simulation domain shaped as a parallelepiped with triclinic symmetry. As
defined by the [region prism](region) command, the parallelepiped has
its \"origin\" at (xlo,ylo,zlo) and is defined by three edge vectors
starting from the origin given by
$\vec a = (x_\text{hi}-x_\text{lo},0,0)$;
$\vec b = (xy,y_\text{hi}-y_\text{lo},0)$; and
$\vec c = (xz,yz,z_\text{hi}-z_\text{lo})$. The parameters *xy*, *xz*,
and *yz* can be 0.0 or positive or negative values and are called \"tilt
factors\" because they are the amount of displacement applied to faces
of an originally orthogonal box to transform it into the parallelepiped.

By default, a *prism* region used with the create_box command must have
tilt factors $(xy,xz,yz)$ that do not skew the box more than half the
distance of the parallel box length. For example, if $x_\text{lo} = 2$
and $x_\text{hi} = 12$, then the $x$ box length is 10 and the $xy$ tilt
factor must be between $-5$ and $5$. Similarly, both $xz$ and $yz$ must
be between $-(x_\text{hi}-x_\text{lo})/2$ and
$+(y_\text{hi}-y_\text{lo})/2$. Note that this is not a limitation,
since if the maximum tilt factor is 5 (as in this example), then
configurations with tilt $= \dots, -15$, $-5$, $5$, $15$, $25, \dots$
are all geometrically equivalent. Simulations with large tilt factors
will run inefficiently, since they require more ghost atoms and thus
more communication. With very large tilt factors, LAMMPS will eventually
produce incorrect trajectories and stop with errors due to lost atoms or
similar.

See the [Howto triclinic](Howto_triclinic) page for a geometric
description of triclinic boxes, as defined by LAMMPS, and how to
transform these parameters to and from other commonly used triclinic
representations.

When a prism region is used, the simulation domain should normally be
periodic in the dimension that the tilt is applied to, which is given by
the second dimension of the tilt factor (e.g., $y$ for $xy$ tilt). This
is so that pairs of atoms interacting across that boundary will have one
of them shifted by the tilt factor. Periodicity is set by the
[boundary](boundary) command. For example, if the $xy$ tilt factor is
non-zero, then the $y$ dimension should be periodic. Similarly, the $z$
dimension should be periodic if $xz$ or $yz$ is non-zero. LAMMPS does
not require this periodicity, but you may lose atoms if this is not the
case.

Note that if your simulation will tilt the box (e.g., via the [fix
deform](fix_deform) command), the simulation box must be set up to be
triclinic, even if the tilt factors are initially 0.0. You can also
change an orthogonal box to a triclinic box or vice versa by using the
[change box](change_box) command with its *ortho* and *triclinic*
options.

:::: note
::: title
Note
:::

If the system is non-periodic (in a dimension), then you should not make
the lo/hi box dimensions (as defined in your [region](region) command)
radically smaller/larger than the extent of the atoms you eventually
plan to create (e.g., via the [create_atoms](create_atoms) command). For
example, if your atoms extend from 0 to 50, you should not specify the
box bounds as $-10000$ and $10000$. This is because as described above,
LAMMPS uses the specified box size to lay out the 3d grid of processors.
A huge (mostly empty) box will be sub-optimal for performance when using
\"fixed\" boundary conditions (see the [boundary](boundary) command).
When using \"shrink-wrap\" boundary conditions (see the
[boundary](boundary) command), a huge (mostly empty) box may cause a
parallel simulation to lose atoms the first time that LAMMPS
shrink-wraps the box around the atoms.
::::

------------------------------------------------------------------------

The optional keywords can be used to create a system that allows for
bond (angle, dihedral, improper) interactions, or for molecules with
special 1\--2, 1\--3, or 1\--4 neighbors to be added later. These
optional keywords serve the same purpose as the analogous keywords that
can be used in a data file which are recognized by the
[read_data](read_data) command when it sets up a system.

Note that if these keywords are not used, then the create_box command
creates an atomic (non-molecular) simulation that does not allow bonds
between pairs of atoms to be defined, or a [bond potential](bond_style)
to be specified, or for molecules with special neighbors to be added to
the system by commands such as [create_atoms mol](create_atoms), [fix
deposit](fix_deposit) or [fix pour](fix_pour).

As an example, see the examples/deposit/in.deposit.molecule script,
which deposits molecules onto a substrate. Initially there are no
molecules in the system, but they are added later by the [fix
deposit](fix_deposit) command. The create_box command in the script uses
the bond/types and extra/bond/per/atom keywords to allow this. If the
added molecule contained more than one special bond (allowed by
default), an extra/special/per/atom keyword would also need to be
specified.

------------------------------------------------------------------------

## Restrictions

An [atom_style](atom_style) and [region](region) must have been
previously defined to use this command.

## Related commands

[read_data](read_data), [create_atoms](create_atoms), [region](region)

## Default

none
