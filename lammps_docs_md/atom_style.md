# atom_style command

## Syntax

``` LAMMPS
atom_style style args
```

-   style = *amoeba* or *angle* or *atomic* or *body* or *bond* or
    *charge* or *dielectric* or *dipole* or *dpd* or *edpd* or
    *electron* or *ellipsoid* or *full* or *line* or *mdpd* or
    *molecular* or *oxdna* or *peri* or *smd* or *sph* or *sphere* or
    *bpm/sphere* or *spin* or *tdpd* or *tri* or *template* or
    *wavepacket* or *hybrid*

        args = none for any style except the following
          *body* args = bstyle bstyle-args
            bstyle = style of body particles
            bstyle-args = additional arguments specific to the bstyle
                          see the `Howto body <Howto_body>`__ doc
                          page for details
          *sphere* arg = 0/1 (optional) for static/dynamic particle radii
          *bpm/sphere* arg = 0/1 (optional) for static/dynamic particle radii
          *tdpd* arg = Nspecies
            Nspecies = # of chemical species
          *template* arg = template-ID
            template-ID = ID of molecule template specified in a separate `molecule <molecule>`__ command
          *hybrid* args = list of one or more sub-styles, each with their args

-   accelerated styles (with same args) = *angle/kk* or *atomic/kk* or
    *bond/kk* or *charge/kk* or *full/kk* or *molecular/kk* or *spin/kk*

## Examples

``` LAMMPS
atom_style atomic
atom_style bond
atom_style full
atom_style body nparticle 2 10
atom_style hybrid charge bond
atom_style hybrid charge body nparticle 2 5
atom_style spin
atom_style template myMols
atom_style hybrid template twomols charge
atom_style tdpd 2
```

## Description

Define what style of atoms to use in a simulation. This determines what
attributes are associated with the atoms. This command must be used
before a simulation is setup via a [read_data](read_data),
[read_restart](read_restart), or [create_box](create_box) command.

:::: note
::: title
Note
:::

Many of the atom styles discussed here are only enabled if LAMMPS was
built with a specific package, as listed below in the Restrictions
section.
::::

Once a style is assigned, it cannot be changed, so use a style general
enough to encompass all attributes. E.g. with style *bond*, angular
terms cannot be used or added later to the model. It is OK to use a
style more general than needed, though it may be slightly inefficient.

The choice of style affects what quantities are stored by each atom,
what quantities are communicated between processors to enable forces to
be computed, and what quantities are listed in the data file read by the
[read_data](read_data) command.

These are the additional attributes of each style and the typical kinds
of physical systems they are used to model. All styles store
coordinates, velocities, atom IDs and types. See the
[read_data](read_data), [create_atoms](create_atoms), and [set](set)
commands for info on how to set these various quantities.

  -------------- ----------------------------------- -------------------------
  *amoeba*       molecular + charge + 1/5 neighbors  AMOEBA/HIPPO polarized
                                                     force fields

  *angle*        bonds and angles                    bead-spring polymers with
                                                     stiffness

  *atomic*       only the default values             coarse-grain liquids,
                                                     solids, metals

  *body*         mass, inertia moments, quaternion,  arbitrary bodies
                 angular momentum                    

  *bond*         bonds                               bead-spring polymers

  *charge*       charge                              atomic system with
                                                     charges

  *dielectric*   normx normy normz area/patch ed em  system with surface
                 epsilon curv                        polarization

  *dipole*       charge and dipole moment            system with dipolar
                                                     particles

  *dpd*          internal temperature and internal   DPD particles
                 energies                            

  *edpd*         temperature and heat capacity       eDPD particles

  *electron*     charge and spin and eradius         electronic force field

  *ellipsoid*    shape, quaternion, angular momentum aspherical particles

  *full*         molecular + charge                  bio-molecules

  *line*         end points, angular velocity        rigid bodies

  *mdpd*         density                             mDPD particles

  *molecular*    bonds, angles, dihedrals, impropers uncharged molecules

  *oxdna*        nucleotide polarity                 coarse-grained DNA and
                                                     RNA models

  *peri*         mass, volume                        mesoscopic Peridynamic
                                                     models

  *smd*          volume, kernel diameter, contact    solid and fluid SPH
                 radius, mass                        particles

  *sph*          rho, esph, cv                       SPH particles

  *sphere*       diameter, mass, angular velocity    granular models

  *bpm/sphere*   diameter, mass, angular velocity,   granular bonded particle
                 quaternion                          models (BPM)

  *spin*         magnetic moment                     system with magnetic
                                                     particles

  *tdpd*         chemical concentration              tDPD particles

  *template*     template index, template atom       small molecules with
                                                     fixed topology

  *tri*          corner points, angular momentum     rigid bodies

  *wavepacket*   charge, spin, eradius, etag, cs_re, AWPMD
                 cs_im                               
  -------------- ----------------------------------- -------------------------

:::: note
::: title
Note
:::

It is possible to add some attributes, such as a molecule ID, to atom
styles that do not have them via the [fix
property/atom](fix_property_atom) command. This command also allows new
custom attributes consisting of extra integer or floating-point values
to be added to atoms. See the [fix property/atom](fix_property_atom)
page for examples of cases where this is useful and details on how to
initialize, access, and output the custom values.
::::

All of the above styles define point particles, except the *sphere*,
*bpm/sphere*, *ellipsoid*, *electron*, *peri*, *wavepacket*, *line*,
*tri*, and *body* styles, which define finite-size particles. See the
[Howto spherical](Howto_spherical) page for an overview of using
finite-size particle models with LAMMPS.

All of the point-particle styles assign mass to particles on a per-type
basis, using the [mass](mass) command, The finite-size particle styles
assign mass to individual particles on a per-particle basis.

For the *sphere* and *bpm/sphere* styles, the particles are spheres and
each stores a per-particle diameter and mass. If the diameter \> 0.0,
the particle is a finite-size sphere. If the diameter = 0.0, it is a
point particle. Note that by use of the *disc* keyword with the [fix
nve/sphere](fix_nve_sphere), [fix nvt/sphere](fix_nvt_sphere), [fix
nph/sphere](fix_nph_sphere), [fix npt/sphere](fix_npt_sphere) commands
for the *sphere* style, spheres can be effectively treated as 2d discs
for a 2d simulation if desired. See also the [set density/disc](set)
command. These styles take an optional 0 or 1 argument. A value of 0
means the radius of each sphere is constant for the duration of the
simulation. A value of 1 means the radii may vary dynamically during the
simulation, e.g. due to use of the [fix adapt](fix_adapt) command.

For the *ellipsoid* style, the particles are ellipsoids and each stores
a flag which indicates whether it is a finite-size ellipsoid or a point
particle. If it is an ellipsoid, it also stores a shape vector with the
3 diameters of the ellipsoid and a quaternion 4-vector with its
orientation.

For the *dielectric* style, each particle can be either a physical
particle (e.g. an ion), or an interface particle representing a boundary
element between two regions of different dielectric constant. For
interface particles, in addition to the properties associated with
atom_style full, each particle also should be assigned a normal unit
vector (defined by normx, normy, normz), an area (area/patch), the
difference and mean of the dielectric constants of two sides of the
interface along the direction of the normal vector (ed and em), the
local dielectric constant at the boundary element (epsilon), and a mean
local curvature (curv). Physical particles must be assigned these
values, as well, but only their local dielectric constants will be used;
see documentation for associated [pair styles](pair_dielectric) and
[fixes](fix_polarize). The distinction between the physical and
interface particles is only meaningful when [fix polarize](fix_polarize)
commands are applied to the interface particles. This style is part of
the DIELECTRIC package.

For the *dipole* style, a point dipole is defined for each point
particle. Note that if you wish the particles to be finite-size spheres
as in a Stockmayer potential for a dipolar fluid, so that the particles
can rotate due to dipole-dipole interactions, then you need to use
atom_style hybrid sphere dipole, which will assign both a diameter and
dipole moment to each particle.

For the *electron* style, the particles representing electrons are 3d
Gaussians with a specified position and bandwidth or uncertainty in
position, which is represented by the eradius = electron size.

For the *peri* style, the particles are spherical and each stores a
per-particle mass and volume.

The *bpm/sphere* style is part of the BPM package.

The *oxdna* style is for coarse-grained nucleotides and stores the
3\'-to-5\' polarity of the nucleotide strand, which is set through the
bond topology in the data file. The first (second) atom in a bond
definition is understood to point towards the 3\'-end (5\'-end) of the
strand. Note that this style is part of the CG-DNA package.

The *dpd* style is for dissipative particle dynamics (DPD) particles.
Note that it is part of the DPD-REACT package, and is not for use with
the [pair_style dpd or dpd/stat](pair_dpd) commands, which can simply
use atom_style atomic. Atom_style dpd extends DPD particle properties
with internal temperature (dpdTheta), internal conductive energy
(uCond), internal mechanical energy (uMech), and internal chemical
energy (uChem).

The *edpd* style is for energy-conserving dissipative particle dynamics
(eDPD) particles which store a temperature (edpd_temp), and heat
capacity(edpd_cv).

The *mdpd* style is for many-body dissipative particle dynamics (mDPD)
particles which store a density (rho) for considering density-dependent
many-body interactions.

The *tdpd* style is for transport dissipative particle dynamics (tDPD)
particles which store a set of chemical concentration. An integer
\"cc_species\" is required to specify the number of chemical species
involved in a tDPD system.

The *sph* style is for smoothed particle hydrodynamics (SPH) particles
which store a density (rho), energy (esph), and heat capacity (cv).

The *smd* style is for a general formulation of Smooth Particle
Hydrodynamics. Both fluids and solids can be modeled. Particles store
the mass and volume of an integration point, a kernel diameter used for
calculating the field variables (e.g. stress and deformation) and a
contact radius for calculating repulsive forces which prevent individual
physical bodies from penetrating each other.

For the *spin* style, a magnetic spin is associated to each atom. Those
spins have a norm (their magnetic moment) and a direction.

The *wavepacket* style is similar to *electron*, but the electrons may
consist of several Gaussian wave packets, summed up with coefficients
cs= (cs_re,cs_im). Each of the wave packets is treated as a separate
particle in LAMMPS, wave packets belonging to the same electron must
have identical *etag* values.

For the *line* style, the particles are idealized line segments and each
stores a per-particle mass and length and orientation (i.e. the end
points of the line segment).

For the *tri* style, the particles are planar triangles and each stores
a per-particle mass and size and orientation (i.e. the corner points of
the triangle).

The *template* style allows molecular topology (bonds,angles,etc) to be
defined via a molecule template using the [molecule](molecule) command.
The template stores one or more molecules with a single copy of the
topology info (bonds,angles,etc) of each. Individual atoms only store a
template index and template atom to identify which molecule and which
atom-within-the-molecule they represent. Using the *template* style
instead of the *bond*, *angle*, *molecular* styles can save memory for
systems comprised of a large number of small molecules, all of a single
type (or small number of types). See the paper by Grime and Voth, in
[(Grime)](Grime), for examples of how this can be advantageous for
large-scale coarse-grained systems. The `examples/template` directory
has a few demo inputs and examples showing the use of the *template*
atom style versus *molecular*.

:::: note
::: title
Note
:::

When using the *template* style with a [molecule template](molecule)
that contains multiple molecules, you should ensure the atom types, bond
types, angle_types, etc in all the molecules are consistent. E.g. if one
molecule represents H2O and another CO2, then you probably do not want
each molecule file to define 2 atom types and a single bond type,
because they will conflict with each other when a mixture system of H2O
and CO2 molecules is defined, e.g. by the [read_data](read_data)
command. Rather the H2O molecule should define atom types 1 and 2, and
bond type 1. And the CO2 molecule should define atom types 3 and 4 (or
atom types 3 and 2 if a single oxygen type is desired), and bond type 2.
::::

For the *body* style, the particles are arbitrary bodies with internal
attributes defined by the \"style\" of the bodies, which is specified by
the *bstyle* argument. Body particles can represent complex entities,
such as surface meshes of discrete points, collections of sub-particles,
deformable objects, etc.

The [Howto body](Howto_body) page describes the body styles LAMMPS
currently supports, and provides more details as to the kind of body
particles they represent. For all styles, each body particle stores
moments of inertia and a quaternion 4-vector, so that its orientation
and position can be time integrated due to forces and torques.

Note that there may be additional arguments required along with the
*bstyle* specification, in the atom_style body command. These arguments
are described on the [Howto body](Howto_body) doc page.

------------------------------------------------------------------------

Typically, simulations require only a single (non-hybrid) atom style. If
some atoms in the simulation do not have all the properties defined by a
particular style, use the simplest style that defines all the needed
properties by any atom. For example, if some atoms in a simulation are
charged, but others are not, use the *charge* style. If some atoms have
bonds, but others do not, use the *bond* style.

The only scenario where the *hybrid* style is needed is if there is no
single style which defines all needed properties of all atoms. For
example, as mentioned above, if you want dipolar particles which will
rotate due to torque, you need to use \"atom_style hybrid sphere
dipole\". When a hybrid style is used, atoms store and communicate the
union of all quantities implied by the individual styles.

When using the *hybrid* style, you cannot combine the *template* style
with another molecular style that stores bond,angle,etc info on a
per-atom basis.

LAMMPS can be extended with new atom styles as well as new body styles;
see the [Modify](Modify) doc page.

------------------------------------------------------------------------

Styles with a *gpu*, *intel*, *kk*, *omp*, or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the [Accelerator packages](Speed_packages)
page. The accelerated styles take the same arguments and should produce
the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, INTEL, KOKKOS, OPENMP, and
OPT packages, respectively. They are only enabled if LAMMPS was built
with those packages. See the [Build package](Build_package) page for
more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the [-suffix command-line
switch](Run_options) when you invoke LAMMPS, or you can use the
[suffix](suffix) command in your input script.

See the [Accelerator packages](Speed_packages) page for more
instructions on how to use the accelerated styles effectively.

------------------------------------------------------------------------

## Restrictions

This command cannot be used after the simulation box is defined by a
[read_data](read_data) or [create_box](create_box) command.

Many of the styles listed above are only enabled if LAMMPS was built
with a specific package, as listed below. See the [Build
package](Build_package) page for more info.

The *amoeba* style is part of the AMOEBA package.

The *angle*, *bond*, *full*, *molecular*, and *template* styles are part
of the MOLECULE package.

The *line* and *tri* styles are part of the ASPHERE package.

The *body* style is part of the BODY package.

The *dipole* style is part of the DIPOLE package.

The *peri* style is part of the PERI package for Peridynamics.

The *oxdna* style is part of the CG-DNA package for coarse-grained
simulation of DNA and RNA.

The *electron* style is part of the EFF package for [electronic force
fields](pair_eff).

The *dpd* style is part of the DPD-REACT package for dissipative
particle dynamics (DPD).

The *edpd*, *mdpd*, and *tdpd* styles are part of the DPD-MESO package
for energy-conserving dissipative particle dynamics (eDPD), many-body
dissipative particle dynamics (mDPD), and transport dissipative particle
dynamics (tDPD), respectively.

The *sph* style is part of the SPH package for smoothed particle
hydrodynamics (SPH). See [this PDF
guide](PDF/SPH_LAMMPS_userguide.pdf)\_ to using SPH in LAMMPS.

The *spin* style is part of the SPIN package.

The *wavepacket* style is part of the AWPMD package for the
[antisymmetrized wave packet MD method](pair_awpmd).

## Related commands

[read_data](read_data), [pair_style](pair_style)

## Default

The default atom style is atomic. If atom_style sphere is used its
default argument is 0.

------------------------------------------------------------------------

::: {#Grime}
**(Grime)** Grime and Voth, to appear in J Chem Theory & Computation
(2014).
:::
