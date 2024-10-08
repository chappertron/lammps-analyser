# fix rigid command

Accelerator Variants: *rigid/omp*

# fix rigid/nve command

Accelerator Variants: *rigid/nve/omp*

# fix rigid/nvt command

Accelerator Variants: *rigid/nvt/omp*

# fix rigid/npt command

Accelerator Variants: *rigid/npt/omp*

# fix rigid/nph command

Accelerator Variants: *rigid/nph/omp*

# fix rigid/small command

Accelerator Variants: *rigid/small/omp*

# fix rigid/nve/small command

# fix rigid/nvt/small command

# fix rigid/npt/small command

# fix rigid/nph/small command

## Syntax

    fix ID group-ID style bodystyle args keyword values ...

-   ID, group-ID are documented in [fix](fix) command

-   style = *rigid* or *rigid/nve* or *rigid/nvt* or *rigid/npt* or
    *rigid/nph* or *rigid/small* or *rigid/nve/small* or
    *rigid/nvt/small* or *rigid/npt/small* or *rigid/nph/small*

-   bodystyle = *single* or *molecule* or *group*

        *single* args = none
        *molecule* args = none
        *custom* args = *i_propname* or *v_varname*
          i_propname = a custom integer vector defined via fix property/atom
          v_varname  = an atom-style or atomfile-style variable
        *group* args = N groupID1 groupID2 ...
          N = # of groups
          groupID1, groupID2, ... = list of N group IDs

-   zero or more keyword/value pairs may be appended

-   keyword = *langevin* or *reinit* or *temp* or *iso* or *aniso* or
    *x* or *y* or *z* or *couple* or *tparam* or *pchain* or *dilate* or
    *force* or *torque* or *infile* or *gravity*

        *langevin* values = Tstart Tstop Tperiod seed
          Tstart,Tstop = desired temperature at start/stop of run (temperature units)
          Tdamp = temperature damping parameter (time units)
          seed = random number seed to use for white noise (positive integer)
        *reinit* value = *yes* or *no*
        *temp* values = Tstart Tstop Tdamp
          Tstart,Tstop = desired temperature at start/stop of run (temperature units)
          Tdamp = temperature damping parameter (time units)
        *iso* or *aniso* values = Pstart Pstop Pdamp
          Pstart,Pstop = scalar external pressure at start/end of run (pressure units)
          Pdamp = pressure damping parameter (time units)
        *x* or *y* or *z* values = Pstart Pstop Pdamp
          Pstart,Pstop = external stress tensor component at start/end of run (pressure units)
          Pdamp = stress damping parameter (time units)
        *couple* value = *none* or *xyz* or *xy* or *yz* or *xz*
        *tparam* values = Tchain Titer Torder
          Tchain = length of Nose/Hoover thermostat chain
          Titer = number of thermostat iterations performed
          Torder = 3 or 5 = Yoshida-Suzuki integration parameters
        *pchain* values = Pchain
          Pchain = length of the Nose/Hoover thermostat chain coupled with the barostat
        *dilate* value = dilate-group-ID
          dilate-group-ID = only dilate atoms in this group due to barostat volume changes
        *force* values = M xflag yflag zflag
          M = which rigid body from 1-Nbody (see asterisk form below)
          xflag,yflag,zflag = off/on if component of center-of-mass force is active
        *torque* values = M xflag yflag zflag
          M = which rigid body from 1-Nbody (see asterisk form below)
          xflag,yflag,zflag = off/on if component of center-of-mass torque is active
        *infile* filename
          filename = file with per-body values of mass, center-of-mass, moments of inertia
        *gravity* values = gravity-ID
          gravity-ID = ID of fix gravity command to add gravitational forces

## Examples

``` LAMMPS
fix 1 clump rigid single reinit yes
fix 1 clump rigid/small molecule
fix 1 clump rigid single force 1 off off on langevin 1.0 1.0 1.0 428984
fix 1 polychains rigid/nvt molecule temp 1.0 1.0 5.0 reinit no
fix 1 polychains rigid molecule force 1*5 off off off force 6*10 off off on
fix 1 polychains rigid/small molecule langevin 1.0 1.0 1.0 428984
fix 2 fluid rigid group 3 clump1 clump2 clump3 torque * off off off
fix 1 rods rigid/npt molecule temp 300.0 300.0 100.0 iso 0.5 0.5 10.0
fix 1 particles rigid/npt molecule temp 1.0 1.0 5.0 x 0.5 0.5 1.0 z 0.5 0.5 1.0 couple xz
fix 1 water rigid/nph molecule iso 0.5 0.5 1.0
fix 1 particles rigid/npt/small molecule temp 1.0 1.0 1.0 iso 0.5 0.5 1.0

variable bodyid atom 1.0*gmask(clump1)+2.0*gmask(clump2)+3.0*gmask(clump3)
fix 1 clump rigid custom v_bodyid

variable bodyid atomfile bodies.txt
fix 1 clump rigid custom v_bodyid

fix 0 all property/atom i_bodyid
read_restart data.rigid fix 0 NULL Bodies
fix 1 clump rigid/small custom i_bodyid
```

## Description

Treat one or more sets of atoms as independent rigid bodies. This means
that each timestep the total force and torque on each rigid body is
computed as the sum of the forces and torques on its constituent
particles. The coordinates, velocities, and orientations of the atoms in
each body are then updated so that the body moves and rotates as a
single entity. This is implemented by creating internal data structures
for each rigid body and performing time integration on these data
structures. Positions, velocities, and orientations of the constituent
particles are regenerated from the rigid body data structures in every
time step. This restricts which operations and fixes can be applied to
rigid bodies. See below for a detailed discussion.

Examples of large rigid bodies are a colloidal particle, or portions of
a biomolecule such as a protein.

Example of small rigid bodies are patchy nanoparticles, such as those
modeled in [this paper](Zhang1) by Sharon Glotzer\'s group, clumps of
granular particles, lipid molecules consisting of one or more point
dipoles connected to other spheroids or ellipsoids, irregular particles
built from line segments (2d) or triangles (3d), and coarse-grain models
of nano or colloidal particles consisting of a small number of
constituent particles. Note that the [fix shake](fix_shake) command can
also be used to rigidify small molecules of 2, 3, or 4 atoms, e.g. water
molecules. That fix treats the constituent atoms as point masses.

These fixes also update the positions and velocities of the atoms in
each rigid body via time integration, in the NVE, NVT, NPT, or NPH
ensemble, as described below.

There are two main variants of this fix, fix rigid and fix rigid/small.
The NVE/NVT/NPT/NHT versions belong to one of the two variants, as their
style names indicate.

:::: note
::: title
Note
:::

Not all of the *bodystyle* options and keyword/value options are
available for both the *rigid* and *rigid/small* variants. See details
below.
::::

The *rigid* styles are typically the best choice for a system with a
small number of large rigid bodies, each of which can extend across the
domain of many processors. It operates by creating a single global list
of rigid bodies, which all processors contribute to. MPI_Allreduce
operations are performed each timestep to sum the contributions from
each processor to the force and torque on all the bodies. This operation
will not scale well in parallel if large numbers of rigid bodies are
simulated.

The *rigid/small* styles are typically best for a system with a large
number of small rigid bodies. Each body is assigned to the atom closest
to the geometrical center of the body. The fix operates using local
lists of rigid bodies owned by each processor and information is
exchanged and summed via local communication between neighboring
processors when ghost atom info is accumulated.

:::: note
::: title
Note
:::

To use the *rigid/small* styles the ghost atom cutoff must be large
enough to span the distance between the atom that owns the body and
every other atom in the body. This distance value is printed out when
the rigid bodies are defined. If the [pair_style](pair_style) cutoff
plus neighbor skin does not span this distance, then you should use the
[comm_modify cutoff](comm_modify) command with a setting epsilon larger
than the distance.
::::

Which of the two variants is faster for a particular problem is hard to
predict. The best way to decide is to perform a short test run. Both
variants should give identical numerical answers for short runs. Long
runs should give statistically similar results, but round-off
differences may accumulate to produce divergent trajectories.

:::: note
::: title
Note
:::

You should not update the atoms in rigid bodies via other
time-integration fixes (e.g. [fix nve](fix_nve), [fix nvt](fix_nh), [fix
npt](fix_nh), [fix move](fix_move)), or you will have conflicting
updates to positions and velocities resulting in unphysical behavior in
most cases. When performing a hybrid simulation with some atoms in rigid
bodies, and some not, a separate time integration fix like [fix
nve](fix_nve) or [fix nvt](fix_nh) should be used for the non-rigid
particles.
::::

:::: note
::: title
Note
:::

These fixes are overkill if you simply want to hold a collection of
atoms stationary or have them move with a constant velocity. A simpler
way to hold atoms stationary is to not include those atoms in your time
integration fix. E.g. use \"fix 1 mobile nve\" instead of \"fix 1 all
nve\", where \"mobile\" is the group of atoms that you want to move. You
can move atoms with a constant velocity by assigning them an initial
velocity (via the [velocity](velocity) command), setting the force on
them to 0.0 (via the [fix setforce](fix_setforce) command), and
integrating them as usual (e.g. via the [fix nve](fix_nve) command).
::::

:::: warning
::: title
Warning
:::

The aggregate properties of each rigid body are calculated at the start
of a simulation run and are maintained in internal data structures. The
properties include the position and velocity of the center-of-mass of
the body, its moments of inertia, and its angular momentum. This is done
using the properties of the constituent atoms of the body at that point
in time (or see the *infile* keyword option). Thereafter, changing these
properties of individual atoms in the body will have no effect on a
rigid body\'s dynamics, unless they effect any computation of per-atom
forces or torques. If the keyword *reinit* is set to *yes* (the
default), the rigid body data structures will be recreated at the
beginning of each *run* command; if the keyword *reinit* is set to *no*,
the rigid body data structures will be built only at the very first
*run* command and maintained for as long as the rigid fix is defined.
For example, you might think you could displace the atoms in a body or
add a large velocity to each atom in a body to make it move in a desired
direction before a second run is performed, using the [set](set) or
[displace_atoms](displace_atoms) or [velocity](velocity) commands. But
these commands will not affect the internal attributes of the body
unless *reinit* is set to *yes*. With *reinit* set to *no* (or using the
*infile* option, which implies *reinit* *no*) the position and velocity
of individual atoms in the body will be reset when time integration
starts again.
::::

------------------------------------------------------------------------

Each rigid body must have two or more atoms. An atom can belong to at
most one rigid body. Which atoms are in which bodies can be defined via
several options.

:::: note
::: title
Note
:::

With the *rigid/small* styles, which require that *bodystyle* be
specified as *molecule* or *custom*, you can define a system that has no
rigid bodies initially. This is useful when you are using the *mol*
keyword in conjunction with another fix that is adding rigid bodies
on-the-fly as molecules, such as [fix deposit](fix_deposit) or [fix
pour](fix_pour).
::::

For bodystyle *single* the entire fix group of atoms is treated as one
rigid body. This option is only allowed for the *rigid* styles.

For bodystyle *molecule*, atoms are grouped into rigid bodies by their
respective molecule IDs: each set of atoms in the fix group with the
same molecule ID is treated as a different rigid body. This option is
allowed for both the *rigid* and *rigid/small* styles. Note that atoms
with a molecule ID = 0 will be treated as a single rigid body. For a
system with atomic solvent (typically this is atoms with molecule ID =
0) surrounding rigid bodies, this may not be what you want. Thus you
should be careful to use a fix group that only includes atoms you want
to be part of rigid bodies.

Bodystyle *custom* is similar to bodystyle *molecule* except that it is
more flexible in using other per-atom properties to define the sets of
atoms that form rigid bodies. A custom per-atom integer vector defined
by the [fix property/atom](fix_property_atom) command can be used. Or an
[atom-style or atomfile-style variable](variable) can be used; the
floating-point value produced by the variable is rounded to an integer.
As with bodystyle *molecule*, each set of atoms in the fix groups with
the same integer value is treated as a different rigid body. Since fix
property/atom custom vectors and atom-style variables produce values for
all atoms, you should be careful to use a fix group that only includes
atoms you want to be part of rigid bodies.

:::: note
::: title
Note
:::

To compute the initial center-of-mass position and other properties of
each rigid body, the image flags for each atom in the body are used to
\"unwrap\" the atom coordinates. Thus you must ensure that these image
flags are consistent so that the unwrapping creates a valid rigid body
(one where the atoms are close together), particularly if the atoms in a
single rigid body straddle a periodic boundary. This means the input
data file or restart file must define the image flags for each atom
consistently or that you have used the [set](set) command to specify
them correctly. If a dimension is non-periodic then the image flag of
each atom must be 0 in that dimension, else an error is generated.
::::

The *force* and *torque* keywords discussed next are only allowed for
the *rigid* styles.

By default, each rigid body is acted on by other atoms which induce an
external force and torque on its center of mass, causing it to translate
and rotate. Components of the external center-of-mass force and torque
can be turned off by the *force* and *torque* keywords. This may be
useful if you wish a body to rotate but not translate, or vice versa, or
if you wish it to rotate or translate continuously unaffected by
interactions with other particles. Note that if you expect a rigid body
not to move or rotate by using these keywords, you must ensure its
initial center-of-mass translational or angular velocity is 0.0.
Otherwise the initial translational or angular momentum the body has
will persist.

An xflag, yflag, or zflag set to *off* means turn off the component of
force of torque in that dimension. A setting of *on* means turn on the
component, which is the default. Which rigid body(s) the settings apply
to is determined by the first argument of the *force* and *torque*
keywords. It can be an integer M from 1 to Nbody, where Nbody is the
number of rigid bodies defined. A wild-card asterisk can be used in
place of, or in conjunction with, the M argument to set the flags for
multiple rigid bodies. This takes the form \"\*\" or \"\*n\" or \"n\*\"
or \"m\*n\". If N = the number of rigid bodies, then an asterisk with no
numeric values means all bodies from 1 to N. A leading asterisk means
all bodies from 1 to n (inclusive). A trailing asterisk means all bodies
from n to N (inclusive). A middle asterisk means all types from m to n
(inclusive). Note that you can use the *force* or *torque* keywords as
many times as you like. If a particular rigid body has its component
flags set multiple times, the settings from the final keyword are used.

:::: note
::: title
Note
:::

For computational efficiency, you may wish to turn off pairwise and bond
interactions within each rigid body, as they no longer contribute to the
motion. The [neigh_modify exclude](neigh_modify) and
[delete_bonds](delete_bonds) commands are used to do this. If the rigid
bodies have strongly overlapping atoms, you may need to turn off these
interactions to avoid numerical problems due to large equal/opposite
intra-body forces swamping the contribution of small inter-body forces.
::::

For computational efficiency, you should typically define one fix rigid
or fix rigid/small command which includes all the desired rigid bodies.
LAMMPS will allow multiple rigid fixes to be defined, but it is more
expensive.

------------------------------------------------------------------------

The constituent particles within a rigid body can be point particles
(the default in LAMMPS) or finite-size particles, such as spheres or
ellipsoids or line segments or triangles. See the [atom_style sphere and
ellipsoid and line and tri](atom_style) commands for more details on
these kinds of particles. Finite-size particles contribute differently
to the moment of inertia of a rigid body than do point particles.
Finite-size particles can also experience torque (e.g. due to
[frictional granular interactions](pair_gran)) and have an orientation.
These contributions are accounted for by these fixes.

Forces between particles within a body do not contribute to the external
force or torque on the body. Thus for computational efficiency, you may
wish to turn off pairwise and bond interactions between particles within
each rigid body. The [neigh_modify exclude](neigh_modify) and
[delete_bonds](delete_bonds) commands are used to do this. For
finite-size particles this also means the particles can be highly
overlapped when creating the rigid body.

------------------------------------------------------------------------

The *rigid*, *rigid/nve*, *rigid/small*, and *rigid/small/nve* styles
perform constant NVE time integration. They are referred to below as the
4 NVE rigid styles. The only difference is that the *rigid* and
*rigid/small* styles use an integration technique based on Richardson
iterations. The *rigid/nve* and *rigid/small/nve* styles uses the
methods described in the paper by [Miller](Miller3), which are thought
to provide better energy conservation than an iterative approach.

The *rigid/nvt* and *rigid/nvt/small* styles performs constant NVT
integration using a Nose/Hoover thermostat with chains as described
originally in [(Hoover)](Hoover) and [(Martyna)](Martyna2), which
thermostats both the translational and rotational degrees of freedom of
the rigid bodies. They are referred to below as the 2 NVT rigid styles.
The rigid-body algorithm used by *rigid/nvt* is described in the paper
by [Kamberaj](Kamberaj).

The *rigid/npt*, *rigid/nph*, *rigid/npt/small*, and *rigid/nph/small*
styles perform constant NPT or NPH integration using a Nose/Hoover
barostat with chains. They are referred to below as the 4 NPT and NPH
rigid styles. For the NPT case, the same Nose/Hoover thermostat is also
used as with *rigid/nvt* and *rigid/nvt/small*.

The barostat parameters are specified using one or more of the *iso*,
*aniso*, *x*, *y*, *z* and *couple* keywords. These keywords give you
the ability to specify 3 diagonal components of the external stress
tensor, and to couple these components together so that the dimensions
they represent are varied together during a constant-pressure
simulation. The effects of these keywords are similar to those defined
in [fix npt/nph](fix_nh)

:::: note
::: title
Note
:::

Currently the *rigid/npt*, *rigid/nph*, *rigid/npt/small*, and
*rigid/nph/small* styles do not support triclinic (non-orthogonal)
boxes.
::::

The target pressures for each of the 6 components of the stress tensor
can be specified independently via the *x*, *y*, *z* keywords, which
correspond to the 3 simulation box dimensions. For each component, the
external pressure or tensor component at each timestep is a ramped value
during the run from *Pstart* to *Pstop*. If a target pressure is
specified for a component, then the corresponding box dimension will
change during a simulation. For example, if the *y* keyword is used, the
y-box length will change. A box dimension will not change if that
component is not specified, although you have the option to change that
dimension via the [fix deform](fix_deform) command.

For all barostat keywords, the *Pdamp* parameter operates like the
*Tdamp* parameter, determining the time scale on which pressure is
relaxed. For example, a value of 10.0 means to relax the pressure in a
timespan of (roughly) 10 time units (e.g. $\tau$ or fs or ps - see the
[units](units) command).

Regardless of what atoms are in the fix group (the only atoms which are
time integrated), a global pressure or stress tensor is computed for all
atoms. Similarly, when the size of the simulation box is changed, all
atoms are re-scaled to new positions, unless the keyword *dilate* is
specified with a *dilate-group-ID* for a group that represents a subset
of the atoms. This can be useful, for example, to leave the coordinates
of atoms in a solid substrate unchanged and controlling the pressure of
a surrounding fluid. Another example is a system consisting of rigid
bodies and point particles where the barostat is only coupled with the
rigid bodies. This option should be used with care, since it can be
unphysical to dilate some atoms and not others, because it can introduce
large, instantaneous displacements between a pair of atoms (one dilated,
one not) that are far from the dilation origin.

The *couple* keyword allows two or three of the diagonal components of
the pressure tensor to be \"coupled\" together. The value specified with
the keyword determines which are coupled. For example, *xz* means the
*Pxx* and *Pzz* components of the stress tensor are coupled. *Xyz* means
all 3 diagonal components are coupled. Coupling means two things: the
instantaneous stress will be computed as an average of the corresponding
diagonal components, and the coupled box dimensions will be changed
together in lockstep, meaning coupled dimensions will be dilated or
contracted by the same percentage every timestep. The *Pstart*, *Pstop*,
*Pdamp* parameters for any coupled dimensions must be identical. *Couple
xyz* can be used for a 2d simulation; the *z* dimension is simply
ignored.

The *iso* and *aniso* keywords are simply shortcuts that are equivalent
to specifying several other keywords together.

The keyword *iso* means couple all 3 diagonal components together when
pressure is computed (hydrostatic pressure), and dilate/contract the
dimensions together. Using \"iso Pstart Pstop Pdamp\" is the same as
specifying these 4 keywords:

    x Pstart Pstop Pdamp
    y Pstart Pstop Pdamp
    z Pstart Pstop Pdamp
    couple xyz

The keyword *aniso* means *x*, *y*, and *z* dimensions are controlled
independently using the *Pxx*, *Pyy*, and *Pzz* components of the stress
tensor as the driving forces, and the specified scalar external
pressure. Using \"aniso Pstart Pstop Pdamp\" is the same as specifying
these 4 keywords:

    x Pstart Pstop Pdamp
    y Pstart Pstop Pdamp
    z Pstart Pstop Pdamp
    couple none

------------------------------------------------------------------------

The keyword/value option pairs are used in the following ways.

The *reinit* keyword determines, whether the rigid body properties are
re-initialized between run commands. With the option *yes* (the default)
this is done, with the option *no* this is not done. Turning off the
re-initialization can be helpful to protect rigid bodies against
unphysical manipulations between runs or when properties cannot be
easily re-computed (e.g. when read from a file). When using the *infile*
keyword, the *reinit* option is automatically set to *no*.

The *langevin* and *temp* and *tparam* keywords perform thermostatting
of the rigid bodies, altering both their translational and rotational
degrees of freedom. What is meant by \"temperature\" of a collection of
rigid bodies and how it can be monitored via the fix output is discussed
below.

The *langevin* keyword applies a Langevin thermostat to the constant NVE
time integration performed by any of the 4 NVE rigid styles: *rigid*,
*rigid/nve*, *rigid/small*, *rigid/small/nve*. It cannot be used with
the 2 NVT rigid styles: *rigid/nvt*, *rigid/small/nvt*. The desired
temperature at each timestep is a ramped value during the run from
*Tstart* to *Tstop*. The *Tdamp* parameter is specified in time units
and determines how rapidly the temperature is relaxed. For example, a
value of 100.0 means to relax the temperature in a timespan of (roughly)
100 time units ($\tau$ or fs or ps - see the [units](units) command).
The random \# *seed* must be a positive integer.

The way that Langevin thermostatting operates is explained on the [fix
langevin](fix_langevin) doc page. If you wish to simply viscously damp
the rotational motion without thermostatting, you can set *Tstart* and
*Tstop* to 0.0, which means only the viscous drag term in the Langevin
thermostat will be applied. See the discussion on the [fix
viscous](fix_viscous) page for details.

:::: note
::: title
Note
:::

When the *langevin* keyword is used with fix rigid versus fix
rigid/small, different dynamics will result for parallel runs. This is
because of the way random numbers are used in the two cases. The
dynamics for the two cases should be statistically similar, but will not
be identical, even for a single timestep.
::::

The *temp* and *tparam* keywords apply a Nose/Hoover thermostat to the
NVT time integration performed by the 2 NVT rigid styles. They cannot be
used with the 4 NVE rigid styles. The desired temperature at each
timestep is a ramped value during the run from *Tstart* to *Tstop*. The
*Tdamp* parameter is specified in time units and determines how rapidly
the temperature is relaxed. For example, a value of 100.0 means to relax
the temperature in a timespan of (roughly) 100 time units (tau or fs or
ps - see the [units](units) command).

Nose/Hoover chains are used in conjunction with this thermostat. The
*tparam* keyword can optionally be used to change the chain settings
used. *Tchain* is the number of thermostats in the Nose Hoover chain.
This value, along with *Tdamp* can be varied to dampen undesirable
oscillations in temperature that can occur in a simulation. As a rule of
thumb, increasing the chain length should lead to smaller oscillations.
The keyword *pchain* specifies the number of thermostats in the chain
thermostatting the barostat degrees of freedom.

:::: note
::: title
Note
:::

There are alternate ways to thermostat a system of rigid bodies. You can
use [fix langevin](fix_langevin) to treat the individual particles in
the rigid bodies as effectively immersed in an implicit solvent, e.g. a
Brownian dynamics model. For hybrid systems with both rigid bodies and
solvent particles, you can thermostat only the solvent particles that
surround one or more rigid bodies by appropriate choice of groups in the
compute and fix commands for temperature and thermostatting. The solvent
interactions with the rigid bodies should then effectively thermostat
the rigid body temperature as well without use of the Langevin or
Nose/Hoover options associated with the fix rigid commands.
::::

------------------------------------------------------------------------

The *mol* keyword can only be used with the *rigid/small* styles. It
must be used when other commands, such as [fix deposit](fix_deposit) or
[fix pour](fix_pour), add rigid bodies on-the-fly during a simulation.
You specify a *template-ID* previously defined using the
[molecule](molecule) command, which reads a file that defines the
molecule. You must use the same *template-ID* that the other fix which
is adding rigid bodies uses. The coordinates, atom types, atom
diameters, center-of-mass, and moments of inertia can be specified in
the molecule file. See the [molecule](molecule) command for details. The
only settings required to be in this file are the coordinates and types
of atoms in the molecule, in which case the molecule command calculates
the other quantities itself.

Note that these other fixes create new rigid bodies, in addition to
those defined initially by this fix via the *bodystyle* setting.

Also note that when using the *mol* keyword, extra restart information
about all rigid bodies is written out whenever a restart file is written
out. See the NOTE in the next section for details.

------------------------------------------------------------------------

The *infile* keyword allows a file of rigid body attributes to be read
in from a file, rather then having LAMMPS compute them. There are 5 such
attributes: the total mass of the rigid body, its center-of-mass
position, its 6 moments of inertia, its center-of-mass velocity, and the
3 image flags of the center-of-mass position. For rigid bodies
consisting of point particles or non-overlapping finite-size particles,
LAMMPS can compute these values accurately. However, for rigid bodies
consisting of finite-size particles which overlap each other, LAMMPS
will ignore the overlaps when computing these 4 attributes. The amount
of error this induces depends on the amount of overlap. To avoid this
issue, the values can be pre-computed (e.g. using Monte Carlo
integration).

The format of the file is as follows. Note that the file does not have
to list attributes for every rigid body integrated by fix rigid. Only
bodies which the file specifies will have their computed attributes
overridden. The file can contain initial blank lines or comment lines
starting with \"#\" which are ignored. The first non-blank, non-comment
line should list N = the number of lines to follow. The N successive
lines contain the following information:

    ID1 masstotal xcm ycm zcm ixx iyy izz ixy ixz iyz vxcm vycm vzcm lx ly lz ixcm iycm izcm
    ID2 masstotal xcm ycm zcm ixx iyy izz ixy ixz iyz vxcm vycm vzcm lx ly lz ixcm iycm izcm
    ...
    IDN masstotal xcm ycm zcm ixx iyy izz ixy ixz iyz vxcm vycm vzcm lx ly lz ixcm iycm izcm

The rigid body IDs are all positive integers. For the *single*
bodystyle, only an ID of 1 can be used. For the *group* bodystyle, IDs
from 1 to Ng can be used where Ng is the number of specified groups. For
the *molecule* bodystyle, use the molecule ID for the atoms in a
specific rigid body as the rigid body ID.

The masstotal and center-of-mass coordinates (xcm,ycm,zcm) are
self-explanatory. The center-of-mass should be consistent with what is
calculated for the position of the rigid body with all its atoms
unwrapped by their respective image flags. If this produces a
center-of-mass that is outside the simulation box, LAMMPS wraps it back
into the box.

The 6 moments of inertia (ixx,iyy,izz,ixy,ixz,iyz) should be the values
consistent with the current orientation of the rigid body around its
center of mass. The values are with respect to the simulation box XYZ
axes, not with respect to the principal axes of the rigid body itself.
LAMMPS performs the latter calculation internally.

The (vxcm,vycm,vzcm) values are the velocity of the center of mass. The
(lx,ly,lz) values are the angular momentum of the body. The
(vxcm,vycm,vzcm) and (lx,ly,lz) values can simply be set to 0 if you
wish the body to have no initial motion.

The (ixcm,iycm,izcm) values are the image flags of the center of mass of
the body. For periodic dimensions, they specify which image of the
simulation box the body is considered to be in. An image of 0 means it
is inside the box as defined. A value of 2 means add 2 box lengths to
get the true value. A value of -1 means subtract 1 box length to get the
true value. LAMMPS updates these flags as the rigid bodies cross
periodic boundaries during the simulation.

:::: note
::: title
Note
:::

If you use the *infile* or *mol* keywords and write restart files during
a simulation, then each time a restart file is written, the fix also
write an auxiliary restart file with the name rfile.rigid, where
\"rfile\" is the name of the restart file, e.g. tmp.restart.10000 and
tmp.restart.10000.rigid. This auxiliary file is in the same format
described above. Thus it can be used in a new input script that restarts
the run and re-specifies a rigid fix using an *infile* keyword and the
appropriate filename. Note that the auxiliary file will contain one line
for every rigid body, even if the original file only listed a subset of
the rigid bodies.
::::

------------------------------------------------------------------------

If you use a [temperature compute](compute) with a group that includes
particles in rigid bodies, the degrees-of-freedom removed by each rigid
body are accounted for in the temperature (and pressure) computation,
but only if the temperature group includes all the particles in a
particular rigid body.

A 3d rigid body has 6 degrees of freedom (3 translational, 3
rotational), except for a collection of point particles lying on a
straight line, which has only 5, e.g a dimer. A 2d rigid body has 3
degrees of freedom (2 translational, 1 rotational).

:::: note
::: title
Note
:::

You may wish to explicitly subtract additional degrees-of-freedom if you
use the *force* and *torque* keywords to eliminate certain motions of
one or more rigid bodies. LAMMPS does not do this automatically.
::::

The rigid body contribution to the pressure of the system (virial) is
also accounted for by this fix.

------------------------------------------------------------------------

If your simulation is a hybrid model with a mixture of rigid bodies and
non-rigid particles (e.g. solvent) there are several ways these rigid
fixes can be used in tandem with [fix nve](fix_nve), [fix nvt](fix_nh),
[fix npt](fix_nh), and [fix nph](fix_nh).

If you wish to perform NVE dynamics (no thermostatting or barostatting),
use one of 4 NVE rigid styles to integrate the rigid bodies, and [fix
nve](fix_nve) to integrate the non-rigid particles.

If you wish to perform NVT dynamics (thermostatting, but no
barostatting), you can use one of the 2 NVT rigid styles for the rigid
bodies, and any thermostatting fix for the non-rigid particles ([fix
nvt](fix_nh), [fix langevin](fix_langevin), [fix
temp/berendsen](fix_temp_berendsen)). You can also use one of the 4 NVE
rigid styles for the rigid bodies and thermostat them using [fix
langevin](fix_langevin) on the group that contains all the particles in
the rigid bodies. The net force added by [fix langevin](fix_langevin) to
each rigid body effectively thermostats its translational center-of-mass
motion. Not sure how well it does at thermostatting its rotational
motion.

If you wish to perform NPT or NPH dynamics (barostatting), you cannot
use both [fix npt](fix_nh) and the NPT or NPH rigid styles. This is
because there can only be one fix which monitors the global pressure and
changes the simulation box dimensions. So you have 3 choices:

1.  Use one of the 4 NPT or NPH styles for the rigid bodies. Use the
    *dilate* all option so that it will dilate the positions of the
    non-rigid particles as well. Use [fix nvt](fix_nh) (or any other
    thermostat) for the non-rigid particles.
2.  Use [fix npt](fix_nh) for the group of non-rigid particles. Use the
    *dilate* all option so that it will dilate the center-of-mass
    positions of the rigid bodies as well. Use one of the 4 NVE or 2 NVT
    rigid styles for the rigid bodies.
3.  Use [fix press/berendsen](fix_press_berendsen) to compute the
    pressure and change the box dimensions. Use one of the 4 NVE or 2
    NVT rigid styles for the rigid bodies. Use [fix nvt](fix_nh) (or any
    other thermostat) for the non-rigid particles.

In all case, the rigid bodies and non-rigid particles both contribute to
the global pressure and the box is scaled the same by any of the
barostatting fixes.

You could even use the second and third options for a non-hybrid
simulation consisting of only rigid bodies, assuming you give [fix
npt](fix_nh) an empty group, though it\'s an odd thing to do. The
barostatting fixes ([fix npt](fix_nh) and [fix
press/berensen](fix_press_berendsen)) will monitor the pressure and
change the box dimensions, but not time integrate any particles. The
integration of the rigid bodies will be performed by fix rigid/nvt.

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

## Restart, fix_modify, output, run start/stop, minimize info

No information about the 4 NVE rigid styles is written to [binary
restart files](restart). The exception is if the *infile* or *mol*
keyword is used, in which case an auxiliary file is written out with
rigid body information each time a restart file is written, as explained
above for the *infile* keyword. For the 2 NVT rigid styles, the state of
the Nose/Hoover thermostat is written to [binary restart
files](restart). Ditto for the 4 NPT and NPH rigid styles, and the state
of the Nose/Hoover barostat. See the [read_restart](read_restart)
command for info on how to re-specify a fix in an input script that
reads a restart file, so that the operation of the fix continues in an
uninterrupted fashion.

The [fix_modify](fix_modify) *temp* and *press* options are supported by
the 4 NPT and NPH rigid styles to change the computes used to calculate
the instantaneous pressure tensor. Note that the 2 NVT rigid fixes do
not use any external compute to compute instantaneous temperature.

The [fix_modify](fix_modify) *bodyforces* option is supported by all
rigid styles to set whether per-body forces and torques are computed
early or late in a timestep, i.e. at the post-force stage or at the
final-integrate stage or the timestep, respectively.

The cumulative energy change in the system imposed by the 6 NVT, NPT,
NPH rigid fixes, via either thermostatting and/or barostatting, is
included in the [thermodynamic output](thermo_style) keywords *ecouple*
and *econserve*. See the [thermo_style](thermo_style) doc page for
details.

The 2 NVE rigid fixes compute a global scalar which can be accessed by
various [output commands](Howto_output). The scalar value calculated by
these fixes is \"intensive\". The scalar is the current temperature of
the collection of rigid bodies. This is averaged over all rigid bodies
and their translational and rotational degrees of freedom. The
translational energy of a rigid body is 1/2 m v\^2, where m = total mass
of the body and v = the velocity of its center of mass. The rotational
energy of a rigid body is 1/2 I w\^2, where I = the moment of inertia
tensor of the body and w = its angular velocity. Degrees of freedom
constrained by the *force* and *torque* keywords are removed from this
calculation, but only for the *rigid* and *rigid/nve* fixes.

The 6 NVT, NPT, NPH rigid fixes compute a global scalar which can be
accessed by various [output commands](Howto_output). The scalar is the
same cumulative energy change due to these fixes described above. The
scalar value calculated by this fix is \"extensive\".

The [fix_modify](fix_modify) *virial* option is supported by these fixes
to add the contribution due to the added forces on atoms to both the
global pressure and per-atom stress of the system via the [compute
pressure](compute_pressure) and [compute
stress/atom](compute_stress_atom) commands. The former can be accessed
by [thermodynamic output](thermo_style). The default setting for this
fix is [fix_modify virial yes](fix_modify).

All of the *rigid* styles (not the *rigid/small* styles) compute a
global array of values which can be accessed by various [output
commands](Howto_output). Similar information about the bodies defined by
the *rigid/small* styles can be accessed via the [compute
rigid/local](compute_rigid_local) command.

The number of rows in the array is equal to the number of rigid bodies.
The number of columns is 15. Thus for each rigid body, 15 values are
stored: the xyz coords of the center of mass (COM), the xyz components
of the COM velocity, the xyz components of the force acting on the COM,
the xyz components of the torque acting on the COM, and the xyz image
flags of the COM.

The center of mass (COM) for each body is similar to unwrapped
coordinates written to a dump file. It will always be inside (or
slightly outside) the simulation box. The image flags have the same
meaning as image flags for atom positions (see the \"dump\" command).
This means you can calculate the unwrapped COM by applying the image
flags to the COM, the same as when unwrapped coordinates are written to
a dump file.

The force and torque values in the array are not affected by the *force*
and *torque* keywords in the fix rigid command; they reflect values
before any changes are made by those keywords.

The ordering of the rigid bodies (by row in the array) is as follows.
For the *single* keyword there is just one rigid body. For the
*molecule* keyword, the bodies are ordered by ascending molecule ID. For
the *group* keyword, the list of group IDs determines the ordering of
bodies.

The array values calculated by these fixes are \"intensive\", meaning
they are independent of the number of atoms in the simulation.

No parameter of these fixes can be used with the *start/stop* keywords
of the [run](run) command. These fixes are not invoked during [energy
minimization](minimize).

------------------------------------------------------------------------

## Restrictions

These fixes are all part of the RIGID package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

Assigning a temperature via the [velocity create](velocity) command to a
system with [rigid bodies](fix_rigid) may not have the desired outcome
for two reasons. First, the velocity command can be invoked before the
rigid-body fix is invoked or initialized and the number of adjusted
degrees of freedom (DOFs) is known. Thus it is not possible to compute
the target temperature correctly. Second, the assigned velocities may be
partially canceled when constraints are first enforced, leading to a
different temperature than desired. A workaround for this is to perform
a [run 0](run) command, which ensures all DOFs are accounted for
properly, and then rescale the temperature to the desired value before
performing a simulation. For example:

``` LAMMPS
velocity all create 300.0 12345
run 0                             # temperature may not be 300K
velocity all scale 300.0          # now it should be
```

## Related commands

[delete_bonds](delete_bonds), [neigh_modify](neigh_modify) exclude, [fix
shake](fix_shake)

## Default

The option defaults are force \* on on on and torque \* on on on,
meaning all rigid bodies are acted on by center-of-mass force and
torque. Also Tchain = Pchain = 10, Titer = 1, Torder = 3, reinit = yes.

------------------------------------------------------------------------

::: {#Hoover}
**(Hoover)** Hoover, Phys Rev A, 31, 1695 (1985).
:::

::: {#Kamberaj}
**(Kamberaj)** Kamberaj, Low, Neal, J Chem Phys, 122, 224114 (2005).
:::

::: {#Martyna2}
**(Martyna)** Martyna, Klein, Tuckerman, J Chem Phys, 97, 2635 (1992);
Martyna, Tuckerman, Tobias, Klein, Mol Phys, 87, 1117.
:::

::: {#Miller3}
**(Miller)** Miller, Eleftheriou, Pattnaik, Ndirango, and Newns, J Chem
Phys, 116, 8649 (2002).
:::

::: {#Zhang1}
**(Zhang)** Zhang, Glotzer, Nanoletters, 4, 1407-1413 (2004).
:::
