# Input script structure

This page describes the structure of a typical LAMMPS input script. The
examples directory in the LAMMPS distribution contains many sample input
scripts; it is discussed on the [Examples](Examples) doc page.

A LAMMPS input script typically has 4 parts:

1.  [Initialization](init)
2.  [System definition](system)
3.  [Simulation settings](settings)
4.  [Run a simulation](run)

The last 2 parts can be repeated as many times as desired. I.e. run a
simulation, change some settings, run some more, etc. Each of the 4
parts is now described in more detail. Remember that almost all commands
need only be used if a non-default value is desired.

## Initialization {#init}

Set parameters that need to be defined before atoms are created or
read-in from a file.

The relevant commands are [units](units), [dimension](dimension),
[newton](newton), [processors](processors), [boundary](boundary),
[atom_style](atom_style), [atom_modify](atom_modify).

If force-field parameters appear in the files that will be read, these
commands tell LAMMPS what kinds of force fields are being used:
[pair_style](pair_style), [bond_style](bond_style),
[angle_style](angle_style), [dihedral_style](dihedral_style),
[improper_style](improper_style).

## System definition {#system}

There are 3 ways to define the simulation cell and reserve space for
force field info and fill it with atoms in LAMMPS. Read them in from (1)
a data file or (2) a restart file via the [read_data](read_data) or
[read_restart](read_restart) commands, respectively. These files can
also contain molecular topology information. Or (3) create a simulation
cell and fill it with atoms on a lattice (with no molecular topology),
using these commands: [lattice](lattice), [region](region),
[create_box](create_box), [create_atoms](create_atoms) or
[read_dump](read_dump).

The entire set of atoms can be duplicated to make a larger simulation
using the [replicate](replicate) command.

## Simulation settings {#settings}

Once atoms and molecular topology are defined, a variety of settings can
be specified: force field coefficients, simulation parameters, output
options, and more.

Force field coefficients are set by these commands (they can also be set
in the read-in files): [pair_coeff](pair_coeff),
[bond_coeff](bond_coeff), [angle_coeff](angle_coeff),
[dihedral_coeff](dihedral_coeff), [improper_coeff](improper_coeff),
[kspace_style](kspace_style), [dielectric](dielectric),
[special_bonds](special_bonds).

Various simulation parameters are set by these commands:
[neighbor](neighbor), [neigh_modify](neigh_modify), [group](group),
[timestep](timestep), [reset_timestep](reset_timestep),
[run_style](run_style), [min_style](min_style),
[min_modify](min_modify).

Fixes impose a variety of boundary conditions, time integration, and
diagnostic options. The [fix](fix) command comes in many flavors.

Various computations can be specified for execution during a simulation
using the [compute](compute), [compute_modify](compute_modify), and
[variable](variable) commands.

Output options are set by the [thermo](thermo), [dump](dump), and
[restart](restart) commands.

## Run a simulation {#run}

A molecular dynamics simulation is run using the [run](run) command.
Energy minimization (molecular statics) is performed using the
[minimize](minimize) command. A parallel tempering (replica-exchange)
simulation can be run using the [temper](temper) command.
