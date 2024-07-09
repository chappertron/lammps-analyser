# LAMMPS features

LAMMPS is a classical molecular dynamics (MD) code with these general
classes of functionality:

1.  [General features](general)
2.  [Particle and model types](particle)
3.  [Interatomic potentials (force fields)](ff)
4.  [Atom creation](create)
5.  [Ensembles, constraints, and boundary conditions](ensemble)
6.  [Integrators](integrate)
7.  [Diagnostics](diag)
8.  [Output](output)
9.  [Multi-replica models](replica1)
10. [Pre- and post-processing](prepost)
11. [Specialized features (beyond MD itself)](special)

------------------------------------------------------------------------

## General features {#general}

-   runs on a single processor or in parallel
-   distributed memory message-passing parallelism (MPI)
-   shared memory multi-threading parallelism (OpenMP)
-   spatial decomposition of simulation domain for MPI parallelism
-   particle decomposition inside spatial decomposition for OpenMP and
    GPU parallelism
-   GPLv2 licensed open-source distribution
-   highly portable C++-11
-   modular code with most functionality in optional packages
-   only depends on MPI library for basic parallel functionality, MPI
    stub for serial compilation
-   other libraries are optional and only required for specific packages
-   GPU (CUDA, OpenCL, HIP, SYCL), Intel Xeon Phi, and OpenMP support
    for many code features
-   easy to extend with new features and functionality
-   runs from an input script
-   syntax for defining and using variables and formulas
-   syntax for looping over runs and breaking out of loops
-   run one or multiple simulations simultaneously (in parallel) from
    one script
-   build as library, invoke LAMMPS through library interface (from C,
    C++, Fortran) or provided Python wrapper or SWIG based wrappers
-   couple with other codes: LAMMPS calls other code, other code calls
    LAMMPS, umbrella code calls both, MDI coupling interface
-   call out to Python for computing forces, time integration, or other
    tasks
-   plugin interface for loading external features at runtime
-   large integrated collection of tests

## Particle and model types {#particle}

(See [atom style](atom_style) command)

-   atoms
-   coarse-grained particles (e.g. bead-spring polymers)
-   united-atom polymers or organic molecules
-   all-atom polymers, organic molecules, proteins, DNA
-   metals
-   metal oxides
-   granular materials
-   coarse-grained mesoscale models
-   finite-size spherical and ellipsoidal particles
-   finite-size line segment (2d) and triangle (3d) particles
-   finite-size rounded polygons (2d) and polyhedra (3d) particles
-   point dipole particles
-   particles with magnetic spin
-   rigid collections of n particles
-   hybrid combinations of these

## Interatomic potentials (force fields) {#ff}

(See [pair style](pair_style), [bond style](bond_style), [angle
style](angle_style), [dihedral style](dihedral_style), [improper
style](improper_style), [kspace style](kspace_style) commands)

-   pairwise potentials: Lennard-Jones, Buckingham, Morse,
    Born-Mayer-Huggins, Yukawa, soft, Class II (COMPASS), hydrogen bond,
    harmonic, gaussian, tabulated, scripted
-   charged pairwise potentials: Coulombic, point-dipole
-   many-body potentials: EAM, Finnis/Sinclair, MEAM, MEAM+SW, EIM,
    EDIP, ADP, Stillinger-Weber, Tersoff, REBO, AIREBO, ReaxFF, COMB,
    Streitz-Mintmire, 3-body polymorphic, BOP, Vashishta
-   machine learning potentials: ACE, AGNI, GAP, Behler-Parrinello
    (N2P2), POD, RANN
-   interfaces to ML potentials distributed by external groups: ANI,
    ChIMES, DeepPot, HIPNN, MTP
-   long-range interactions for charge, point-dipoles, and LJ
    dispersion: Ewald, Wolf, PPPM (similar to particle-mesh Ewald), MSM,
    ScaFaCoS
-   polarization models: [QEq](fix_qeq), [core/shell
    model](Howto_coreshell), [Drude dipole model](Howto_drude)
-   charge equilibration (QEq via dynamic, point, shielded, Slater
    methods)
-   coarse-grained potentials: DPD, GayBerne, REsquared, colloidal,
    DLVO, oxDNA / oxRNA, SPICA
-   mesoscopic potentials: granular, Peridynamics, SPH, mesoscopic
    tubular potential (MESONT)
-   semi-empirical potentials: multi-ion generalized pseudopotential
    theory (MGPT), second moment tight binding + QEq (SMTB-Q)
-   electron force field (eFF, AWPMD)
-   bond potentials: harmonic, FENE, Morse, nonlinear, Class II
    (COMPASS), quartic (breakable), tabulated, scripted
-   angle potentials: harmonic, CHARMM, cosine, cosine/squared,
    cosine/periodic, Class II (COMPASS), tabulated, scripted
-   dihedral potentials: harmonic, CHARMM, multi-harmonic, helix, Class
    II (COMPASS), OPLS, tabulated, scripted
-   improper potentials: harmonic, cvff, umbrella, Class II (COMPASS),
    tabulated
-   polymer potentials: all-atom, united-atom, bead-spring, breakable
-   water potentials: TIP3P, TIP4P, SPC, SPC/E and variants
-   interlayer potentials for graphene and analogues, hetero-junctions
-   metal-organic framework potentials (QuickFF, MO-FF)
-   implicit solvent potentials: hydrodynamic lubrication, Debye
-   force-field compatibility with CHARMM, AMBER, DREIDING, OPLS,
    GROMACS, Class II (COMPASS), UFF, ClayFF, DREIDING, AMOEBA,
    INTERFACE
-   access to the [OpenKIM Repository](https://openkim.org)\_ of
    potentials via the [kim command](kim_commands)
-   hybrid potentials: multiple pair, bond, angle, dihedral, improper
    potentials can be used in one simulation
-   overlaid potentials: superposition of multiple pair potentials
    (including many-body) with optional scale factor

## Atom creation {#create}

(See [read_data](read_data), [lattice](lattice),
[create_atoms](create_atoms), [delete_atoms](delete_atoms),
[displace_atoms](displace_atoms), [replicate](replicate) commands)

-   read in atom coordinates from files
-   create atoms on one or more lattices (e.g. grain boundaries)
-   delete geometric or logical groups of atoms (e.g. voids)
-   replicate existing atoms multiple times
-   displace atoms

## Ensembles, constraints, and boundary conditions {#ensemble}

(See [fix](fix) command)

-   2d or 3d systems
-   orthogonal or non-orthogonal (triclinic symmetry) simulation domains
-   constant NVE, NVT, NPT, NPH, Parrinello/Rahman integrators
-   thermostatting options for groups and geometric regions of atoms
-   pressure control via Nose/Hoover or Berendsen barostatting in 1 to 3
    dimensions
-   simulation box deformation (tensile and shear)
-   harmonic (umbrella) constraint forces
-   rigid body constraints
-   SHAKE / RATTLE bond and angle constraints
-   motion constraints to manifold surfaces
-   Monte Carlo bond breaking, formation, swapping, template based
    reaction modeling
-   atom/molecule insertion and deletion
-   walls of various kinds, static and moving
-   non-equilibrium molecular dynamics (NEMD)
-   variety of additional boundary conditions and constraints

## Integrators {#integrate}

(See [run](run), [run_style](run_style), [minimize](minimize) commands)

-   velocity-Verlet integrator
-   Brownian dynamics
-   rigid body integration
-   energy minimization via conjugate gradient, steepest descent
    relaxation, or damped dynamics (FIRE, Quickmin)
-   rRESPA hierarchical timestepping
-   fixed or adaptive time step
-   rerun command for post-processing of dump files

## Diagnostics {#diag}

-   see various flavors of the [fix](fix) and [compute](compute)
    commands
-   introspection command for system, simulation, and compile time
    settings and configurations

## Output

([dump](dump), [restart](restart) commands)

-   log file of thermodynamic info
-   text dump files of atom coordinates, velocities, other per-atom
    quantities
-   dump output on fixed and variable intervals, based timestep or
    simulated time
-   binary restart files
-   parallel I/O of dump and restart files
-   per-atom quantities (energy, stress, centro-symmetry parameter, CNA,
    etc.)
-   user-defined system-wide (log file) or per-atom (dump file)
    calculations
-   custom partitioning (chunks) for binning, and static or dynamic
    grouping of atoms for analysis
-   spatial, time, and per-chunk averaging of per-atom quantities
-   time averaging and histogramming of system-wide quantities
-   atom snapshots in native, XYZ, XTC, DCD, CFG, NetCDF, HDF5, ADIOS2,
    YAML formats
-   on-the-fly compression of output and decompression of read in files

## Multi-replica models {#replica1}

-   [nudged elastic band](neb)
-   [hyperdynamics](hyper)
-   [parallel replica dynamics](prd)
-   [temperature accelerated dynamics](tad)
-   [parallel tempering](temper)
-   path-integral MD: [first variant](fix_pimd), [second
    variant](fix_ipi)
-   multi-walker collective variables with [Colvars](fix_colvars) and
    [Plumed](fix_plumed)

## Pre- and post-processing {#prepost}

-   A handful of pre- and post-processing tools are packaged with
    LAMMPS, some of which can convert input and output files to/from
    formats used by other codes; see the [Tools](Tools) page.
-   Our group has also written and released a separate toolkit called
    [Pizza.py](https://lammps.github.io/pizza)\_ which provides tools
    for doing setup, analysis, plotting, and visualization for LAMMPS
    simulations. Pizza.py is written in
    [Python](https://www.python.org)\_ and is available for download
    from [the Pizza.py WWW site](https://lammps.github.io/pizza)\_.

## Specialized features {#special}

LAMMPS can be built with optional packages which implement a variety of
additional capabilities. See the [Optional Packages](Packages) page for
details.

These are LAMMPS capabilities which you may not think of as typical
classical MD options:

-   [static](balance) and [dynamic load-balancing](fix_balance),
    optional with recursive bisectioning decomposition
-   [generalized aspherical particles](Howto_body)
-   [stochastic rotation dynamics (SRD)](fix_srd)
-   [real-time visualization and interactive MD](fix_imd), [built-in
    renderer for images and movies](dump_image)
-   calculate [virtual diffraction patterns](compute_xrd)
-   calculate [finite temperature phonon dispersion](fix_phonon) and the
    [dynamical matrix of minimized structures](dynamical_matrix)
-   [atom-to-continuum coupling](fix_atc) with finite elements
-   coupled rigid body integration via the [POEMS](fix_poems) library
-   [QM/MM coupling](fix_qmmm)
-   Monte Carlo via [GCMC](fix_gcmc) and [tfMC](fix_tfmc) and [atom
    swapping](fix_atom_swap)
-   [path-integral molecular dynamics (PIMD)](fix_ipi) and [this as
    well](fix_pimd)
-   [Direct Simulation Monte Carlo](pair_dsmc) for low-density fluids
-   [Peridynamics modeling](pair_peri)
-   [Lattice Boltzmann fluid](fix_lb_fluid)
-   [targeted](fix_tmd) and [steered](fix_smd) molecular dynamics
-   [two-temperature electron model](fix_ttm)
