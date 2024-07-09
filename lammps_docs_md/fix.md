# fix command

## Syntax

``` LAMMPS
fix ID group-ID style args
```

-   ID = user-assigned name for the fix
-   group-ID = ID of the group of atoms to apply the fix to
-   style = one of a long list of possible style names (see below)
-   args = arguments used by a particular style

## Examples

``` LAMMPS
fix 1 all nve
fix 3 all nvt temp 300.0 300.0 0.01
fix mine top setforce 0.0 NULL 0.0
```

## Description

Set a fix that will be applied to a group of atoms. In LAMMPS, a \"fix\"
is any operation that is applied to the system during timestepping or
minimization. Examples include updating of atom positions and velocities
due to time integration, controlling temperature, applying constraint
forces to atoms, enforcing boundary conditions, computing diagnostics,
etc. There are hundreds of fixes defined in LAMMPS and new ones can be
added; see the [Modify](Modify) page for details.

Fixes perform their operations at different stages of the timestep. If
two or more fixes operate at the same stage of the timestep, they are
invoked in the order they were specified in the input script.

The ID of a fix can only contain alphanumeric characters and
underscores.

Fixes can be deleted with the [unfix](unfix) command.

:::: note
::: title
Note
:::

The [unfix](unfix) command is the only way to turn off a fix; simply
specifying a new fix with a similar style will not turn off the first
one. This is especially important to realize for integration fixes. For
example, using a [fix nve](fix_nve) command for a second run after using
a [fix nvt](fix_nh) command for the first run will not cancel out the
NVT time integration invoked by the \"fix nvt\" command. Thus, two time
integrators would be in place!
::::

If you specify a new fix with the same ID and style as an existing fix,
the old fix is deleted and the new one is created (presumably with new
settings). This is the same as if an \"unfix\" command were first
performed on the old fix, except that the new fix is kept in the same
order relative to the existing fixes as the old one originally was. Note
that this operation also wipes out any additional changes made to the
old fix via the [fix_modify](fix_modify) command.

The [fix modify](fix_modify) command allows settings for some fixes to
be reset. See the page for individual fixes for details.

Some fixes store an internal \"state\" which is written to binary
restart files via the [restart](restart) or
[write_restart](write_restart) commands. This allows the fix to continue
on with its calculations in a restarted simulation. See the
[read_restart](read_restart) command for info on how to re-specify a fix
in an input script that reads a restart file. See the doc pages for
individual fixes for info on which ones can be restarted.

------------------------------------------------------------------------

Some fixes calculate one or more of four styles of quantities: global,
per-atom, local, or per-grid, which can be used by other commands or
output as described below. A global quantity is one or more system-wide
values, e.g. the energy of a wall interacting with particles. A per-atom
quantity is one or more values per atom, e.g. the displacement vector
for each atom since time 0. Per-atom values are set to 0.0 for atoms not
in the specified fix group. Local quantities are calculated by each
processor based on the atoms it owns, but there may be zero or more per
atoms. Per-grid quantities are calculated on a regular 2d or 3d grid
which overlays a 2d or 3d simulation domain. The grid points and the
data they store are distributed across processors; each processor owns
the grid points which fall within its subdomain.

Note that a single fix typically produces either global or per-atom or
local or per-grid values (or none at all). It does not produce both
global and per-atom. It can produce local or per-grid values in tandem
with global or per-atom values. The fix doc page will explain the
details.

Global, per-atom, local, and per-grid quantities come in three kinds: a
single scalar value, a vector of values, or a 2d array of values. The
doc page for each fix describes the style and kind of values it
produces, e.g. a per-atom vector. Some fixes produce more than one kind
of a single style, e.g. a global scalar and a global vector.

When a fix quantity is accessed, as in many of the output commands
discussed below, it can be referenced via the following bracket
notation, where ID is the ID of the fix:

  ---------------- --------------------------------------------
  f_ID             entire scalar, vector, or array

  f_ID\[I\]        one element of vector, one column of array

  f_ID\[I\]\[J\]   one element of array
  ---------------- --------------------------------------------

In other words, using one bracket reduces the dimension of the quantity
once (vector $\to$ scalar, array $\to$ vector). Using two brackets
reduces the dimension twice (array $\to$ scalar). Thus, a command that
uses scalar fix values as input can also process elements of a vector or
array.

Note that commands and [variables](variable) that use fix quantities
typically do not allow for all kinds (e.g., a command may require a
vector of values, not a scalar), and even if they do, the context in
which they are called can be used to resolve which output is being
requested. This means there is no ambiguity about referring to a fix
quantity as f_ID even if it produces, for example, both a scalar and
vector. The doc pages for various commands explain the details.

------------------------------------------------------------------------

In LAMMPS, the values generated by a fix can be used in several ways:

-   Global values can be output via the [thermo_style
    custom](thermo_style) or [fix ave/time](fix_ave_time) command.
    Alternatively, the values can be referenced in an [equal-style
    variable](variable) command.
-   Per-atom values can be output via the [dump custom](dump) command,
    or they can be time-averaged via the [fix ave/atom](fix_ave_atom)
    command or reduced by the [compute reduce](compute_reduce) command.
    Alternatively, per-atom values can be referenced in an [atom-style
    variable](variable).
-   Local values can be reduced by the [compute reduce](compute_reduce)
    command or histogrammed by the [fix ave/histo](fix_ave_histo)
    command. They can also be output by the [dump local](dump) command.

See the [Howto output](Howto_output) page for a summary of various
LAMMPS output options, many of which involve fixes.

The results of fixes that calculate global quantities can be either
\"intensive\" or \"extensive\" values. Intensive means the value is
independent of the number of atoms in the simulation (e.g.,
temperature). Extensive means the value scales with the number of atoms
in the simulation (e.g., total rotational kinetic energy).
[Thermodynamic output](thermo_style) will normalize extensive values by
the number of atoms in the system, depending on the \"thermo_modify
norm\" setting. It will not normalize intensive values. If a fix value
is accessed in another way (e.g., by a [variable](variable)), you may
want to know whether it is an intensive or extensive value. See the page
for individual fix styles for further info.

------------------------------------------------------------------------

Each fix style has its own page that describes its arguments and what it
does, as listed below. Here is an alphabetical list of fix styles
available in LAMMPS. They are also listed in more compact form on the
[Commands fix](Commands_fix) doc page.

There are also additional accelerated fix styles included in the LAMMPS
distribution for faster performance on CPUs, GPUs, and KNLs. The
individual style names on the [Commands fix](Commands_fix) doc page are
followed by one or more of (g,i,k,o,t) to indicate which accelerated
styles exist.

-   [accelerate/cos](fix_accelerate_cos) - apply cosine-shaped
    acceleration to atoms
-   [acks2/reaxff](fix_acks2_reaxff) - apply ACKS2 charge equilibration
-   [adapt](fix_adapt) - change a simulation parameter over time
-   [adapt/fep](fix_adapt_fep) - enhanced version of fix adapt
-   [addforce](fix_addforce) - add a force to each atom
-   [addtorque](fix_addtorque) - add a torque to a group of atoms
-   [alchemy](fix_alchemy) - perform an \"alchemical transformation\"
    between two partitions
-   [amoeba/bitorsion](fix_amoeba_bitorsion) - torsion/torsion terms in
    AMOEBA force field
-   [amoeba/pitorsion](fix_amoeba_pitorsion) - 6-body terms in AMOEBA
    force field
-   [append/atoms](fix_append_atoms) - append atoms to a running
    simulation
-   [atc](fix_atc) - initiates a coupled MD/FE simulation
-   [atom/swap](fix_atom_swap) - Monte Carlo atom type swapping
-   [ave/atom](fix_ave_atom) - compute per-atom time-averaged quantities
-   [ave/chunk](fix_ave_chunk) - compute per-chunk time-averaged
    quantities
-   [ave/correlate](fix_ave_correlate) - compute/output time
    correlations
-   [ave/correlate/long](fix_ave_correlate_long) - alternative to
    [ave/correlate](fix_ave_correlate) that allows efficient calculation
    over long time windows
-   [ave/grid](fix_ave_grid) - compute per-grid time-averaged quantities
-   [ave/histo](fix_ave_histo) - compute/output time-averaged histograms
-   [ave/histo/weight](fix_ave_histo) - weighted version of fix
    ave/histo
-   [ave/time](fix_ave_time) - compute/output global time-averaged
    quantities
-   [aveforce](fix_aveforce) - add an averaged force to each atom
-   [balance](fix_balance) - perform dynamic load-balancing
-   [brownian](fix_brownian) - overdamped translational brownian motion
-   [brownian/asphere](fix_brownian) - overdamped translational and
    rotational brownian motion for ellipsoids
-   [brownian/sphere](fix_brownian) - overdamped translational and
    rotational brownian motion for spheres
-   [bocs](fix_bocs) - NPT style time integration with pressure
    correction
-   [bond/break](fix_bond_break) - break bonds on the fly
-   [bond/create](fix_bond_create) - create bonds on the fly
-   [bond/create/angle](fix_bond_create) - create bonds on the fly with
    angle constraints
-   [bond/react](fix_bond_react) - apply topology changes to model
    reactions
-   [bond/swap](fix_bond_swap) - Monte Carlo bond swapping
-   [box/relax](fix_box_relax) - relax box size during energy
    minimization
-   [charge/regulation](fix_charge_regulation) - Monte Carlo sampling of
    charge regulation
-   [cmap](fix_cmap) - CMAP torsion/torsion terms in CHARMM force field
-   [colvars](fix_colvars) - interface to the collective variables
    \"Colvars\" library
-   [controller](fix_controller) - apply control loop feedback mechanism
-   [damping/cundall](fix_damping_cundall) - Cundall non-viscous damping
    for granular simulations
-   [deform](fix_deform) - change the simulation box size/shape
-   [deposit](fix_deposit) - add new atoms above a surface
-   [dpd/energy](fix_dpd_energy) - constant energy dissipative particle
    dynamics
-   [drag](fix_drag) - drag atoms towards a defined coordinate
-   [drude](fix_drude) - part of Drude oscillator polarization model
-   [drude/transform/direct](fix_drude_transform) - part of Drude
    oscillator polarization model
-   [drude/transform/inverse](fix_drude_transform) - part of Drude
    oscillator polarization model
-   [dt/reset](fix_dt_reset) - reset the timestep based on velocity,
    forces
-   [edpd/source](fix_dpd_source) - add heat source to eDPD simulations
-   [efield](fix_efield) - impose electric field on system
-   [efield/tip4p](fix_efield) - impose electric field on system with
    TIP4P molecules
-   [ehex](fix_ehex) - enhanced heat exchange algorithm
-   [electrode/conp](fix_electrode) - impose electric potential
-   [electrode/conq](fix_electrode) - impose total electric charge
-   [electrode/thermo](fix_electrode) - apply thermo-potentiostat
-   [electron/stopping](fix_electron_stopping) - electronic stopping
    power as a friction force
-   [electron/stopping/fit](fix_electron_stopping) - electronic stopping
    power as a friction force
-   [enforce2d](fix_enforce2d) - zero out *z*-dimension velocity and
    force
-   [eos/cv](fix_eos_cv) - applies a mesoparticle equation of state to
    relate the particle internal energy to the particle internal
    temperature
-   [eos/table](fix_eos_table) - applies a tabulated mesoparticle
    equation of state to relate the particle internal energy to the
    particle internal temperature
-   [eos/table/rx](fix_eos_table_rx) - applies a tabulated mesoparticle
    equation of state to relate the concentration-dependent particle
    internal energy to the particle internal temperature
-   [evaporate](fix_evaporate) - remove atoms from simulation
    periodically
-   [external](fix_external) - callback to an external driver program
-   [ffl](fix_ffl) - apply a Fast-Forward Langevin equation thermostat
-   [filter/corotate](fix_filter_corotate) - implement corotation filter
    to allow larger timesteps with r-RESPA
-   [flow/gauss](fix_flow_gauss) - Gaussian dynamics for constant mass
    flux
-   [freeze](fix_freeze) - freeze atoms in a granular simulation
-   [gcmc](fix_gcmc) - grand canonical insertions/deletions
-   [gld](fix_gld) - generalized Langevin dynamics integrator
-   [gle](fix_gle) - generalized Langevin equation thermostat
-   [gravity](fix_gravity) - add gravity to atoms in a granular
    simulation
-   [grem](fix_grem) - implements the generalized replica exchange
    method
-   [halt](fix_halt) - terminate a dynamics run or minimization
-   [heat](fix_heat) - add/subtract momentum-conserving heat
-   [heat/flow](fix_heat_flow) - plain time integration of heat flow
    with per-atom temperature updates
-   [hyper/global](fix_hyper_global) - global hyperdynamics
-   [hyper/local](fix_hyper_local) - local hyperdynamics
-   [imd](fix_imd) - implements the \"Interactive MD\" (IMD) protocol
-   [indent](fix_indent) - impose force due to an indenter
-   [ipi](fix_ipi) - enable LAMMPS to run as a client for i-PI
    path-integral simulations
-   [langevin](fix_langevin) - Langevin temperature control
-   [langevin/drude](fix_langevin_drude) - Langevin temperature control
    of Drude oscillators
-   [langevin/eff](fix_langevin_eff) - Langevin temperature control for
    the electron force field model
-   [langevin/spin](fix_langevin_spin) - Langevin temperature control
    for a spin or spin-lattice system
-   [lb/fluid](fix_lb_fluid) - lattice-Boltzmann fluid on a uniform mesh
-   [lb/momentum](fix_lb_momentum) - [fix momentum](fix_momentum)
    replacement for use with a lattice-Boltzmann fluid
-   [lb/viscous](fix_lb_viscous) - [fix viscous](fix_viscous)
    replacement for use with a lattice-Boltzmann fluid
-   [lineforce](fix_lineforce) - constrain atoms to move in a line
-   [manifoldforce](fix_manifoldforce) - restrain atoms to a manifold
    during minimization
-   [mdi/qm](fix_mdi_qm) - LAMMPS operates as a client for a quantum
    code via the MolSSI Driver Interface (MDI)
-   [mdi/qmmm](fix_mdi_qmmm) - LAMMPS operates as client for QM/MM
    simulation with a quantum code via the MolSSI Driver Interface (MDI)
-   [meso/move](fix_meso_move) - move mesoscopic SPH/SDPD particles in a
    prescribed fashion
-   [mol/swap](fix_mol_swap) - Monte Carlo atom type swapping with a
    molecule
-   [momentum](fix_momentum) - zero the linear and/or angular momentum
    of a group of atoms
-   [momentum/chunk](fix_momentum) - zero the linear and/or angular
    momentum of a chunk of atoms
-   [move](fix_move) - move atoms in a prescribed fashion
-   [mscg](fix_mscg) - apply MSCG method for force-matching to generate
    coarse grain models
-   [msst](fix_msst) - multi-scale shock technique (MSST) integration
-   [mvv/dpd](fix_mvv_dpd) - DPD using the modified velocity-Verlet
    integration algorithm
-   [mvv/edpd](fix_mvv_dpd) - constant energy DPD using the modified
    velocity-Verlet algorithm
-   [mvv/tdpd](fix_mvv_dpd) - constant temperature DPD using the
    modified velocity-Verlet algorithm
-   [neb](fix_neb) - nudged elastic band (NEB) spring forces
-   [neb/spin](fix_neb_spin) - nudged elastic band (NEB) spring forces
    for spins
-   [nph](fix_nh) - constant NPH time integration via Nose/Hoover
-   [nph/asphere](fix_nph_asphere) - NPH for aspherical particles
-   [nph/body](fix_nph_body) - NPH for body particles
-   [nph/eff](fix_nh_eff) - NPH for nuclei and electrons in the electron
    force field model
-   [nph/sphere](fix_nph_sphere) - NPH for spherical particles
-   [nphug](fix_nphug) - constant-stress Hugoniostat integration
-   [npt](fix_nh) - constant NPT time integration via Nose/Hoover
-   [npt/asphere](fix_npt_asphere) - NPT for aspherical particles
-   [npt/body](fix_npt_body) - NPT for body particles
-   [npt/cauchy](fix_npt_cauchy) - NPT with Cauchy stress
-   [npt/eff](fix_nh_eff) - NPT for nuclei and electrons in the electron
    force field model
-   [npt/sphere](fix_npt_sphere) - NPT for spherical particles
-   [npt/uef](fix_nh_uef) - NPT style time integration with diagonal
    flow
-   [numdiff](fix_numdiff) - numerically approximate atomic forces using
    finite energy differences
-   [numdiff/virial](fix_numdiff_virial) - numerically approximate
    virial stress tensor using finite energy differences
-   [nve](fix_nve) - constant NVE time integration
-   [nve/asphere](fix_nve_asphere) - NVE for aspherical particles
-   [nve/asphere/noforce](fix_nve_asphere_noforce) - NVE for aspherical
    particles without forces
-   [nve/awpmd](fix_nve_awpmd) - NVE for the Antisymmetrized Wave Packet
    Molecular Dynamics model
-   [nve/body](fix_nve_body) - NVE for body particles
-   [nve/dot](fix_nve_dot) - rigid body constant energy time integrator
    for coarse grain models
-   [nve/dotc/langevin](fix_nve_dotc_langevin) - Langevin style rigid
    body time integrator for coarse grain models
-   [nve/eff](fix_nve_eff) - NVE for nuclei and electrons in the
    electron force field model
-   [nve/limit](fix_nve_limit) - NVE with limited step length
-   [nve/line](fix_nve_line) - NVE for line segments
-   [nve/manifold/rattle](fix_nve_manifold_rattle) - NVE time
    integration for atoms constrained to a curved surface (manifold)
-   [nve/noforce](fix_nve_noforce) - NVE without forces (update
    positions only)
-   [nve/sphere](fix_nve_sphere) - NVE for spherical particles
-   [nve/bpm/sphere](fix_nve_bpm_sphere) - NVE for spherical particles
    used in the BPM package
-   [nve/spin](fix_nve_spin) - NVE for a spin or spin-lattice system
-   [nve/tri](fix_nve_tri) - NVE for triangles
-   [nvk](fix_nvk) - constant kinetic energy time integration
-   [nvt](fix_nh) - NVT time integration via Nose/Hoover
-   [nvt/asphere](fix_nvt_asphere) - NVT for aspherical particles
-   [nvt/body](fix_nvt_body) - NVT for body particles
-   [nvt/eff](fix_nh_eff) - NVE for nuclei and electrons in the electron
    force field model
-   [nvt/manifold/rattle](fix_nvt_manifold_rattle) - NVT time
    integration for atoms constrained to a curved surface (manifold)
-   [nvt/sllod](fix_nvt_sllod) - NVT for NEMD with SLLOD equations
-   [nvt/sllod/eff](fix_nvt_sllod_eff) - NVT for NEMD with SLLOD
    equations for the electron force field model
-   [nvt/sphere](fix_nvt_sphere) - NVT for spherical particles
-   [nvt/uef](fix_nh_uef) - NVT style time integration with diagonal
    flow
-   [oneway](fix_oneway) - constrain particles on move in one direction
-   [orient/bcc](fix_orient) - add grain boundary migration force for
    BCC
-   [orient/fcc](fix_orient) - add grain boundary migration force for
    FCC
-   [orient/eco](fix_orient_eco) - add generalized grain boundary
    migration force
-   [pafi](fix_pafi) - constrained force averages on hyper-planes to
    compute free energies (PAFI)
-   [pair](fix_pair) - access per-atom info from pair styles
-   [phonon](fix_phonon) - calculate dynamical matrix from MD
    simulations
-   [pimd/langevin](fix_pimd) - Feynman path-integral molecular dynamics
    with stochastic thermostat
-   [pimd/nvt](fix_pimd) - Feynman path-integral molecular dynamics with
    Nose-Hoover thermostat
-   [planeforce](fix_planeforce) - constrain atoms to move in a plane
-   [plumed](fix_plumed) - wrapper on PLUMED free energy library
-   [poems](fix_poems) - constrain clusters of atoms to move as coupled
    rigid bodies
-   [polarize/bem/gmres](fix_polarize) - compute induced charges at the
    interface between impermeable media with different dielectric
    constants with generalized minimum residual (GMRES)
-   [polarize/bem/icc](fix_polarize) - compute induced charges at the
    interface between impermeable media with different dielectric
    constants with the successive over-relaxation algorithm
-   [polarize/functional](fix_polarize) - compute induced charges at the
    interface between impermeable media with different dielectric
    constants with the energy variational approach
-   [pour](fix_pour) - pour new atoms/molecules into a granular
    simulation domain
-   [precession/spin](fix_precession_spin) - apply a precession torque
    to each magnetic spin
-   [press/berendsen](fix_press_berendsen) - pressure control by
    Berendsen barostat
-   [print](fix_print) - print text and variables during a simulation
-   [propel/self](fix_propel_self) - model self-propelled particles
-   [property/atom](fix_property_atom) - add customized per-atom values
-   [python/invoke](fix_python_invoke) - call a Python function during a
    simulation
-   [python/move](fix_python_move) - move particles using a Python
    function during a simulation run
-   [qbmsst](fix_qbmsst) - quantum bath multi-scale shock technique time
    integrator
-   [qeq/comb](fix_qeq_comb) - charge equilibration for COMB potential
-   [qeq/dynamic](fix_qeq) - charge equilibration via dynamic method
-   [qeq/fire](fix_qeq) - charge equilibration via FIRE minimizer
-   [qeq/point](fix_qeq) - charge equilibration via point method
-   [qeq/reaxff](fix_qeq_reaxff) - charge equilibration for ReaxFF
    potential
-   [qeq/shielded](fix_qeq) - charge equilibration via shielded method
-   [qeq/slater](fix_qeq) - charge equilibration via Slater method
-   [qmmm](fix_qmmm) - functionality to enable a quantum
    mechanics/molecular mechanics coupling
-   [qtb](fix_qtb) - implement quantum thermal bath scheme
-   [rattle](fix_shake) - RATTLE constraints on bonds and/or angles
-   [reaxff/bonds](fix_reaxff_bonds) - write out ReaxFF bond information
-   [reaxff/species](fix_reaxff_species) - write out ReaxFF molecule
    information
-   [recenter](fix_recenter) - constrain the center-of-mass position of
    a group of atoms
-   [restrain](fix_restrain) - constrain a bond, angle, dihedral
-   [rhok](fix_rhok) - add bias potential for long-range ordered systems
-   [rigid](fix_rigid) - constrain one or more clusters of atoms to move
    as a rigid body with NVE integration
-   [rigid/meso](fix_rigid_meso) - constrain clusters of mesoscopic
    SPH/SDPD particles to move as a rigid body
-   [rigid/nph](fix_rigid) - constrain one or more clusters of atoms to
    move as a rigid body with NPH integration
-   [rigid/nph/small](fix_rigid) - constrain many small clusters of
    atoms to move as a rigid body with NPH integration
-   [rigid/npt](fix_rigid) - constrain one or more clusters of atoms to
    move as a rigid body with NPT integration
-   [rigid/npt/small](fix_rigid) - constrain many small clusters of
    atoms to move as a rigid body with NPT integration
-   [rigid/nve](fix_rigid) - constrain one or more clusters of atoms to
    move as a rigid body with alternate NVE integration
-   [rigid/nve/small](fix_rigid) - constrain many small clusters of
    atoms to move as a rigid body with alternate NVE integration
-   [rigid/nvt](fix_rigid) - constrain one or more clusters of atoms to
    move as a rigid body with NVT integration
-   [rigid/nvt/small](fix_rigid) - constrain many small clusters of
    atoms to move as a rigid body with NVT integration
-   [rigid/small](fix_rigid) - constrain many small clusters of atoms to
    move as a rigid body with NVE integration
-   [rx](fix_rx) - solve reaction kinetic ODEs for a defined reaction
    set
-   [saed/vtk](fix_saed_vtk) - time-average the intensities from
    [compute saed](compute_saed)
-   [setforce](fix_setforce) - set the force on each atom
-   [setforce/spin](fix_setforce) - set magnetic precession vectors on
    each atom
-   [sgcmc](fix_sgcmc) - fix for hybrid semi-grand canonical MD/MC
    simulations
-   [shake](fix_shake) - SHAKE constraints on bonds and/or angles
-   [shardlow](fix_shardlow) - integration of DPD equations of motion
    using the Shardlow splitting
-   [smd](fix_smd) - applied a steered MD force to a group
-   [smd/adjust_dt](fix_smd_adjust_dt) - calculate a new stable time
    increment for use with SMD integrators
-   [smd/integrate_tlsph](fix_smd_integrate_tlsph) - explicit time
    integration with total Lagrangian SPH pair style
-   [smd/integrate_ulsph](fix_smd_integrate_ulsph) - explicit time
    integration with updated Lagrangian SPH pair style
-   [smd/move_tri_surf](fix_smd_move_triangulated_surface) - update
    position and velocity near rigid surfaces using SPH integrators
-   [smd/setvel](fix_smd_setvel) - sets each velocity component,
    ignoring forces, for Smooth Mach Dynamics
-   [smd/wall_surface](fix_smd_wall_surface) - create a rigid wall with
    a triangulated surface for use in Smooth Mach Dynamics
-   [sph](fix_sph) - time integration for SPH/DPDE particles
-   [sph/stationary](fix_sph_stationary) - update energy and density but
    not position or velocity in Smooth Particle Hydrodynamics
-   [spring](fix_spring) - apply harmonic spring force to group of atoms
-   [spring/chunk](fix_spring_chunk) - apply harmonic spring force to
    each chunk of atoms
-   [spring/rg](fix_spring_rg) - spring on radius of gyration of group
    of atoms
-   [spring/self](fix_spring_self) - spring from each atom to its origin
-   [srd](fix_srd) - stochastic rotation dynamics (SRD)
-   [store/force](fix_store_force) - store force on each atom
-   [store/state](fix_store_state) - store attributes for each atom
-   [tdpd/source](fix_dpd_source) - add external concentration source
-   [temp/berendsen](fix_temp_berendsen) - temperature control by
    Berendsen thermostat
-   [temp/csld](fix_temp_csvr) - canonical sampling thermostat with
    Langevin dynamics
-   [temp/csvr](fix_temp_csvr) - canonical sampling thermostat with
    Hamiltonian dynamics
-   [temp/rescale](fix_temp_rescale) - temperature control by velocity
    rescaling
-   [temp/rescale/eff](fix_temp_rescale_eff) - temperature control by
    velocity rescaling in the electron force field model
-   [tfmc](fix_tfmc) - perform force-bias Monte Carlo with time-stamped
    method
-   [tgnvt/drude](fix_tgnh_drude) - NVT time integration for Drude
    polarizable model via temperature-grouped Nose-Hoover
-   [tgnpt/drude](fix_tgnh_drude) - NPT time integration for Drude
    polarizable model via temperature-grouped Nose-Hoover
-   [thermal/conductivity](fix_thermal_conductivity) - Mueller-Plathe
    kinetic energy exchange for thermal conductivity calculation
-   [ti/spring](fix_ti_spring) - perform thermodynamic integration
    between a solid and an Einstein crystal
-   [tmd](fix_tmd) - guide a group of atoms to a new configuration
-   [ttm](fix_ttm) - two-temperature model for electronic/atomic
    coupling (replicated grid)
-   [ttm/grid](fix_ttm) - two-temperature model for electronic/atomic
    coupling (distributed grid)
-   [ttm/mod](fix_ttm) - enhanced two-temperature model with additional
    options
-   [tune/kspace](fix_tune_kspace) - auto-tune $k$-space parameters
-   [vector](fix_vector) - accumulate a global vector every *N*
    timesteps
-   [viscosity](fix_viscosity) - Mueller-Plathe momentum exchange for
    viscosity calculation
-   [viscous](fix_viscous) - viscous damping for granular simulations
-   [viscous/sphere](fix_viscous_sphere) - viscous damping on angular
    velocity for granular simulations
-   [wall/body/polygon](fix_wall_body_polygon) - time integration for
    body particles of style [rounded/polygon](Howto_body)
-   [wall/body/polyhedron](fix_wall_body_polyhedron) - time integration
    for body particles of style [rounded/polyhedron](Howto_body)
-   [wall/colloid](fix_wall) - Lennard-Jones wall interacting with
    finite-size particles
-   [wall/ees](fix_wall_ees) - wall for ellipsoidal particles
-   [wall/gran](fix_wall_gran) - frictional wall(s) for granular
    simulations
-   [wall/gran/region](fix_wall_gran_region) - [fix
    wall/region](fix_wall_region) equivalent for use with granular
    particles
-   [wall/harmonic](fix_wall) - harmonic spring wall
-   [wall/lj1043](fix_wall) - Lennard-Jones 10\--4\--3 wall
-   [wall/lj126](fix_wall) - Lennard-Jones 12\--6 wall
-   [wall/lj93](fix_wall) - Lennard-Jones 9\--3 wall
-   [wall/lepton](fix_wall) - Custom Lepton expression wall
-   [wall/morse](fix_wall) - Morse potential wall
-   [wall/piston](fix_wall_piston) - moving reflective piston wall
-   [wall/reflect](fix_wall_reflect) - reflecting wall(s)
-   [wall/reflect/stochastic](fix_wall_reflect_stochastic) - reflecting
    wall(s) with finite temperature
-   [wall/region](fix_wall_region) - use region surface as wall
-   [wall/region/ees](fix_wall_ees) - use region surface as wall for
    ellipsoidal particles
-   [wall/srd](fix_wall_srd) - slip/no-slip wall for SRD particles
-   [wall/table](fix_wall) - Tabulated potential wall wall
-   [widom](fix_widom) - Widom insertions of atoms or molecules

## Restrictions

Some fix styles are part of specific packages. They are only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info. The doc pages for individual
fixes tell if it is part of a package.

## Related commands

[unfix](unfix), [fix_modify](fix_modify)

## Default

none
