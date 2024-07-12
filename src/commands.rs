/// A macro for constructing the CommandName enum.
/// Takes a list of the variant Identifiers, actual command names and the minimimum number of
/// arguments needed for the command.
macro_rules! commands {
    ($(($variant:tt, $lit:literal, $nargs:literal) ),+) => {

        #[allow(non_camel_case_types)]
        #[derive(Debug,Eq,PartialEq,Clone,Hash)]
        /// Enum of all possible valid LAMMPS commands.
        pub enum CommandName {
            $($variant,)+
            InvalidCommand(String),
        }

        impl CommandName {
            /// Returns the minimum number of positional arguments the compute style takes
            /// Note: This does not include the required arguments common to all computes, i.e.
            /// the `compute` keyword, the compute id, the group id and the compute style name
            pub const fn n_positional_args(&self) -> usize {
                match self {
                    $(CommandName::$variant => $nargs,)+
                    CommandName::InvalidCommand(_) => 0,
                }

            }
        }

        // TODO: Decide if try from or From is more appropriate
        impl From<&str> for CommandName {
            fn from(value: &str) -> Self {
                match value {
                    $($lit => CommandName::$variant,)+
                    s => CommandName::InvalidCommand(s.to_owned()),
                }
            }
        }
        impl std::fmt::Display for CommandName {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match self {
                    $(CommandName::$variant => write!(f, $lit),)+
                    CommandName::InvalidCommand(s) => write!(f, "Invalid Command: {s}"),
                }
            }
        }
    };
}

commands!(
    // Variant Name, Command Name, Min number of args.
    (AngleCoeff, "angle_coeff", 1), // Assumes some angle_style exists with zero arguments
    (AngleStyle, "angle_style", 1),
    (AngleWrite, "angle_write", 4),
    (AtomModify, "atom_modify", 2), // One or more kword value pairs.
    (AtomStyle, "atom_style", 1),
    (Balance, "balance", 2),
    (BondCoeff, "bond_coeff", 1),
    (BondStyle, "bond_style", 1),
    (BondWrite, "bond_write", 7),
    (Boundary, "boundary", 3),
    (ChangeBox, "change_box", 2), // Likely needs a lot more than two
    (Clear, "clear", 0),
    (CommModify, "comm_modify", 2),
    (CommStyle, "comm_style", 1),
    (Compute, "compute", 3), //  ID, Group-ID, style
    (ComputeModify, "compute_modify", 2),
    (CreateAtoms, "create_atoms", 2), // Requires type and style.
    (CreateBonds, "create_bonds", 1),
    (CreateBox, "create_box", 2),
    (DeleteAtoms, "delete_atoms", 1),
    (DeleteBonds, "delete_bonds", 2),
    (Dielectric, "dielectric", 1),
    (DihedralCoeff, "dihedral_coeff", 1),
    (DihedralStyle, "dihedral_style", 1),
    (DihedralWrite, "dihedral_write", 4),
    (Dimension, "dimension", 1),
    (DisplaceAtoms, "displace_atoms", 3),
    (Dump, "dump", 5),
    (DumpModify, "dump_modify", 2),
    (DynamicalMatrix, "dynamical_matrix", 3),
    (Echo, "echo", 1),
    (Fix, "fix", 3),
    (FixModify, "fix_modify", 3),
    (Fitpod, "fitpod", 2),
    (Group, "group", 2),
    (Group2Ndx, "group2ndx", 1),
    (Ndx2Group, "ndx2group", 1),
    (Hyper, "hyper", 4),
    (If, "if", 1), // Then and following command may be compulsory.
    (ImproperCoeff, "improper_coeff", 1),
    (ImproperStyle, "improper_style", 1),
    (Include, "include", 1), // label
    (Info, "info", 1),
    (Jump, "jump", 2), // file and label
    (Kim, "kim", 1),
    (KspaceModify, "kspace_modify", 2),
    (KspaceStyle, "kspace_style", 1),
    (Label, "label", 1),
    (Labelmap, "labelmap", 1),
    (Lattice, "lattice", 2), // Style and scale are mandatory
    (Log, "log", 1),         // File + optional append
    (Mass, "mass", 2),
    (Mdi, "mdi", 1),
    (MinModify, "min_modify", 2),
    (MinStyle, "min_style", 1),
    (Minimize, "minimize", 4), // etol ftol maxiter maxeval
    (Molecule, "molecule", 2), // ID file1 keyword values ... file2 keyword values ... fileN ...
    (Neb, "neb", 6),           // NOTE: May require more?
    (NebSpin, "neb/spin", 6),
    (NeighModify, "neigh_modify", 2),
    (Neighbor, "neighbor", 2),
    (Newton, "newton", 1),
    (Next, "next", 1),
    (Package, "package", 1),
    (PairCoeff, "pair_coeff", 2),
    (PairModify, "pair_modify", 2),
    (PairStyle, "pair_style", 1),
    (PairWrite, "pair_write", 8),
    (Partition, "partition", 3),
    (Plugin, "plugin", 1),
    (Prd, "prd", 7),
    (Print, "print", 2),
    (Processors, "processors", 3),
    (Python, "python", 1),
    (Quit, "quit", 0),
    (ReadData, "read_data", 1),
    (ReadDump, "read_dump", 2),
    (ReadRestart, "read_restart", 1),
    (Region, "region", 2),
    (Replicate, "replicate", 3),
    (Rerun, "rerun", 1),
    (ResetAtoms, "reset_atoms", 1), // Requires at least the property field.
    (ResetTimestep, "reset_timestep", 1),
    (Restart, "restart", 1),
    (Run, "run", 1),
    (RunStyle, "run_style", 1),
    (Set, "set", 2), // At least style and ID May need more.
    (Shell, "shell", 1),
    (SpecialBonds, "special_bonds", 1),
    (Suffix, "suffix", 1),
    (Tad, "tad", 7),
    (Temper, "temper", 6),
    (TemperGrem, "temper/grem", 7),
    (TemperNpt, "temper/npt", 7),
    (Thermo, "thermo", 1),
    (ThermoModify, "thermo_modify", 2),
    (ThermoStyle, "thermo_style", 1),
    (ThirdOrder, "third_order", 3),
    (Timer, "timer", 1),
    (Timestep, "timestep", 1),
    (Uncompute, "uncompute", 1),
    (Undump, "undump", 1),
    (Unfix, "unfix", 1),
    (Units, "units", 1),
    (Variable, "variable", 2),
    (Velocity, "velocity", 2),
    (WriteCoeff, "write_coeff", 1),
    (WriteData, "write_data", 1),
    (WriteDump, "write_dump", 3),
    (WriteRestart, "write_restart", 1)
);

impl Default for CommandName {
    fn default() -> Self {
        CommandName::InvalidCommand("".to_owned())
    }
}
