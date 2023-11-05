pub struct InvalidFix(String);

macro_rules! fix_styles {
    ($(($variant:tt, $lit:literal, $nargs:literal) ),+) => {

        #[allow(non_camel_case_types)]
        #[derive(Debug,Default,Eq,PartialEq,Clone,Copy)]
        pub enum FixStyle {
            $($variant,)+
            #[default] // Because no other varitant can be default with macro def
            InvalidFixStyle,
        }

        impl FixStyle {
            pub const fn n_positional_args(&self) -> usize {
                match self {
                    $(FixStyle::$variant => $nargs,)+
                    FixStyle::InvalidFixStyle => 0,
                }

            }
        }

        // TODO: Decide if try from or From is more appropriate

        // impl TryFrom<&str> for FixStyle {
        //    type Error = InvalidFix;
        //         fn try_from(value: &str) -> Result<Self, Self::Error> {
        //             match value {
        //                 $($lit => Ok(FixStyle::$variant),)+
        //                 value => Err(InvalidFix(value.into())),
        //             }
        //         }
        // }

        impl From<&str> for FixStyle {
            fn from(value: &str) -> Self {
                match value {
                    $($lit => FixStyle::$variant,)+
                    _ => FixStyle::InvalidFixStyle,
                }
            }
        }
        impl std::fmt::Display for FixStyle {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match self {
                    $(FixStyle::$variant => write!(f, $lit),)+
                    FixStyle::InvalidFixStyle => write!(f, "INVALID FIX STYLE"),
                }
            }
        }
    };
}

// Defines an enum of all the possible fix styles.
// (EnumVariant, fix_name, number_of_positional_args)
fix_styles!(
    (AccelerateCos, "accelerate/cos", 1),
    (Acks2Reaxff, "acks2/reaxff", 5),
    (Adapt, "adapt", 4), // All attributes have at least one arg
    (AdaptFep, "adapt/fep", 4),
    (Addforce, "addforce", 3),
    (Addtorque, "addtorque", 3),
    (Alchemy, "alchemy", 1),
    (AmoebaBitorsion, "amoeba/bitorsion", 1),
    (AmoebaPitorsion, "amoeba/pitorsion", 0),
    (AppendAtoms, "append/atoms", 1),
    (Atc, "atc", 2),
    (AtomSwap, "atom/swap", 4),
    (AveAtom, "ave/atom", 4),
    (AveChunk, "ave/chunk", 5),
    (AveCorrelate, "ave/correlate", 4),
    (AveCorrelateLong, "ave/correlate/long", 4),
    (AveGrid, "ave/grid", 7),
    (AveHisto, "ave/histo", 7),
    (AveHistoWeight, "ave/histo/weight", 7),
    // DEPRECATED
    (AveSpatial, "ave/spatial", 0),
    (AveTime, "ave/time", 4),
    (Aveforce, "aveforce", 3),
    (Balance, "balance", 3), // Some styles have more args
    (Bocs, "bocs", 0),
    (BondBreak, "bond/break", 3),
    (BondCreate, "bond/create", 6),
    (BondCreateAngle, "bond/create/angle", 6),
    (BondReact, "bond/react", 0),
    (BondSwap, "bond/swap", 4),
    (BoxRelax, "box/relax", 0),
    (Brownian, "brownian", 2),
    (BrownianSphere, "brownian/sphere", 2),
    (BrownianAsphere, "brownian/asphere", 2),
    (ChargeRegulation, "charge/regulation", 2),
    (Cmap, "cmap", 1),
    (Colvars, "colvars", 1),
    (Controller, "controller", 8),
    (DampingCundall, "damping/cundall", 2),
    (Deform, "deform", 1),
    (Deposit, "deposit", 4),
    (DpdEnergy, "dpd/energy", 0),
    (EdpdSource, "edpd/source", 0),
    (TdpdSource, "tdpd/source", 1),
    (Drag, "drag", 5),
    (Drude, "drude", 1), // Should have between 1 and N arguments,
    // N being the number of atom types
    (DrudeTransformDirect, "drude/transform/direct", 0),
    (DrudeTransformInverse, "drude/transform/inverse", 0),
    (DtReset, "dt/reset", 4),
    (Efield, "efield", 3),
    (EfieldTip4p, "efield/tip4p", 3),
    (Ehex, "ehex", 2),
    (ElectrodeConp, "electrode/conp", 2),
    (ElectrodeConq, "electrode/conq", 2),
    (ElectrodeThermo, "electrode/thermo", 5),
    (ElectronStopping, "electron/stopping", 2),
    (ElectronStoppingFit, "electron/stopping/fit", 3),
    (Enforce2d, "enforce2d", 0),
    (EosCv, "eos/cv", 1),
    (EosTable, "eos/table", 3),
    (EosTableRx, "eos/table/rx", 3),
    (Evaporate, "evaporate", 4),
    (External, "external", 2),
    (Ffl, "ffl", 4),
    (FilterCorotate, "filter/corotate", 0),
    (FlowGauss, "flow/gauss", 3),
    (Freeze, "freeze", 0),
    (Gcmc, "gcmc", 8),
    (Gld, "gld", 7),
    (Gle, "gle", 5),
    (Gravity, "gravity", 2),
    (Grem, "grem", 4),
    (Halt, "halt", 4),
    (Heat, "heat", 2),
    (HeatFlow, "heat/flow", 2),
    (HyperGlobal, "hyper/global", 4),
    (HyperLocal, "hyper/local", 7),
    (Imd, "imd", 2),
    (Indent, "indent", 1),
    (Ipi, "ipi", 2),
    (Langevin, "langevin", 4),
    (LangevinDrude, "langevin/drude", 6),
    (LangevinEff, "langevin/eff", 4),
    (LangevinSpin, "langevin/spin", 3),
    (LbFluid, "lb/fluid", 3),
    (LbMomentum, "lb/momentum", 1),
    (LbViscous, "lb/viscous", 0),
    (Lineforce, "lineforce", 3),
    (Manifoldforce, "manifoldforce", 2),
    (MdiQm, "mdi/qm", 0),
    (MdiQmmm, "mdi/qmmm", 1),
    (MesoMove, "meso/move", 4), // Each style has at least 3 args
    (MolSwap, "mol/swap", 6),
    (Momentum, "momentum", 1),
    (MomentumChunk, "momentum/chunk", 2),
    (Move, "move", 4), // Each style has at least 3 args
    (Mscg, "mscg", 1),
    (Msst, "msst", 2),
    (MvvDpd, "mvv/dpd", 0),
    (MvvEdpd, "mvv/edpd", 0),
    (MvvTdpd, "mvv/tdpd", 0),
    (Neb, "neb", 1),
    (NebSpin, "neb/spin", 1),
    (Nvt, "nvt", 0),
    (Npt, "npt", 0),
    (Nph, "nph", 0),
    (NvtEff, "nvt/eff", 0),
    (NptEff, "npt/eff", 0),
    (NphEff, "nph/eff", 0),
    (NvtUef, "nvt/uef", 0),
    (NptUef, "npt/uef", 0),
    (NphAsphere, "nph/asphere", 0),
    (NphBody, "nph/body", 0),
    (NphSphere, "nph/sphere", 0),
    (Nphug, "nphug", 0),
    (NptAsphere, "npt/asphere", 0),
    (NptBody, "npt/body", 0),
    (NptCauchy, "npt/cauchy", 0),
    (NptSphere, "npt/sphere", 0),
    (Numdiff, "numdiff", 2),
    (NumdiffVirial, "numdiff/virial", 2),
    (Nve, "nve", 0),
    (NveAsphere, "nve/asphere", 0),
    (NveAsphereNoforce, "nve/asphere/noforce", 0),
    (NveAwpmd, "nve/awpmd", 0),
    (NveBody, "nve/body", 0),
    (NveBpmSphere, "nve/bpm/sphere", 0),
    (NveDot, "nve/dot", 0),
    (NveDotcLangevin, "nve/dotc/langevin", 4),
    (NveEff, "nve/eff", 0),
    (NveLimit, "nve/limit", 1),
    (NveLine, "nve/line", 0),
    (NveManifoldRattle, "nve/manifold/rattle", 3),
    (NveNoforce, "nve/noforce", 0),
    (NveSphere, "nve/sphere", 0),
    (NveSpin, "nve/spin", 0),
    (NveTri, "nve/tri", 0),
    (Nvk, "nvk", 0),
    (NvtAsphere, "nvt/asphere", 0),
    (NvtBody, "nvt/body", 0),
    (NvtManifoldRattle, "nvt/manifold/rattle", 3),
    (NvtSllod, "nvt/sllod", 0),
    (NvtSllodEff, "nvt/sllod/eff", 0),
    (NvtSphere, "nvt/sphere", 0),
    (Oneway, "oneway", 3),
    (OrientFcc, "orient/fcc", 8),
    (OrientBcc, "orient/bcc", 8),
    (OrientEco, "orient/eco", 4),
    (Pafi, "pafi", 4),
    (Pair, "pair", 4),
    (Phonon, "phonon", 5),
    (PimdLangevin, "pimd/langevin", 0),
    (PimdNvt, "pimd/nvt", 0),
    (Planeforce, "planeforce", 3),
    (Plumed, "plumed", 0),
    (Poems, "poems", 0),
    (PolarizeBemGmres, "polarize/bem/gmres", 2),
    (PolarizeBemIcc, "polarize/bem/icc", 2),
    (PolarizeFunctional, "polarize/functional", 2),
    (Pour, "pour", 3),
    (PrecessionSpin, "precession/spin", 5),
    (PressBerendsen, "press/berendsen", 0),
    (Print, "print", 2),
    (PropelSelf, "propel/self", 2),
    (PropertyAtom, "property/atom", 1),
    (PythonInvoke, "python/invoke", 3),
    (PythonMove, "python/move", 1),
    (Qbmsst, "qbmsst", 2),
    (QeqPoint, "qeq/point", 5),
    (QeqShielded, "qeq/shielded", 5),
    (QeqSlater, "qeq/slater", 5),
    (QeqDynamic, "qeq/dynamic", 5),
    (QeqFire, "qeq/fire", 5),
    (QeqComb, "qeq/comb", 2),
    (QeqReaxff, "qeq/reaxff", 5),
    (Qmmm, "qmmm", 0),
    (Qtb, "qtb", 0),
    (ReaxffBonds, "reaxff/bonds", 2),
    (ReaxffSpecies, "reaxff/species", 4),
    (Recenter, "recenter", 3),
    (Restrain, "restrain", 0),
    (Rhok, "rhok", 5),
    (Rigid, "rigid", 1),
    (RigidNve, "rigid/nve", 1),
    (RigidNvt, "rigid/nvt", 1),
    (RigidNpt, "rigid/npt", 1),
    (RigidNph, "rigid/nph", 1),
    (RigidSmall, "rigid/small", 1),
    (RigidNveSmall, "rigid/nve/small", 1),
    (RigidNvtSmall, "rigid/nvt/small", 1),
    (RigidNptSmall, "rigid/npt/small", 1),
    (RigidNphSmall, "rigid/nph/small", 1),
    (RigidMeso, "rigid/meso", 2),
    (Rx, "rx", 4),
    (SaedVtk, "saed/vtk", 4),
    (Setforce, "setforce", 3),
    (SetforceSpin, "setforce/spin", 3),
    (Sgcmc, "sgcmc", 4),
    (Shake, "shake", 5),
    (Rattle, "rattle", 5),
    (Shardlow, "shardlow", 0),
    (Smd, "smd", 2),
    (SmdAdjust_dt, "smd/adjust_dt", 1),
    (SmdIntegrate_tlsph, "smd/integrate_tlsph", 0),
    (SmdIntegrate_ulsph, "smd/integrate_ulsph", 0),
    (SmdMove_tri_surf, "smd/move_tri_surf", 0),
    (SmdSetvel, "smd/setvel", 3),
    (SmdWall_surface, "smd/wall_surface", 3),
    (Sph, "sph", 0),
    (SphStationary, "sph/stationary", 0),
    (Spring, "spring", 6), // TODO: Check that zero isn't just acceptable?
    (SpringChunk, "spring/chunk", 3),
    (SpringRg, "spring/rg", 2),
    (SpringSelf, "spring/self", 2),
    (Srd, "srd", 5),
    (StoreForce, "store/force", 0),
    (StoreState, "store/state", 2),
    (TempBerendsen, "temp/berendsen", 3),
    (TempCsvr, "temp/csvr", 4),
    (TempCsld, "temp/csld", 4),
    (TempRescale, "temp/rescale", 5),
    (TempRescaleEff, "temp/rescale/eff", 5),
    (Tfmc, "tfmc", 3),
    (TgnvtDrude, "tgnvt/drude", 0),
    (TgnptDrude, "tgnpt/drude", 0),
    (ThermalConductivity, "thermal/conductivity", 3),
    (TiSpring, "ti/spring", 3),
    (Tmd, "tmd", 4),
    (Ttm, "ttm", 10),
    (TtmGrid, "ttm/grid", 10),
    (TtmMod, "ttm/mod", 5),
    (TuneKspace, "tune/kspace", 1),
    (Vector, "vector", 3),
    (Viscosity, "viscosity", 4),
    (Viscous, "viscous", 1),
    (ViscousSphere, "viscous/sphere", 1),
    // TODO: Really Not sure about these wall arg counts...
    (WallLj93, "wall/lj93", 2),
    (WallLj126, "wall/lj126", 2),
    (WallLj1043, "wall/lj1043", 2),
    (WallColloid, "wall/colloid", 2),
    (WallHarmonic, "wall/harmonic", 2),
    (WallLepton, "wall/lepton", 2),
    (WallMorse, "wall/morse", 2),
    (WallTable, "wall/table", 4),
    (WallBodyPolygon, "wall/body/polygon", 5),
    (WallBodyPolyhedron, "wall/body/polyhedron", 5),
    (WallEes, "wall/ees", 5),
    (WallRegionEes, "wall/region/ees", 4),
    (WallGran, "wall/gran", 4),
    (WallGranRegion, "wall/gran/region", 4),
    (WallPiston, "wall/piston", 1),
    (WallReflect, "wall/reflect", 2),
    (WallReflectStochastic, "wall/reflect/stochastic", 11), // TODO: Check...
    (WallRegion, "wall/region", 5),
    (WallSrd, "wall/srd", 2),
    (Widom, "widom", 5)
);
