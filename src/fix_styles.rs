#[allow(non_camel_case_types)] // For some types

/// Enum for all the possible styles of a fix.
/// TODO Seperate by package
#[derive(Default, Debug, Eq, PartialEq, Clone, Copy)]
pub enum FixStyle {
    AccelerateCos,
    Acks2Reaxff,
    Adapt,
    AdaptFep,
    Addforce,
    Addtorque,
    Alchemy,
    AmoebaBitorsion,
    AmoebaPitorsion,
    AppendAtoms,
    Atc,
    AtomSwap,
    AveAtom,
    AveChunk,
    AveCorrelate,
    AveCorrelateLong,
    AveGrid,
    AveHisto,
    AveHistoWeight,
    /// Deprecated
    AveSpatial,
    AveTime,
    Aveforce,
    Balance,
    Bocs,
    BondBreak,
    BondCreate,
    BondCreateAngle,
    BondReact,
    BondSwap,
    BoxRelax,
    Brownian,
    BrownianSphere,
    BrownianAsphere,
    ChargeRegulation,
    Cmap,
    Colvars,
    Controller,
    DampingCundall,
    Deform,
    Deposit,
    DpdEnergy,
    EdpdSource,
    TdpdSource,
    Drag,
    Drude,
    DrudeTransformDirect,
    DrudeTransformInverse,
    DtReset,
    Efield,
    EfieldTip4p,
    Ehex,
    ElectrodeConp,
    ElectrodeConq,
    ElectrodeThermo,
    ElectronStopping,
    ElectronStoppingFit,
    Enforce2d,
    EosCv,
    EosTable,
    EosTableRx,
    Evaporate,
    External,
    Ffl,
    FilterCorotate,
    FlowGauss,
    Freeze,
    Gcmc,
    Gld,
    Gle,
    Gravity,
    Grem,
    Halt,
    Heat,
    HeatFlow,
    HyperGlobal,
    HyperLocal,
    Imd,
    Indent,
    Ipi,
    Langevin,
    LangevinDrude,
    LangevinEff,
    LangevinSpin,
    LbFluid,
    LbMomentum,
    LbViscous,
    Lineforce,
    Manifoldforce,
    MdiQm,
    MdiQmmm,
    MesoMove,
    MolSwap,
    Momentum,
    MomentumChunk,
    Move,
    Mscg,
    Msst,
    MvvDpd,
    MvvEdpd,
    MvvTdpd,
    Neb,
    NebSpin,
    Nvt,
    Npt,
    Nph,
    NvtEff,
    NptEff,
    NphEff,
    NvtUef,
    NptUef,
    NphAsphere,
    NphBody,
    NphSphere,
    Nphug,
    NptAsphere,
    NptBody,
    NptCauchy,
    NptSphere,
    Numdiff,
    NumdiffVirial,
    #[default]
    Nve,
    NveAsphere,
    NveAsphereNoforce,
    NveAwpmd,
    NveBody,
    NveBpmSphere,
    NveDot,
    NveDotcLangevin,
    NveEff,
    NveLimit,
    NveLine,
    NveManifoldRattle,
    NveNoforce,
    NveSphere,
    NveSpin,
    NveTri,
    Nvk,
    NvtAsphere,
    NvtBody,
    NvtManifoldRattle,
    NvtSllod,
    NvtSllodEff,
    NvtSphere,
    Oneway,
    OrientFcc,
    OrientBcc,
    OrientEco,
    Pafi,
    Pair,
    Phonon,
    PimdLangevin,
    PimdNvt,
    Planeforce,
    Plumed,
    Poems,
    PolarizeBemGmres,
    PolarizeBemIcc,
    PolarizeFunctional,
    Pour,
    PrecessionSpin,
    PressBerendsen,
    Print,
    PropelSelf,
    PropertyAtom,
    PythonInvoke,
    PythonMove,
    Qbmsst,
    QeqPoint,
    QeqShielded,
    QeqSlater,
    QeqDynamic,
    QeqFire,
    QeqComb,
    QeqReaxff,
    Qmmm,
    Qtb,
    ReaxffBonds,
    ReaxffSpecies,
    Recenter,
    Restrain,
    Rhok,
    Rigid,
    RigidNve,
    RigidNvt,
    RigidNpt,
    RigidNph,
    RigidSmall,
    RigidNveSmall,
    RigidNvtSmall,
    RigidNptSmall,
    RigidNphSmall,
    RigidMeso,
    Rx,
    SaedVtk,
    Setforce,
    SetforceSpin,
    Sgcmc,
    Shake,
    Rattle,
    Shardlow,
    Smd,
    SmdAdjust_dt,
    SmdIntegrate_tlsph,
    SmdIntegrate_ulsph,
    SmdMove_tri_surf,
    SmdSetvel,
    SmdWall_surface,
    Sph,
    SphStationary,
    Spring,
    SpringChunk,
    SpringRg,
    SpringSelf,
    Srd,
    StoreForce,
    StoreState,
    TempBerendsen,
    TempCsvr,
    TempCsld,
    TempRescale,
    TempRescaleEff,
    Tfmc,
    TgnvtDrude,
    TgnptDrude,
    ThermalConductivity,
    TiSpring,
    Tmd,
    Ttm,
    TtmGrid,
    TtmMod,
    TuneKspace,
    Vector,
    Viscosity,
    Viscous,
    ViscousSphere,
    WallLj93,
    WallLj126,
    WallLj1043,
    WallColloid,
    WallHarmonic,
    WallLepton,
    WallMorse,
    WallTable,
    WallBodyPolygon,
    WallBodyPolyhedron,
    WallEes,
    WallRegionEes,
    WallGran,
    WallGranRegion,
    WallPiston,
    WallReflect,
    WallReflectStochastic,
    WallRegion,
    WallSrd,
    Widom,
    InvalidFixStyle,
}

impl From<&str> for FixStyle {
    fn from(value: &str) -> Self {
        use FixStyle::*;
        match value {
            "accelerate/cos" => AccelerateCos,
            "acks2/reaxff" => Acks2Reaxff,
            "adapt" => Adapt,
            "adapt/fep" => AdaptFep,
            "addforce" => Addforce,
            "addtorque" => Addtorque,
            "alchemy" => Alchemy,
            "amoeba/bitorsion" => AmoebaBitorsion,
            "amoeba/pitorsion" => AmoebaPitorsion,
            "append/atoms" => AppendAtoms,
            "atc" => Atc,
            "atom/swap" => AtomSwap,
            "ave/atom" => AveAtom,
            "ave/chunk" => AveChunk,
            "ave/correlate" => AveCorrelate,
            "ave/correlate/long" => AveCorrelateLong,
            "ave/grid" => AveGrid,
            "ave/histo" => AveHisto,
            "ave/histo/weight" => AveHistoWeight,
            "ave/spatial" => AveSpatial,
            "ave/time" => AveTime,
            "aveforce" => Aveforce,
            "balance" => Balance,
            "bocs" => Bocs,
            "bond/break" => BondBreak,
            "bond/create" => BondCreate,
            "bond/create/angle" => BondCreateAngle,
            "bond/react" => BondReact,
            "bond/swap" => BondSwap,
            "box/relax" => BoxRelax,
            "brownian" => Brownian,
            "brownian/sphere" => BrownianSphere,
            "brownian/asphere" => BrownianAsphere,
            "charge/regulation" => ChargeRegulation,
            "cmap" => Cmap,
            "colvars" => Colvars,
            "controller" => Controller,
            "damping/cundall" => DampingCundall,
            "deform" => Deform,
            "deposit" => Deposit,
            "dpd/energy" => DpdEnergy,
            "edpd/source" => EdpdSource,
            "tdpd/source" => TdpdSource,
            "drag" => Drag,
            "drude" => Drude,
            "drude/transform/direct" => DrudeTransformDirect,
            "drude/transform/inverse" => DrudeTransformInverse,
            "dt/reset" => DtReset,
            "efield" => Efield,
            "efield/tip4p" => EfieldTip4p,
            "ehex" => Ehex,
            "electrode/conp" => ElectrodeConp,
            "electrode/conq" => ElectrodeConq,
            "electrode/thermo" => ElectrodeThermo,
            "electron/stopping" => ElectronStopping,
            "electron/stopping/fit" => ElectronStoppingFit,
            "enforce2d" => Enforce2d,
            "eos/cv" => EosCv,
            "eos/table" => EosTable,
            "eos/table/rx" => EosTableRx,
            "evaporate" => Evaporate,
            "external" => External,
            "ffl" => Ffl,
            "filter/corotate" => FilterCorotate,
            "flow/gauss" => FlowGauss,
            "freeze" => Freeze,
            "gcmc" => Gcmc,
            "gld" => Gld,
            "gle" => Gle,
            "gravity" => Gravity,
            "grem" => Grem,
            "halt" => Halt,
            "heat" => Heat,
            "heat/flow" => HeatFlow,
            "hyper/global" => HyperGlobal,
            "hyper/local" => HyperLocal,
            "imd" => Imd,
            "indent" => Indent,
            "ipi" => Ipi,
            "langevin" => Langevin,
            "langevin/drude" => LangevinDrude,
            "langevin/eff" => LangevinEff,
            "langevin/spin" => LangevinSpin,
            "lb/fluid" => LbFluid,
            "lb/momentum" => LbMomentum,
            "lb/viscous" => LbViscous,
            "lineforce" => Lineforce,
            "manifoldforce" => Manifoldforce,
            "mdi/qm" => MdiQm,
            "mdi/qmmm" => MdiQmmm,
            "meso/move" => MesoMove,
            "mol/swap" => MolSwap,
            "momentum" => Momentum,
            "momentum/chunk" => MomentumChunk,
            "move" => Move,
            "mscg" => Mscg,
            "msst" => Msst,
            "mvv/dpd" => MvvDpd,
            "mvv/edpd" => MvvEdpd,
            "mvv/tdpd" => MvvTdpd,
            "neb" => Neb,
            "neb/spin" => NebSpin,
            "nvt" => Nvt,
            "npt" => Npt,
            "nph" => Nph,
            "nvt/eff" => NvtEff,
            "npt/eff" => NptEff,
            "nph/eff" => NphEff,
            "nvt/uef" => NvtUef,
            "npt/uef" => NptUef,
            "nph/asphere" => NphAsphere,
            "nph/body" => NphBody,
            "nph/sphere" => NphSphere,
            "nphug" => Nphug,
            "npt/asphere" => NptAsphere,
            "npt/body" => NptBody,
            "npt/cauchy" => NptCauchy,
            "npt/sphere" => NptSphere,
            "numdiff" => Numdiff,
            "numdiff/virial" => NumdiffVirial,
            "nve" => Nve,
            "nve/asphere" => NveAsphere,
            "nve/asphere/noforce" => NveAsphereNoforce,
            "nve/awpmd" => NveAwpmd,
            "nve/body" => NveBody,
            "nve/bpm/sphere" => NveBpmSphere,
            "nve/dot" => NveDot,
            "nve/dotc/langevin" => NveDotcLangevin,
            "nve/eff" => NveEff,
            "nve/limit" => NveLimit,
            "nve/line" => NveLine,
            "nve/manifold/rattle" => NveManifoldRattle,
            "nve/noforce" => NveNoforce,
            "nve/sphere" => NveSphere,
            "nve/spin" => NveSpin,
            "nve/tri" => NveTri,
            "nvk" => Nvk,
            "nvt/asphere" => NvtAsphere,
            "nvt/body" => NvtBody,
            "nvt/manifold/rattle" => NvtManifoldRattle,
            "nvt/sllod" => NvtSllod,
            "nvt/sllod/eff" => NvtSllodEff,
            "nvt/sphere" => NvtSphere,
            "oneway" => Oneway,
            "orient/fcc" => OrientFcc,
            "orient/bcc" => OrientBcc,
            "orient/eco" => OrientEco,
            "pafi" => Pafi,
            "pair" => Pair,
            "phonon" => Phonon,
            "pimd/langevin" => PimdLangevin,
            "pimd/nvt" => PimdNvt,
            "planeforce" => Planeforce,
            "plumed" => Plumed,
            "poems" => Poems,
            "polarize/bem/gmres" => PolarizeBemGmres,
            "polarize/bem/icc" => PolarizeBemIcc,
            "polarize/functional" => PolarizeFunctional,
            "pour" => Pour,
            "precession/spin" => PrecessionSpin,
            "press/berendsen" => PressBerendsen,
            "print" => Print,
            "propel/self" => PropelSelf,
            "property/atom" => PropertyAtom,
            "python/invoke" => PythonInvoke,
            "python/move" => PythonMove,
            "qbmsst" => Qbmsst,
            "qeq/point" => QeqPoint,
            "qeq/shielded" => QeqShielded,
            "qeq/slater" => QeqSlater,
            "qeq/dynamic" => QeqDynamic,
            "qeq/fire" => QeqFire,
            "qeq/comb" => QeqComb,
            "qeq/reaxff" => QeqReaxff,
            "qmmm" => Qmmm,
            "qtb" => Qtb,
            "reaxff/bonds" => ReaxffBonds,
            "reaxff/species" => ReaxffSpecies,
            "recenter" => Recenter,
            "restrain" => Restrain,
            "rhok" => Rhok,
            "rigid" => Rigid,
            "rigid/nve" => RigidNve,
            "rigid/nvt" => RigidNvt,
            "rigid/npt" => RigidNpt,
            "rigid/nph" => RigidNph,
            "rigid/small" => RigidSmall,
            "rigid/nve/small" => RigidNveSmall,
            "rigid/nvt/small" => RigidNvtSmall,
            "rigid/npt/small" => RigidNptSmall,
            "rigid/nph/small" => RigidNphSmall,
            "rigid/meso" => RigidMeso,
            "rx" => Rx,
            "saed/vtk" => SaedVtk,
            "setforce" => Setforce,
            "setforce/spin" => SetforceSpin,
            "sgcmc" => Sgcmc,
            "shake" => Shake,
            "rattle" => Rattle,
            "shardlow" => Shardlow,
            "smd" => Smd,
            "smd/adjust_dt" => SmdAdjust_dt,
            "smd/integrate_tlsph" => SmdIntegrate_tlsph,
            "smd/integrate_ulsph" => SmdIntegrate_ulsph,
            "smd/move_tri_surf" => SmdMove_tri_surf,
            "smd/setvel" => SmdSetvel,
            "smd/wall_surface" => SmdWall_surface,
            "sph" => Sph,
            "sph/stationary" => SphStationary,
            "spring" => Spring,
            "spring/chunk" => SpringChunk,
            "spring/rg" => SpringRg,
            "spring/self" => SpringSelf,
            "srd" => Srd,
            "store/force" => StoreForce,
            "store/state" => StoreState,
            "temp/berendsen" => TempBerendsen,
            "temp/csvr" => TempCsvr,
            "temp/csld" => TempCsld,
            "temp/rescale" => TempRescale,
            "temp/rescale/eff" => TempRescaleEff,
            "tfmc" => Tfmc,
            "tgnvt/drude" => TgnvtDrude,
            "tgnpt/drude" => TgnptDrude,
            "thermal/conductivity" => ThermalConductivity,
            "ti/spring" => TiSpring,
            "tmd" => Tmd,
            "ttm" => Ttm,
            "ttm/grid" => TtmGrid,
            "ttm/mod" => TtmMod,
            "tune/kspace" => TuneKspace,
            "vector" => Vector,
            "viscosity" => Viscosity,
            "viscous" => Viscous,
            "viscous/sphere" => ViscousSphere,
            "wall/lj93" => WallLj93,
            "wall/lj126" => WallLj126,
            "wall/lj1043" => WallLj1043,
            "wall/colloid" => WallColloid,
            "wall/harmonic" => WallHarmonic,
            "wall/lepton" => WallLepton,
            "wall/morse" => WallMorse,
            "wall/table" => WallTable,
            "wall/body/polygon" => WallBodyPolygon,
            "wall/body/polyhedron" => WallBodyPolyhedron,
            "wall/ees" => WallEes,
            "wall/region/ees" => WallRegionEes,
            "wall/gran" => WallGran,
            "wall/gran/region" => WallGranRegion,
            "wall/piston" => WallPiston,
            "wall/reflect" => WallReflect,
            "wall/reflect/stochastic" => WallReflectStochastic,
            "wall/region" => WallRegion,
            "wall/srd" => WallSrd,
            "widom" => Widom,
            _ => InvalidFixStyle,
        }
    }
}
