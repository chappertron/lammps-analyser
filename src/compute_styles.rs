macro_rules! compute_styles {
    ($(($variant:tt, $lit:literal, $nargs:literal) ),+) => {

        #[allow(non_camel_case_types)]
        #[derive(Debug,Default,Eq,PartialEq,Clone,Copy)]
        pub enum ComputeStyle {
            $($variant,)+
            #[default] // Because no other varitant can be default with macro def
            InvalidComputeStyle,
        }

        impl ComputeStyle {
            /// Returns the minimum number of positional arguments the compute style takes
            /// Note: This does not include the required arguments common to all computes, i.e.
            /// the `compute` keyword, the compute id, the group id and the compute stlye name
            pub const fn n_positional_args(&self) -> usize {
                match self {
                    $(ComputeStyle::$variant => $nargs,)+
                    ComputeStyle::InvalidComputeStyle => 0,
                }

            }
        }

        // TODO: Decide if try from or From is more appropriate
        impl From<&str> for ComputeStyle {
            fn from(value: &str) -> Self {
                match value {
                    $($lit => ComputeStyle::$variant,)+
                    _ => ComputeStyle::InvalidComputeStyle,
                }
            }
        }
        impl std::fmt::Display for ComputeStyle {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                match self {
                    $(ComputeStyle::$variant => write!(f, $lit),)+
                    ComputeStyle::InvalidComputeStyle => write!(f, "Invalid Compute Style"),
                }
            }
        }
    };
}

compute_styles!(
    (AcklandAtom, "ackland/atom", 0),
    (Adf, "adf", 0),
    (Angle, "angle", 0),
    (AngleLocal, "angle/local", 0),
    (AngmomChunk, "angmom/chunk", 0),
    (AveSphereAtom, "ave/sphere/atom", 0),
    (BasalAtom, "basal/atom", 0),
    (BodyLocal, "body/local", 0),
    (Bond, "bond", 0),
    (BondLocal, "bond/local", 0),
    (BornMatrix, "born/matrix", 0),
    (CentroAtom, "centro/atom", 0),
    (ChunkAtom, "chunk/atom", 0),
    (ChunkSpreadAtom, "chunk/spread/atom", 0),
    (ClusterAtom, "cluster/atom", 0),
    (FragmentAtom, "fragment/atom", 0),
    (AggregateAtom, "aggregate/atom", 0),
    (CnaAtom, "cna/atom", 0),
    (CnpAtom, "cnp/atom", 0),
    (Com, "com", 0),
    (ComChunk, "com/chunk", 0),
    (ContactAtom, "contact/atom", 0),
    (CoordAtom, "coord/atom", 0),
    (CountType, "count/type", 0),
    (DamageAtom, "damage/atom", 0),
    (Dihedral, "dihedral", 0),
    (DihedralLocal, "dihedral/local", 0),
    (DilatationAtom, "dilatation/atom", 0),
    (Dipole, "dipole", 0),
    (DipoleTip4p, "dipole/tip4p", 0),
    (DipoleChunk, "dipole/chunk", 0),
    (DipoleTip4pChunk, "dipole/tip4p/chunk", 0),
    (DisplaceAtom, "displace/atom", 0),
    (Dpd, "dpd", 0),
    (DpdAtom, "dpd/atom", 0),
    (EdpdTempAtom, "edpd/temp/atom", 0),
    (EfieldAtom, "efield/atom", 0),
    (EfieldWolfAtom, "efield/wolf/atom", 0),
    (EntropyAtom, "entropy/atom", 0),
    (ErotateAsphere, "erotate/asphere", 0),
    (ErotateRigid, "erotate/rigid", 0),
    (ErotateSphere, "erotate/sphere", 0),
    (ErotateSphereAtom, "erotate/sphere/atom", 0),
    (EventDisplace, "event/displace", 0),
    (Fabric, "fabric", 0),
    (Fep, "fep", 0),
    (FepTa, "fep/ta", 0),
    (GlobalAtom, "global/atom", 0),
    (GroupGroup, "group/group", 0),
    (Gyration, "gyration", 0),
    (GyrationChunk, "gyration/chunk", 0),
    (GyrationShape, "gyration/shape", 0),
    (GyrationShapeChunk, "gyration/shape/chunk", 0),
    (HeatFlux, "heat/flux", 0),
    (HexorderAtom, "hexorder/atom", 0),
    (Hma, "hma", 0),
    (Improper, "improper", 0),
    (ImproperLocal, "improper/local", 0),
    (InertiaChunk, "inertia/chunk", 0),
    (Ke, "ke", 0),
    (KeAtom, "ke/atom", 0),
    (KeAtomEff, "ke/atom/eff", 0),
    (KeEff, "ke/eff", 0),
    (KeRigid, "ke/rigid", 0),
    (Mliap, "mliap", 0),
    (Modify, "modify", 0),
    (Momentum, "momentum", 0),
    (Msd, "msd", 0),
    (MsdChunk, "msd/chunk", 0),
    (MsdNongauss, "msd/nongauss", 0),
    (NbondAtom, "nbond/atom", 0),
    (OmegaChunk, "omega/chunk", 0),
    (OrientorderAtom, "orientorder/atom", 0),
    (Pair, "pair", 0),
    (PairLocal, "pair/local", 0),
    (Pe, "pe", 0),
    (PeAtom, "pe/atom", 0),
    (PlasticityAtom, "plasticity/atom", 0),
    (Pressure, "pressure", 0),
    (PressureAlchemy, "pressure/alchemy", 0),
    (PressureUef, "pressure/uef", 0),
    (PropertyAtom, "property/atom", 0),
    (PropertyChunk, "property/chunk", 0),
    (PropertyGrid, "property/grid", 0),
    (PropertyLocal, "property/local", 0),
    (PtmAtom, "ptm/atom", 0),
    (Rdf, "rdf", 0),
    (Reduce, "reduce", 0),
    (ReduceRegion, "reduce/region", 0),
    (ReduceChunk, "reduce/chunk", 0),
    (RigidLocal, "rigid/local", 0),
    (Saed, "saed", 0),
    (Slice, "slice", 0),
    (SmdContactRadius, "smd/contact/radius", 0),
    (SmdDamage, "smd/damage", 0),
    (SmdHourglassError, "smd/hourglass/error", 0),
    (SmdInternalEnergy, "smd/internal/energy", 0),
    (SmdPlasticStrain, "smd/plastic/strain", 0),
    (SmdPlasticStrainRate, "smd/plastic/strain/rate", 0),
    (SmdRho, "smd/rho", 0),
    (SmdTlsphDefgrad, "smd/tlsph/defgrad", 0),
    (SmdTlsphDt, "smd/tlsph/dt", 0),
    (SmdTlsphNumNeighs, "smd/tlsph/num/neighs", 0),
    (SmdTlsphShape, "smd/tlsph/shape", 0),
    (SmdTlsphStrain, "smd/tlsph/strain", 0),
    (SmdTlsphStrainRate, "smd/tlsph/strain/rate", 0),
    (SmdTlsphStress, "smd/tlsph/stress", 0),
    (SmdTriangleVertices, "smd/triangle/vertices", 0),
    (SmdUlsphEffm, "smd/ulsph/effm", 0),
    (SmdUlsphNumNeighs, "smd/ulsph/num/neighs", 0),
    (SmdUlsphStrain, "smd/ulsph/strain", 0),
    (SmdUlsphStrainRate, "smd/ulsph/strain/rate", 0),
    (SmdUlsphStress, "smd/ulsph/stress", 0),
    (SmdVol, "smd/vol", 0),
    (SnaAtom, "sna/atom", 0),
    (SnadAtom, "snad/atom", 0),
    (SnavAtom, "snav/atom", 0),
    (Snap, "snap", 0),
    (SnaGrid, "sna/grid", 0),
    (SnaGridLocal, "sna/grid/local", 0),
    (SphEAtom, "sph/e/atom", 0),
    (SphRhoAtom, "sph/rho/atom", 0),
    (SphTAtom, "sph/t/atom", 0),
    (Spin, "spin", 0),
    (StressAtom, "stress/atom", 0),
    (CentroidStressAtom, "centroid/stress/atom", 0),
    (StressCartesian, "stress/cartesian", 0),
    (StressCylinder, "stress/cylinder", 0),
    (StressSpherical, "stress/spherical", 0),
    (StressMop, "stress/mop", 0),
    (StressMopProfile, "stress/mop/profile", 0),
    (ForceTally, "force/tally", 0),
    (HeatFluxTally, "heat/flux/tally", 0),
    (HeatFluxVirialTally, "heat/flux/virial/tally", 0),
    (PeTally, "pe/tally", 0),
    (PeMolTally, "pe/mol/tally", 0),
    (StressTally, "stress/tally", 0),
    (TdpdCcAtom, "tdpd/cc/atom", 0),
    (Temp, "temp", 0),
    (TempAsphere, "temp/asphere", 0),
    (TempBody, "temp/body", 0),
    (TempChunk, "temp/chunk", 0),
    (TempCom, "temp/com", 0),
    (TempCs, "temp/cs", 0),
    (TempDeform, "temp/deform", 0),
    (TempDeformEff, "temp/deform/eff", 0),
    (TempDrude, "temp/drude", 0),
    (TempEff, "temp/eff", 0),
    (TempPartial, "temp/partial", 0),
    (TempProfile, "temp/profile", 0),
    (TempRamp, "temp/ramp", 0),
    (TempRegion, "temp/region", 0),
    (TempRegionEff, "temp/region/eff", 0),
    (TempRotate, "temp/rotate", 0),
    (TempSphere, "temp/sphere", 0),
    (TempUef, "temp/uef", 0),
    (Ti, "ti", 0),
    (TorqueChunk, "torque/chunk", 0),
    (Vacf, "vacf", 0),
    (VcmChunk, "vcm/chunk", 0),
    (ViscosityCos, "viscosity/cos", 0),
    (VoronoiAtom, "voronoi/atom", 0),
    (Xrd, "xrd", 0)
);
