macro_rules! compute_styles {
    ($(($variant:tt, $lit:literal, $nargs:literal) ),+) => {

        #[allow(non_camel_case_types)]
        #[derive(Debug,Default,Eq,PartialEq,Clone,Copy,Hash)]
        pub enum ComputeStyle {
            $($variant,)+
            #[default] // Because no other varitant can be default with macro def
            InvalidComputeStyle,
        }

        impl ComputeStyle {
            /// Returns the minimum number of positional arguments the compute style takes
            /// Note: This does not include the required arguments common to all computes, i.e.
            /// the `compute` keyword, the compute id, the group id and the compute style name
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
    (Adf, "adf", 1),
    (Angle, "angle", 0),
    (AngleLocal, "angle/local", 1), // One or more values may be appended.
    (AngmomChunk, "angmom/chunk", 1),
    (AveSphereAtom, "ave/sphere/atom", 0),
    (BasalAtom, "basal/atom", 0),
    (BodyLocal, "body/local", 1),
    (Bond, "bond", 0),
    (BondLocal, "bond/local", 1),
    (BornMatrix, "born/matrix", 0),
    (CentroAtom, "centro/atom", 1),
    (ChunkAtom, "chunk/atom", 1), // Requires at least the style argument.
    (ChunkSpreadAtom, "chunk/spread/atom", 1), // requires at least the compute chunk/atom id.
    (ClusterAtom, "cluster/atom", 1), // Cut-off
    (FragmentAtom, "fragment/atom", 0),
    (AggregateAtom, "aggregate/atom", 1), // Cut-off
    (CnaAtom, "cna/atom", 1),             // Cut-off
    (CnpAtom, "cnp/atom", 1),             // Cut-off
    (Com, "com", 0),
    (ComChunk, "com/chunk", 1), // chunkID
    (ContactAtom, "contact/atom", 0),
    (CompositionAtom, "composition/atom", 0), // Added in 21Nov2023
    (CoordAtom, "coord/atom", 1),             // Requires a style argument.
    (CountType, "count/type", 1),
    (DamageAtom, "damage/atom", 0),
    (Dihedral, "dihedral", 0),
    (DihedralLocal, "dihedral/local", 1),
    (DilatationAtom, "dilatation/atom", 0),
    (Dipole, "dipole", 0),
    (DipoleTip4p, "dipole/tip4p", 0),
    (DipoleChunk, "dipole/chunk", 1),            // chunkID
    (DipoleTip4pChunk, "dipole/tip4p/chunk", 1), // chunkID
    (DisplaceAtom, "displace/atom", 0),
    (Dpd, "dpd", 0),
    (DpdAtom, "dpd/atom", 0),
    (EdpdTempAtom, "edpd/temp/atom", 0),
    (EfieldAtom, "efield/atom", 0),
    (EfieldWolfAtom, "efield/wolf/atom", 1), // alpha is required.
    (EntropyAtom, "entropy/atom", 2),        // sigma and cut-off are required.
    (ErotateAsphere, "erotate/asphere", 0),
    (ErotateRigid, "erotate/rigid", 1), // rigid fix-ID
    (ErotateSphere, "erotate/sphere", 0),
    (ErotateSphereAtom, "erotate/sphere/atom", 0),
    (EventDisplace, "event/displace", 1), // threshold
    (Fabric, "fabric", 1),                // cutoff argument
    (Fep, "fep", 1),                      // temp
    (FepTa, "fep/ta", 3),                 // temp plane scale_factor
    (GlobalAtom, "global/atom", 2),
    (GroupGroup, "group/group", 1), // group2-ID
    (Gyration, "gyration", 0),
    (GyrationChunk, "gyration/chunk", 1),            // chunk-ID
    (GyrationShape, "gyration/shape", 1),            // compute-ID for compute gyration
    (GyrationShapeChunk, "gyration/shape/chunk", 1), // compute-ID for compute gyration/chunk
    (HeatFlux, "heat/flux", 3),                      // ke-ID pe-ID stress-ID
    (HexorderAtom, "hexorder/atom", 0),
    (Hma, "hma", 1), // temp-ID
    (Improper, "improper", 0),
    (ImproperLocal, "improper/local", 0), // 1 or more chi
    (InertiaChunk, "inertia/chunk", 1),   // chunk-ID
    (Ke, "ke", 0),
    (KeAtom, "ke/atom", 0),
    (KeAtomEff, "ke/atom/eff", 0),
    (KeEff, "ke/eff", 0),
    (KeRigid, "ke/rigid", 1), // fix-ID
    (Mliap, "mliap", 4),      // requires 2 or more keyword-value pairs.
    (Momentum, "momentum", 0),
    (Msd, "msd", 0),
    (MsdChunk, "msd/chunk", 1), // chunkID
    (MsdNongauss, "msd/nongauss", 0),
    (NbondAtom, "nbond/atom", 0),
    (OmegaChunk, "omega/chunk", 1),
    (OrientorderAtom, "orientorder/atom", 0),
    (Pace, "pace", 1), // filename
    (Pair, "pair", 1), // pstyle
    (PairLocal, "pair/local", 1),
    (Pe, "pe", 0),
    (PeAtom, "pe/atom", 0),
    (PlasticityAtom, "plasticity/atom", 0),
    (Pressure, "pressure", 1), // Temp compute ID. Can be NULL
    (PressureAlchemy, "pressure/alchemy", 1), // Temp compute ID. Can be NULL
    (PressureUef, "pressure/uef", 1), // Temp compute ID. Can be NULL
    (PropertyAtom, "property/atom", 1), // One or more properties.
    (PropertyChunk, "property/chunk", 2), // chunkID + One or more properties.
    (PropertyGrid, "property/grid", 4), // Nx Ny Nz + One or more properties.
    (PropertyLocal, "property/local", 1), // One or more properties.
    (PtmAtom, "ptm/atom", 3),  // structures threshold group2-ID
    (RattlersAtom, "rattlers/atom", 3), // cutoff zmin ntries
    (Rdf, "rdf", 1),           // Nbin
    (ReaxffAtom, "reaxff/atom", 1), // takes at least an attribute,
    (Reduce, "reduce", 0),
    (ReduceRegion, "reduce/region", 1), // region-ID
    (ReduceChunk, "reduce/chunk", 2),   // chunk-ID mode
    (RigidLocal, "rigid/local", 2),     // rigidID one or more properties
    (Saed, "saed", 2),                  // lambda + one or more types
    (SlcsaAtom, "slcsa/atom", 7),
    (Slice, "slice", 4), // Nstart Nstop Nskip + 1 or more inputs
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
    (SmdVol, "smd/vol", 2),
    // sna* styles all seem to require at least two arguments
    // The commands all being documented on asingle page makes it hard to understand.
    (SnaAtom, "sna/atom", 2),
    (SnadAtom, "snad/atom", 2),
    (SnavAtom, "snav/atom", 2),
    (Snap, "snap", 2),
    (SnaGrid, "sna/grid", 3),
    (SnaGridLocal, "sna/grid/local", 3),
    (SphEAtom, "sph/e/atom", 0),
    (SphRhoAtom, "sph/rho/atom", 0),
    (SphTAtom, "sph/t/atom", 0),
    (Spin, "spin", 0),
    (StressAtom, "stress/atom", 1), // temp-ID  compute-ID
    (CentroidStressAtom, "centroid/stress/atom", 1), // temp-ID  compute-ID
    (StressCartesian, "stress/cartesian", 4), // dim1 bw1 dim2 bw2
    (StressCylinder, "stress/cylinder", 4),
    (StressSpherical, "stress/spherical", 5),
    (StressMop, "stress/mop", 2),
    (StressMopProfile, "stress/mop/profile", 3),
    // All /tally styles take a group2-ID argument.
    (ForceTally, "force/tally", 1),
    (HeatFluxTally, "heat/flux/tally", 1),
    (HeatFluxVirialTally, "heat/flux/virial/tally", 1),
    (PeTally, "pe/tally", 1),
    (PeMolTally, "pe/mol/tally", 1),
    (StressTally, "stress/tally", 1),
    (TdpdCcAtom, "tdpd/cc/atom", 1), // index
    (Temp, "temp", 0),
    (TempAsphere, "temp/asphere", 0),
    (TempBody, "temp/body", 0),
    (TempChunk, "temp/chunk", 1),
    (TempCom, "temp/com", 0),
    (TempCs, "temp/cs", 2),
    (TempDeform, "temp/deform", 0),
    (TempDeformEff, "temp/deform/eff", 0),
    (TempDrude, "temp/drude", 0),
    (TempEff, "temp/eff", 0),
    (TempPartial, "temp/partial", 3), // 3 flags
    (TempProfile, "temp/profile", 5), // xflag yflag zflag + style + atleast one arg for each style
    (TempRamp, "temp/ramp", 6),
    (TempRegion, "temp/region", 1),        // region-ID
    (TempRegionEff, "temp/region/eff", 1), // region-ID
    (TempRotate, "temp/rotate", 0),
    (TempSphere, "temp/sphere", 0),
    (TempUef, "temp/uef", 0),
    (Ti, "ti", 0),
    (TorqueChunk, "torque/chunk", 1), // chunk-ID
    (Vacf, "vacf", 0),
    (VcmChunk, "vcm/chunk", 1),
    (ViscosityCos, "viscosity/cos", 0),
    (VoronoiAtom, "voronoi/atom", 0),
    (Xrd, "xrd", 2) // lambda +  1 or more types and keywords.
);

impl ComputeStyle {
    pub fn is_invalid(&self) -> bool {
        matches!(self, ComputeStyle::InvalidComputeStyle)
    }
}

#[cfg(test)]
mod tests {
    use crate::compute_styles::ComputeStyle;

    #[test]
    fn all_computes_defined() {
        // Read all the computes from the txt file ensure they all can be converted and have all
        // been covered

        let computes = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/docs_extract/computes.txt"
        ));

        for compute in computes.lines() {
            dbg!(compute);

            // Ensure that call computes get parsed
            assert_ne!(
                ComputeStyle::from(compute),
                ComputeStyle::InvalidComputeStyle
            )
        }
    }
}
