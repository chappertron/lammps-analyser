# Overview of LAMMPS

LAMMPS is a classical molecular dynamics (MD) code that models ensembles
of particles in a liquid, solid, or gaseous state. It can model atomic,
polymeric, biological, solid-state (metals, ceramics, oxides), granular,
coarse-grained, or macroscopic systems using a variety of interatomic
potentials (force fields) and boundary conditions. It can model 2d or 3d
systems with sizes ranging from only a few particles up to billions.

LAMMPS can be built and run on single laptop or desktop machines, but is
designed for parallel computers. It will run in serial and on any
parallel machine that supports the
[MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface)\_
message-passing library. This includes shared-memory multicore,
multi-CPU servers and distributed-memory clusters and supercomputers.
Parts of LAMMPS also support [OpenMP
multi-threading](https://www.openmp.org)\_, vectorization, and GPU
acceleration.

LAMMPS is written in C++ and requires a compiler that is at least
compatible with the C++-11 standard. Earlier versions were written in
F77, F90, and C++-98. See the [History
page](https://www.lammps.org/history.html)\_ of the website for details.
All versions can be downloaded as source code from the [LAMMPS
website](https://www.lammps.org)\_.

LAMMPS is designed to be easy to modify or extend with new capabilities,
such as new force fields, atom types, boundary conditions, or
diagnostics. See the [Modify]{.title-ref} section of for more details.

In the most general sense, LAMMPS integrates Newton\'s equations of
motion for a collection of interacting particles. A single particle can
be an atom or molecule or electron, a coarse-grained cluster of atoms,
or a mesoscopic or macroscopic clump of material. The interaction models
that LAMMPS includes are mostly short-ranged in nature; some long-range
models are included as well.

LAMMPS uses neighbor lists to keep track of nearby particles. The lists
are optimized for systems with particles that are repulsive at short
distances, so that the local density of particles never becomes too
large. This is in contrast to methods used for modeling plasma or
gravitational bodies (like galaxy formation).

On parallel machines, LAMMPS uses spatial-decomposition techniques with
MPI parallelization to partition the simulation domain into subdomains
of equal computational cost, one of which is assigned to each processor.
Processors communicate and store \"ghost\" atom information for atoms
that border their subdomain. Multi-threading parallelization and GPU
acceleration with particle-decomposition can be used in addition.
