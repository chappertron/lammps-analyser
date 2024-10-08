# compute smd/plastic/strain/rate command

## Syntax

    compute ID group-ID smd/plastic/strain/rate

-   ID, group-ID are documented in [compute](compute) command
-   smd/plastic/strain/rate = style name of this compute command

## Examples

``` LAMMPS
compute 1 all smd/plastic/strain/rate
```

## Description

Define a computation that outputs the time rate of the equivalent
plastic strain. This command is only meaningful if a material model with
plasticity is defined.

See [this PDF guide](PDF/MACHDYN_LAMMPS_userguide.pdf)\_ to use Smooth
Mach Dynamics in LAMMPS.

**Output Info:**

This compute calculates a per-particle vector, which can be accessed by
any command that uses per-particle values from a compute as input. See
the [Howto output](Howto_output) page for an overview of LAMMPS output
options.

The per-particle values will be given in [units](units) of one over
time.

## Restrictions

This compute is part of the MACHDYN package. It is only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info. This compute can only be
used for particles which interact via the updated Lagrangian or total
Lagrangian SPH pair styles.

## Related commands

[smd/plastic/strain](compute_smd_plastic_strain),
[smd/tlsph/strain/rate](compute_smd_tlsph_strain_rate),
[smd/tlsph/strain](compute_smd_tlsph_strain)

## Default

none
