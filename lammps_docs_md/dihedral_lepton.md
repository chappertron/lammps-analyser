# dihedral_style lepton command

Accelerator Variants: *lepton/omp*

## Syntax

``` LAMMPS
dihedral_style lepton
```

## Examples

``` LAMMPS
dihedral_style lepton

dihedral_coeff  1 "k*(1 + d*cos(n*phi)); k=75.0; d=1; n=2"
dihedral_coeff  2 "45*(1-cos(4*phi))"
dihedral_coeff  2 "k2*cos(phi) + k3*cos(phi)^2; k2=100.0"
dihedral_coeff  3 "k*(phi-phi0)^2; k=85.0; phi0=120.0"
```

## Description

::: versionadded
8Feb2023
:::

Dihedral style *lepton* computes dihedral interactions between four
atoms forming a dihedral angle with a custom potential function. The
potential function must be provided as an expression string using
\"phi\" as the dihedral angle variable. For example
[\"200.0\*(phi-120.0)\^2\"]{.title-ref} represents a [quadratic
dihedral](dihedral_quadratic) potential around a 120 degree dihedral
angle with a force constant *K* of 200.0 energy units:

$$U_{dihedral,i} = K (\phi_i - \phi_0)^2$$

The [Lepton library](https://simtk.org/projects/lepton)\_, that the
*lepton* dihedral style interfaces with, evaluates this expression
string at run time to compute the pairwise energy. It also creates an
analytical representation of the first derivative of this expression
with respect to \"phi\" and then uses that to compute the force between
the dihedral atoms as defined by the topology data.

The potential function expression for each dihedral type is provided via
the [dihedral_coeff](dihedral_coeff) command as in the example above, or
in the data file or restart files read by the [read_data](read_data) or
[read_restart](read_restart) commands. The expression is in energy
units.

The Lepton expression must be either enclosed in quotes or must not
contain any whitespace so that LAMMPS recognizes it as a single keyword.
More on valid Lepton expressions below. Dihedral angles are internally
computed in radians and thus the expression must assume \"phi\" is in
radians.

------------------------------------------------------------------------

Lepton expression syntax and features
\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"\"

Lepton supports the following operators in expressions:

The following mathematical functions are available:

Numbers may be given in either decimal or exponential form. All of the
following are valid numbers: [5]{.title-ref}, [-3.1]{.title-ref},
[1e6]{.title-ref}, and [3.12e-2]{.title-ref}.

As an extension to the standard Lepton syntax, it is also possible to
use LAMMPS [variables](variable) in the format \"v_name\". Before
evaluating the expression, \"v_name\" will be replaced with the value of
the variable \"name\". This is compatible with all kinds of scalar
variables, but not with vectors, arrays, local, or per-atom variables.
If necessary, a custom scalar variable needs to be defined that can
access the desired (single) item from a non-scalar variable. As an
example, the following lines will instruct LAMMPS to ramp the force
constant for a harmonic bond from 100.0 to 200.0 during the next run:

``` LAMMPS
variable fconst equal ramp(100.0, 200)
bond_style lepton
bond_coeff 1 1.5 "v_fconst * (r^2)"
```

An expression may be followed by definitions for intermediate values
that appear in the expression. A semicolon \";\" is used as a delimiter
between value definitions. For example, the expression:

``` C
a^2+a*b+b^2; a=a1+a2; b=b1+b2
```

is exactly equivalent to

``` C
(a1+a2)^2+(a1+a2)*(b1+b2)+(b1+b2)^2
```

The definition of an intermediate value may itself involve other
intermediate values. Whitespace and quotation characters (\'\'\' and
\'\"\') are ignored. All uses of a value must appear *before* that
value\'s definition. For efficiency reasons, the expression string is
parsed, optimized, and then stored in an internal, pre-parsed
representation for evaluation.

Evaluating a Lepton expression is typically between 2.5 and 5 times
slower than the corresponding compiled and optimized C++ code. If
additional speed or GPU acceleration (via GPU or KOKKOS) is required,
the interaction can be represented as a table. Suitable table files can
be created either internally using the [pair_write](pair_write) or
[bond_write](bond_write) command or through the Python scripts in the
[tools/tabulate](tabulate) folder.

------------------------------------------------------------------------

Styles with a *gpu*, *intel*, *kk*, *omp*, or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the [Accelerator packages](Speed_packages)
page. The accelerated styles take the same arguments and should produce
the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, INTEL, KOKKOS, OPENMP, and
OPT packages, respectively. They are only enabled if LAMMPS was built
with those packages. See the [Build package](Build_package) page for
more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the [-suffix command-line
switch](Run_options) when you invoke LAMMPS, or you can use the
[suffix](suffix) command in your input script.

See the [Accelerator packages](Speed_packages) page for more
instructions on how to use the accelerated styles effectively.

------------------------------------------------------------------------

## Restrictions

This dihedral style is part of the LEPTON package and only enabled if
LAMMPS was built with this package. See the [Build
package](Build_package) page for more info.

## Related commands

[dihedral_coeff](dihedral_coeff), [dihedral_style
table](dihedral_table), [bond_style lepton](bond_lepton), [angle_style
lepton](angle_lepton)

## Default

none
