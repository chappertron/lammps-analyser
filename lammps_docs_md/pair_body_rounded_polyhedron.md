# pair_style body/rounded/polyhedron command

## Syntax

``` LAMMPS
pair_style body/rounded/polyhedron c_n c_t mu delta_ua cutoff
```

    c_n = normal damping coefficient
    c_t = tangential damping coefficient
    mu = normal friction coefficient during gross sliding
    delta_ua = multiple contact scaling factor
    cutoff = global separation cutoff for interactions (distance units), see below for definition

## Examples

``` LAMMPS
pair_style body/rounded/polyhedron 20.0 5.0 0.0 1.0 0.5
pair_coeff * * 100.0 1.0
pair_coeff 1 1 100.0 1.0
```

## Description

Style *body/rounded/polygon* is for use with 3d models of body particles
of style *rounded/polyhedron*. It calculates pairwise body/body
interactions which can include body particles modeled as 1-vertex
spheres with a specified diameter. See the [Howto body](Howto_body) page
for more details on using body rounded/polyhedron particles.

This pairwise interaction between the rounded polyhedra is described in
[Wang](pair-Wang), where a polyhedron does not have sharp corners and
edges, but is rounded at its vertices and edges by spheres centered on
each vertex with a specified diameter. The edges if the polyhedron are
defined between pairs of adjacent vertices. Its faces are defined by a
loop of edges. The sphere diameter for each polygon is specified in the
data file read by the [read data](read_data) command. This is a discrete
element model (DEM) which allows for multiple contact points.

Note that when two particles interact, the effective surface of each
polyhedron particle is displaced outward from each of its vertices,
edges, and faces by half its sphere diameter. The interaction forces and
energies between two particles are defined with respect to the
separation of their respective rounded surfaces, not by the separation
of the vertices, edges, and faces themselves.

This means that the specified cutoff in the pair_style command is the
cutoff distance, $r_c$, for the surface separation, $\delta_n$ (see
figure below). This is the distance at which two particles no longer
interact. If $r_c$ is specified as 0.0, then it is a contact-only
interaction. I.e. the two particles must overlap in order to exert a
repulsive force on each other. If $r_c > 0.0$, then the force between
two particles will be attractive for surface separations from 0 to
$r_c$, and repulsive once the particles overlap.

Note that unlike for other pair styles, the specified cutoff is not the
distance between the centers of two particles at which they stop
interacting. This center-to-center distance depends on the shape and
size of the two particles and their relative orientation. LAMMPS takes
that into account when computing the surface separation distance and
applying the $r_c$ cutoff.

The forces between vertex-vertex, vertex-edge, vertex-face, edge-edge,
and edge-face overlaps are given by:

$$\begin{aligned}
F_n &= \begin{cases}
       k_n \delta_n - c_n v_n,    & \delta_n \le 0 \\
      -k_{na} \delta_n - c_n v_n  & 0 < \delta_n \le r_c \\
       0                          & \delta_n > r_c \\
       \end{cases} \\
F_t &= \begin{cases}
       \mu k_n \delta_n - c_t v_t & \delta_n \le 0 \\
       0                          & \delta_n > 0
       \end{cases}
\end{aligned}$$

![image](JPG/pair_body_rounded.jpg){.align-center}

In [Wang](pair-Wang), the tangential friction force between two
particles that are in contact is modeled differently prior to gross
sliding (i.e. static friction) and during gross-sliding (kinetic
friction). The latter takes place when the tangential deformation
exceeds the Coulomb frictional limit. In the current implementation,
however, we do not take into account frictional history, i.e. we do not
keep track of how many time steps the two particles have been in contact
nor calculate the tangential deformation. Instead, we assume that gross
sliding takes place as soon as two particles are in contact.

The following coefficients must be defined for each pair of atom types
via the [pair_coeff](pair_coeff) command as in the examples above, or in
the data file read by the [read_data](read_data) command:

-   $k_n$ (energy/distance\^2 units)
-   $k_{na}$ (energy/distance\^2 units)

Effectively, $k_n$ and $k_{na}$ are the slopes of the red lines in the
plot above for force versus surface separation, for $\delta_n$ \< 0 and
$0 < \delta_n < r_c$ respectively.

## Mixing, shift, table, tail correction, restart, rRESPA info

This pair style does not support the [pair_modify](pair_modify) mix,
shift, table, and tail options.

This pair style does not write its information to [binary restart
files](restart). Thus, you need to re-specify the pair_style and
pair_coeff commands in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
[run_style respa](run_style) command. It does not support the *inner*,
*middle*, *outer* keywords.

## Restrictions

These pair styles are part of the BODY package. They are only enabled if
LAMMPS was built with that package. See the [Build
package](Build_package) page for more info.

This pair style requires the [newton](newton) setting to be \"on\" for
pair interactions.

## Related commands

[pair_coeff](pair_coeff)

## Default

none

::: {#pair-Wang}
**(Wang)** J. Wang, H. S. Yu, P. A. Langston, F. Y. Fraige, Granular
Matter, 13, 1 (2011).
:::
