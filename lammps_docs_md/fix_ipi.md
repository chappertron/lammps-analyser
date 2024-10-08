# fix ipi command

## Syntax

    fix ID group-ID ipi address port [unix] [reset]

-   ID, group-ID are documented in [fix](fix) command

-   ipi = style name of this fix command

-   address = internet address (FQDN or IP), or UNIX socket name

-   port = port number (ignored for UNIX sockets)

-   zero or more keywords may be appended

-   keyword = *unix* or *reset*

        *unix* args = none = use a unix socket
        *reset* args = none = reset electrostatics at each call

## Examples

``` LAMMPS
fix 1 all ipi my.server.com 12345
fix 1 all ipi mysocket 666 unix reset
```

## Description

This fix enables LAMMPS to be run as a client for the i-PI Python
wrapper [(IPI)](ipihome) for performing a path integral molecular
dynamics (PIMD) simulation. The philosophy behind i-PI is described in
the following publication [(IPI-CPC)](IPICPC).

A version of the i-PI package, containing only files needed for use with
LAMMPS, is provided in the tools/i-pi directory. See the
tools/i-pi/manual.pdf for an introduction to i-PI. The
examples/PACKAGES/i-pi directory contains example scripts for using i-PI
with LAMMPS.

In brief, the path integral molecular dynamics is performed by the
Python wrapper, while the client (LAMMPS in this case) simply computes
forces and energy for each configuration. The communication between the
two components takes place using sockets, and is reduced to the bare
minimum. All the parameters of the dynamics are specified in the input
of i-PI, and all the parameters of the force field must be specified as
LAMMPS inputs, preceding the *fix ipi* command.

The server address must be specified by the *address* argument, and can
be either the IP address, the fully-qualified name of the server, or the
name of a UNIX socket for local, faster communication. In the case of
internet sockets, the *port* argument specifies the port number on which
i-PI is listening, while the *unix* optional switch specifies that the
socket is a UNIX socket.

Note that there is no check of data integrity, or that the atomic
configurations make sense. It is assumed that the species in the i-PI
input are listed in the same order as in the data file of LAMMPS. The
initial configuration is ignored, as it will be substituted with the
coordinates received from i-PI before forces are ever evaluated.

A note of caution when using potentials that contain long-range
electrostatics, or that contain parameters that depend on box size: all
of these options will be initialized based on the cell size in the
LAMMPS-side initial configuration and kept constant during the run. This
is required to e.g. obtain reproducible and conserved forces. If the
cell varies too wildly, it may be advisable to re-initialize these
interactions at each call. This behavior can be requested by setting the
*reset* switch.

## Restart, fix_modify, output, run start/stop, minimize info

There is no restart information associated with this fix, since all the
dynamical parameters are dealt with by i-PI.

## Restrictions

Using this fix on anything other than all atoms requires particular
care, since i-PI will know nothing on atoms that are not those whose
coordinates are transferred. However, one could use this strategy to
define an external potential acting on the atoms that are moved by i-PI.

Since the i-PI code uses atomic units internally, this fix needs to
convert LAMMPS data to and from its [specified units](units) accordingly
when communicating with i-PI. This is not possible for reduced units
(\"units lj\") and thus *fix ipi* will stop with an error in this case.

This fix is part of the MISC package. It is only enabled if LAMMPS was
built with that package. See the [Build package](Build_package) page for
more info. Because of the use of UNIX domain sockets, this fix will only
work in a UNIX environment.

## Related commands

[fix nve](fix_nve)

------------------------------------------------------------------------

::: {#IPICPC}
**(IPI-CPC)** Ceriotti, More and Manolopoulos, Comp Phys Comm, 185,
1019-1026 (2014).
:::

::: {#ipihome}
**(IPI)** <https://ipi-code.org>\_
:::
