# Download an executable for Linux

Binaries are available for different versions of Linux:

-   [Pre-built static Linux x86_64 executables](static)
-   [Pre-built Ubuntu and Debian Linux executables](ubuntu)
-   [Pre-built Fedora Linux executables](fedora)
-   [Pre-built EPEL Linux executables (RHEL, CentOS)](epel)
-   [Pre-built OpenSuse Linux executables](opensuse)
-   [Gentoo Linux executable](gentoo)
-   [Arch Linux build-script](arch)

:::: note
::: title
Note
:::

If you have questions about these pre-compiled LAMMPS executables, you
need to contact the people preparing those packages. The LAMMPS
developers have no control over how they configure and build their
packages and when they update them. They may only provide packages for
stable release versions and not always update the packages in a timely
fashion after a new LAMMPS release is made.
::::

------------------------------------------------------------------------

## Pre-built static Linux x86_64 executables {#static}

Pre-built LAMMPS executables for Linux, that are statically linked and
compiled for 64-bit x86 CPUs (x86_64 or AMD64) are available for
download at <https://download.lammps.org/static/>\_. Because of that
static linkage (and unlike the Linux distribution specific packages
listed below), they do not depend on any installed software and thus
should run on *any* 64-bit x86 machine with *any* Linux version.

These executable include most of the available packages and multi-thread
parallelization (via INTEL, KOKKOS, or OPENMP package). They are **not**
compatible with MPI. Several of the LAMMPS tools executables (e.g.
`msi2lmp`) and the `lammps-shell` program are included as well. Because
of the static linkage, there is no `liblammps.so` library file and thus
also the LAMMPS python module, which depends on it, is not included.

The compressed tar archives available for download have names following
the pattern `lammps-linux-x86_64-<version>.tar.gz` and will all unpack
into a `lammps-static` folder. The executables are then in the
`lammps-static/bin/` folder. Since they do not depend on any other
software, they may be freely moved or copied around.

------------------------------------------------------------------------

## Pre-built Ubuntu and Debian Linux executables {#ubuntu}

A pre-built LAMMPS executable, suitable for running on the latest Ubuntu
and Debian Linux versions, can be downloaded as a Debian package. This
allows you to install LAMMPS with a single command, and stay (mostly)
up-to-date with the current stable version of LAMMPS by simply updating
your operating system.

To install LAMMPS do the following once:

``` bash
sudo apt-get install lammps
```

This downloads an executable named `lmp` to your box and multiple
packages with supporting data, examples and libraries as well as any
missing dependencies. For example, the LAMMPS binary in this package is
built with the [KIM package](kim) enabled, which results in the above
command also installing the `kim-api` binaries when LAMMPS is installed,
unless they were installed already. In order to use potentials from
[openkim.org](https://openkim.org)\_, you can also install the
`openkim-models` package:

``` bash
sudo apt-get install openkim-models
```

Or use the [KIM-API
commands](https://openkim.org/doc/usage/obtaining-models/#installing_api)\_
to download and install individual models.

This LAMMPS executable can then be used in the usual way to run input
scripts:

``` bash
lmp -in in.lj
```

To update LAMMPS to the latest packaged version, do the following:

``` bash
sudo apt-get update
```

This will also update other packages on your system.

To uninstall LAMMPS, do the following:

``` bash
sudo apt-get remove lammps
```

Please use `lmp -help` to see which compilation options, packages, and
styles are included in the binary.

Thanks to Anton Gladky (gladky.anton at gmail.com) for setting up this
Ubuntu package capability.

------------------------------------------------------------------------

## Pre-built Fedora Linux executables {#fedora}

Pre-built [LAMMPS packages for stable
releases](https://packages.fedoraproject.org/pkgs/lammps/)\_ are
available in the Fedora Linux distribution since Fedora version 28. The
packages can be installed via the dnf package manager. There are 3 basic
varieties (lammps = no MPI, lammps-mpich = MPICH MPI library,
lammps-openmpi = OpenMPI MPI library) and for each support for linking
to the C library interface (lammps-devel, lammps-mpich-devel,
lammps-openmpi-devel), the header for compiling programs using the C
library interface (lammps-headers), and the LAMMPS python module for
Python 3. All packages can be installed at the same time and the name of
the LAMMPS executable is `lmp` and `lmp_openmpi` or `lmp_mpich`
respectively. By default, `lmp` will refer to the serial executable,
unless one of the MPI environment modules is loaded
(`module load mpi/mpich-x86_64` or `module load mpi/openmpi-x86_64`).
Then the corresponding parallel LAMMPS executable can be used. The same
mechanism applies when loading the LAMMPS python module.

To install LAMMPS with OpenMPI and run an input `in.lj` with 2 CPUs do:

``` bash
dnf install lammps-openmpi
module load mpi/openmpi-x86_64
mpirun -np 2 lmp -in in.lj
```

The `dnf install` command is needed only once. In case of a new LAMMPS
stable release, `dnf update` will automatically update to the newer
version as soon as the RPM files are built and uploaded to the download
mirrors. The `module load` command is needed once per (shell) session or
shell terminal instance, unless it is automatically loaded from the
shell profile.

The LAMMPS binary is built with the [KIM package](kim) which results in
the above command also installing the [kim-api]{.title-ref} binaries
when LAMMPS is installed. In order to use potentials from
[openkim.org](https://openkim.org)\_, you can install the
[openkim-models]{.title-ref} package

``` bash
dnf install openkim-models
```

Please use `lmp -help` to see which compilation options, packages, and
styles are included in the binary.

Thanks to Christoph Junghans (LANL) for making LAMMPS available in
Fedora.

------------------------------------------------------------------------

## Pre-built EPEL Linux executable {#epel}

Pre-built LAMMPS (and KIM) packages for stable releases are available in
the [Extra Packages for Enterprise Linux (EPEL)
repository](https://docs.fedoraproject.org/en-US/epel/)\_ for use with
Red Hat Enterprise Linux (RHEL) or CentOS version 7.x and compatible
Linux distributions. Names of packages, executable, and content are the
same as described above for Fedora Linux. But RHEL/CentOS 7.x uses the
`yum` package manager instead of `dnf` in Fedora 28.

Please use `lmp -help` to see which compilation options, packages, and
styles are included in the binary.

Thanks to Christoph Junghans (LANL) for making LAMMPS available in EPEL.

------------------------------------------------------------------------

## Pre-built OpenSuse Linux executable {#opensuse}

A pre-built LAMMPS package for stable releases is available in OpenSuse
as of Leap 15.0. You can install the package with:

``` bash
zypper install lammps
```

This includes support for OpenMPI. The name of the LAMMPS executable is
`lmp`. To run an input in parallel on 2 CPUs you would do:

``` bash
mpirun -np 2 lmp -in in.lj
```

Please use `lmp -help` to see which compilation options, packages, and
styles are included in the binary.

The LAMMPS binary is built with the [KIM package](kim) which results in
the above command also installing the [kim-api]{.title-ref} binaries
when LAMMPS is installed. In order to use potentials from
[openkim.org](https://openkim.org)\_, you can install the
[openkim-models]{.title-ref} package

``` bash
zypper install openkim-models
```

Thanks to Christoph Junghans (LANL) for making LAMMPS available in
OpenSuse.

------------------------------------------------------------------------

## Gentoo Linux executable {#gentoo}

LAMMPS is part of [Gentoo\'s main package
tree](https://packages.gentoo.org/packages/sci-physics/lammps)\_ and can
be installed by typing:

``` bash
emerge --ask lammps
```

Note that in Gentoo the LAMMPS source code is downloaded and the package
is then compiled and installed on your machine.

Certain LAMMPS packages can be enabled via USE flags, type

``` bash
equery uses lammps
```

for details.

Thanks to Nicolas Bock and Christoph Junghans (LANL) for setting up this
Gentoo capability.

------------------------------------------------------------------------

## Archlinux build-script {#arch}

LAMMPS is available via Arch\'s unofficial Arch User repository (AUR).
There are three scripts available, named
[lammps](https://aur.archlinux.org/packages/lammps)\_,
[lammps-beta](https://aur.archlinux.org/packages/lammps)\_ and
[lammps-git](https://aur.archlinux.org/packages/lammps)\_. They
respectively package the stable, feature, and git releases.

To install, you will need to have the git package installed. You may use
any of the above names in-place of lammps.

``` bash
git clone https://aur.archlinux.org/lammps.git
cd lammps
makepkg -s
makepkg -i
```

To update LAMMPS, you may repeat the above, or change into the cloned
directory, and execute the following, after which, if there are any
changes, you may use makepkg as above.

``` bash
git pull
```

Alternatively, you may use an AUR helper to install these packages.

Note that the AUR provides build-scripts that download the source code
and then build and install the package on your machine.
