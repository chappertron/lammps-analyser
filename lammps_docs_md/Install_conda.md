# Download an executable for Linux or macOS via Conda

Pre-compiled LAMMPS binaries are available for macOS and Linux via the
[Conda](https://docs.conda.io/en/latest/index.html)\_ package management
system.

First, one must set up the Conda package manager on your system. Follow
the instructions to install
[Miniconda](https://docs.conda.io/en/latest/miniconda.html)\_, then
create a conda environment (named [my-lammps-env]{.title-ref} or
whatever you prefer) for your LAMMPS install:

``` bash
conda config --add channels conda-forge
conda create -n my-lammps-env
```

Then, you can install LAMMPS on your system with the following command:

``` bash
conda activate my-lammps-env
conda install lammps
```

The LAMMPS binary is built with the [KIM package](kim), which results in
Conda also installing the [kim-api]{.title-ref} binaries when LAMMPS is
installed. In order to use potentials from
[openkim.org](https://openkim.org)\_, you can install the
[openkim-models]{.title-ref} package

``` bash
conda install openkim-models
```

If you have problems with the installation, you can post issues to [this
link](https://github.com/conda-forge/lammps-feedstock/issues)\_. Thanks
to Jan Janssen (Max-Planck-Institut fuer Eisenforschung) for setting up
the Conda capability.

:::: note
::: title
Note
:::

If you have questions about these pre-compiled LAMMPS executables, you
need to contact the people preparing those packages. The LAMMPS
developers have no control over their choices of how they configure and
build their packages and when they update them.
::::
