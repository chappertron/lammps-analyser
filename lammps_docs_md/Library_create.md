# Creating or deleting a LAMMPS object

This section documents the following functions:

-   `lammps_open`{.interpreted-text role="cpp:func"}
-   `lammps_open_no_mpi`{.interpreted-text role="cpp:func"}
-   `lammps_open_fortran`{.interpreted-text role="cpp:func"}
-   `lammps_close`{.interpreted-text role="cpp:func"}
-   `lammps_mpi_init`{.interpreted-text role="cpp:func"}
-   `lammps_mpi_finalize`{.interpreted-text role="cpp:func"}
-   `lammps_kokkos_finalize`{.interpreted-text role="cpp:func"}
-   `lammps_python_finalize`{.interpreted-text role="cpp:func"}
-   `lammps_error`{.interpreted-text role="cpp:func"}

------------------------------------------------------------------------

The `lammps_open`{.interpreted-text role="cpp:func"} and
`lammps_open_no_mpi`{.interpreted-text role="cpp:func"} functions are
used to create and initialize a `LAMMPS`{.interpreted-text
role="cpp:func"} instance. They return a reference to this instance as a
`void *` pointer to be used as the \"handle\" argument in subsequent
function calls until that instance is destroyed by calling
`lammps_close`{.interpreted-text role="cpp:func"}. Here is a simple
example demonstrating its use:

``` c
#include "library.h"
#include <stdio.h>

int main(int argc, char **argv)
{
  void *handle;
  int version;
  const char *lmpargv[] = { "liblammps", "-log", "none"};
  int lmpargc = sizeof(lmpargv)/sizeof(const char *);

  /* create LAMMPS instance */
  handle = lammps_open_no_mpi(lmpargc, (char **)lmpargv, NULL);
  if (handle == NULL) {
    printf("LAMMPS initialization failed");
    lammps_mpi_finalize();
    return 1;
  }

  /* get and print numerical version code */
  version = lammps_version(handle);
  printf("LAMMPS Version: %d\n",version);

  /* delete LAMMPS instance and shut down MPI */
  lammps_close(handle);
  lammps_mpi_finalize();
  return 0;
}
```

The LAMMPS library uses the MPI library it was compiled with and will
either run on all processors in the `MPI_COMM_WORLD` communicator or on
the set of processors in the communicator passed as the `comm` argument
of `lammps_open`{.interpreted-text role="cpp:func"}. This means the
calling code can run LAMMPS on all or a subset of processors. For
example, a wrapper code might decide to alternate between LAMMPS and
another code, allowing them both to run on all the processors. Or it
might allocate part of the processors to LAMMPS and the rest to the
other code by creating a custom communicator with `MPI_Comm_split()` and
running both codes concurrently before syncing them up periodically. Or
it might instantiate multiple instances of LAMMPS to perform different
calculations and either alternate between them, run them concurrently on
split communicators, or run them one after the other. The
`lammps_open`{.interpreted-text role="cpp:func"} function may be called
multiple times for this latter purpose.

The `lammps_close`{.interpreted-text role="cpp:func"} function is used
to shut down the `LAMMPS <LAMMPS_NS::LAMMPS>`{.interpreted-text
role="cpp:class"}\_\_ class pointed to by the handle passed as an
argument and free all its memory. This has to be called for every
instance created with one of the `lammps_open`{.interpreted-text
role="cpp:func"} functions. It will, however, **not** call
`MPI_Finalize()`, since that may only be called once. See
`lammps_mpi_finalize`{.interpreted-text role="cpp:func"} for an
alternative to invoking `MPI_Finalize()` explicitly from the calling
program.

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_open
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_open_no_mpi
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_open_fortran
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_close
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_mpi_init
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_mpi_finalize
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_kokkos_finalize
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_python_finalize
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_error
:::
