# Scatter/gather operations

This section has functions which gather per-atom data from one or more
processors into a contiguous global list ordered by atom ID. The same
list is returned to all calling processors. It also contains functions
which scatter per-atom data from a contiguous global list across the
processors that own those atom IDs. It also has a create_atoms()
function which can create new atoms by scattering them appropriately to
owning processors in the LAMMPS spatial decomposition.

It documents the following functions:

-   `lammps_gather_atoms`{.interpreted-text role="cpp:func"}
-   `lammps_gather_atoms_concat`{.interpreted-text role="cpp:func"}
-   `lammps_gather_atoms_subset`{.interpreted-text role="cpp:func"}
-   `lammps_scatter_atoms`{.interpreted-text role="cpp:func"}
-   `lammps_scatter_atoms_subset`{.interpreted-text role="cpp:func"}
-   `lammps_gather_bonds`{.interpreted-text role="cpp:func"}
-   `lammps_gather_angles`{.interpreted-text role="cpp:func"}
-   `lammps_gather_dihedrals`{.interpreted-text role="cpp:func"}
-   `lammps_gather_impropers`{.interpreted-text role="cpp:func"}
-   `lammps_gather`{.interpreted-text role="cpp:func"}
-   `lammps_gather_concat`{.interpreted-text role="cpp:func"}
-   `lammps_gather_subset`{.interpreted-text role="cpp:func"}
-   `lammps_scatter`{.interpreted-text role="cpp:func"}
-   `lammps_scatter_subset`{.interpreted-text role="cpp:func"}
-   `lammps_create_atoms`{.interpreted-text role="cpp:func"}

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_atoms
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_atoms_concat
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_atoms_subset
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_scatter_atoms
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_scatter_atoms_subset
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_bonds
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_angles
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_dihedrals
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_impropers
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_concat
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_gather_subset
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_scatter
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_scatter_subset
:::

------------------------------------------------------------------------

::: {.doxygenfunction project="progguide"}
lammps_create_atoms(void *handle, int n, const int*id, const int *type,
const double*x, const double *v, const int*image, int bexpand)
:::
