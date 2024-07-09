# The `lammps` Python module

The LAMMPS Python interface is implemented as a module called
`lammps`{.interpreted-text role="py:mod"} which is defined in the
`lammps` package in the `python` folder of the LAMMPS source code
distribution. After compilation of LAMMPS, the module can be installed
into a Python system folder or a user folder with `make install-python`.
Components of the module can then loaded into a Python session with the
`import` command.

There are multiple Python interface classes in the
`lammps`{.interpreted-text role="py:mod"} module:

-   the `lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_
    class. This is a wrapper around the C-library interface and its
    member functions try to replicate the [C-library API](lammps_c_api)
    closely. This is the most feature-complete Python API.
-   the `PyLammps <lammps.PyLammps>`{.interpreted-text
    role="py:class"}\_\_ class. This is a more high-level and more
    Python style class implemented on top of the
    `lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_
    class.
-   the `IPyLammps <lammps.IPyLammps>`{.interpreted-text
    role="py:class"}\_\_ class is derived from
    `PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
    and adds embedded graphics features to conveniently include LAMMPS
    into [Jupyter](https://jupyter.org/)\_ notebooks.

::: {.admonition .note}
Version check

The `lammps`{.interpreted-text role="py:mod"} module stores the version
number of the LAMMPS version it is installed from. When initializing the
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_ class,
this version is checked to be the same as the result from
`lammps.version`{.interpreted-text role="py:func"}, the version of the
LAMMPS shared library that the module interfaces to. If the they are not
the same an AttributeError exception is raised since a mismatch of
versions (e.g. due to incorrect use of the `LD_LIBRARY_PATH` or
`PYTHONPATH` environment variables can lead to crashes or data
corruption and otherwise incorrect behavior.
:::

::: {.automodule members="" noindex=""}
lammps
:::

------------------------------------------------------------------------

## The `lammps` class API

The `lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_
class is the core of the LAMMPS Python interfaces. It is a wrapper
around the [LAMMPS C library API](lammps_c_api) using the [Python ctypes
module](https://docs.python.org/3/library/ctypes.html)\_ and a shared
library compiled from the LAMMPS sources code. The individual methods in
this class try to closely follow the corresponding C functions. The
handle argument that needs to be passed to the C functions is stored
internally in the class and automatically added when calling the C
library functions. Below is a detailed documentation of the API.

::: {.autoclass members=""}
lammps.lammps
:::

::: {.autoclass members=""}
lammps.numpy_wrapper::numpy_wrapper
:::

------------------------------------------------------------------------

## The `PyLammps` class API

The `PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
class is a wrapper that creates a simpler, more \"Pythonic\" interface
to common LAMMPS functionality. LAMMPS data structures are exposed
through objects and properties. This makes Python scripts shorter and
more concise. See the [PyLammps Tutorial](Howto_pylammps) for an
introduction on how to use this interface.

::: {.autoclass members=""}
lammps.PyLammps
:::

::: {.autoclass members=""}
lammps.AtomList
:::

::: {.autoclass members=""}
lammps.Atom
:::

::: {.autoclass members=""}
lammps.Atom2D
:::

------------------------------------------------------------------------

## The `IPyLammps` class API

The `IPyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
class is an extension of `PyLammps <lammps.PyLammps>`{.interpreted-text
role="py:class"}\_\_, adding additional functions to quickly display
visualizations such as images and videos inside of IPython. See the
[PyLammps Tutorial](Howto_pylammps) for examples.

::: {.autoclass members=""}
lammps.IPyLammps
:::

------------------------------------------------------------------------

## Additional components of the `lammps` module

The `lammps`{.interpreted-text role="py:mod"} module additionally
contains several constants and the
`NeighList <lammps.NeighList>`{.interpreted-text role="py:class"}\_\_
class:

### Data Types {#py_datatype_constants}

### Style Constants {#py_style_constants}

### Type Constants {#py_type_constants}

### Variable Type Constants {#py_vartype_constants}

### Classes representing internal objects

::: {.autoclass members="" no-undoc-members=""}
lammps.NeighList
:::

::: {.autoclass members="" no-undoc-members=""}
lammps.numpy_wrapper::NumPyNeighList
:::
