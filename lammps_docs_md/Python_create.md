# Creating or deleting a LAMMPS object {#python_create_lammps}

With the Python interface the creation of a
`LAMMPS <LAMMPS_NS::LAMMPS>`{.interpreted-text role="cpp:class"}\_\_
instance is included in the constructors for the
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_,
`PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_, and
`IPyLammps <lammps.IPyLammps>`{.interpreted-text role="py:class"}\_\_
classes. Internally it will call either `lammps_open`{.interpreted-text
role="cpp:func"} or `lammps_open_no_mpi`{.interpreted-text
role="cpp:func"} from the C library API to create the class instance.

All arguments are optional. The *name* argument allows loading a LAMMPS
shared library that is named `liblammps_machine.so` instead of the
default name of `liblammps.so`. In most cases the latter will be
installed or used. The *ptr* argument is for use of the
`lammps`{.interpreted-text role="py:mod"} module from inside a LAMMPS
instance, e.g. with the [python](python) command, where a pointer to the
already existing `LAMMPS <LAMMPS_NS::LAMMPS>`{.interpreted-text
role="cpp:class"}\_\_ class instance can be passed to the Python class
and used instead of creating a new instance. The *comm* argument may be
used in combination with the [mpi4py](https://mpi4py.readthedocs.io/)\_
module to pass an MPI communicator to LAMMPS and thus it is possible to
run the Python module like the library interface on a subset of the MPI
ranks after splitting the communicator.

Here are simple examples using all three Python interfaces:

:::::: tabs
::: tab
lammps API

``` python
from lammps import lammps

# NOTE: argv[0] is set by the lammps class constructor
args = ["-log", "none"]

# create LAMMPS instance
lmp = lammps(cmdargs=args)

# get and print numerical version code
print("LAMMPS Version: ", lmp.version())

# explicitly close and delete LAMMPS instance (optional)
lmp.close()
```
:::

::: tab
PyLammps API

The `PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
class is a wrapper around the `lammps <lammps.lammps>`{.interpreted-text
role="py:class"}\_\_ class and all of its lower level functions. By
default, it will create a new instance of
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_ passing
along all arguments to the constructor of
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_.

``` python
from lammps import PyLammps

# NOTE: argv[0] is set by the lammps class constructor
args = ["-log", "none"]

# create LAMMPS instance
L = PyLammps(cmdargs=args)

# get and print numerical version code
print("LAMMPS Version: ", L.version())

# explicitly close and delete LAMMPS instance (optional)
L.close()
```

`PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
objects can also be created on top of an existing
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_ object:

``` python
from lammps import lammps, PyLammps
...
# create LAMMPS instance
lmp = lammps(cmdargs=args)

# create PyLammps instance using previously created LAMMPS instance
L = PyLammps(ptr=lmp)
```

This is useful if you have to create the
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_ instance
is a specific way, but want to take advantage of the
`PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
interface.
:::

::: tab
IPyLammps API

The `IPyLammps <lammps.IPyLammps>`{.interpreted-text
role="py:class"}\_\_ class is an extension of the
`PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
class. It has the same construction behavior. By default, it will create
a new instance of `lammps`{.interpreted-text role="py:class"} passing
along all arguments to the constructor of `lammps`{.interpreted-text
role="py:class"}.

``` python
from lammps import IPyLammps

# NOTE: argv[0] is set by the lammps class constructor
args = ["-log", "none"]

# create LAMMPS instance
L = IPyLammps(cmdargs=args)

# get and print numerical version code
print("LAMMPS Version: ", L.version())

# explicitly close and delete LAMMPS instance (optional)
L.close()
```

You can also initialize IPyLammps on top of an existing
`lammps`{.interpreted-text role="py:class"} or
`PyLammps`{.interpreted-text role="py:class"} object:

``` python
from lammps import lammps, IPyLammps
...
# create LAMMPS instance
lmp = lammps(cmdargs=args)

# create PyLammps instance using previously created LAMMPS instance
L = PyLammps(ptr=lmp)
```

This is useful if you have to create the
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_ instance
is a specific way, but want to take advantage of the
`IPyLammps <lammps.IPyLammps>`{.interpreted-text role="py:class"}\_\_
interface.
:::
::::::

In all of the above cases, same as with the [C library
API](lammps_c_api), this will use the `MPI_COMM_WORLD` communicator for
the MPI library that LAMMPS was compiled with.

The `lmp.close() <lammps.lammps.close()>`{.interpreted-text
role="py:func"}\_\_ call is optional since the LAMMPS class instance
will also be deleted automatically during the
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_ class
destructor. Instead of
`lmp.close() <lammps.lammps.close()>`{.interpreted-text
role="py:func"}\_\_ it is also possible to call
`lmp.finalize() <lammps.lammps.finalize()>`{.interpreted-text
role="py:func"}\_\_; this will destruct the LAMMPS instance, but also
finalized and release the MPI and/or Kokkos environment if enabled and
active.

Note that you can create multiple LAMMPS objects in your Python script,
and coordinate and run multiple simulations, e.g.

``` python
from lammps import lammps
lmp1 = lammps()
lmp2 = lammps()
lmp1.file("in.file1")
lmp2.file("in.file2")
```
