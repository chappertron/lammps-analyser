# System properties

Similar to what is described in [Library_properties]{.title-ref}, the
instances of `lammps <lammps.lammps>`{.interpreted-text
role="py:class"}\_\_, `PyLammps <lammps.PyLammps>`{.interpreted-text
role="py:class"}\_\_, or
`IPyLammps <lammps.IPyLammps>`{.interpreted-text role="py:class"}\_\_
can be used to extract different kinds of information about the active
LAMMPS instance and also to modify some of it. The main difference
between the interfaces is how the information is exposed.

While the `lammps <lammps.lammps>`{.interpreted-text
role="py:class"}\_\_ is just a thin layer that wraps C API calls,
`PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_ and
`IPyLammps <lammps.IPyLammps>`{.interpreted-text role="py:class"}\_\_
expose information as objects and properties.

In some cases the data returned is a direct reference to the original
data inside LAMMPS cast to `ctypes` pointers. Where possible, the
wrappers will determine the `ctypes` data type and cast pointers
accordingly. If `numpy` is installed arrays can also be extracted as
numpy arrays, which will access the C arrays directly and have the
correct dimensions to protect against invalid accesses.

:::: warning
::: title
Warning
:::

When accessing per-atom data, please note that this data is the
per-processor local data and indexed accordingly. These arrays can
change sizes and order at every neighbor list rebuild and atom sort
event as atoms are migrating between subdomains.
::::

::::: tabs
::: tab
lammps API

``` python
from lammps import lammps

lmp = lammps()
lmp.file("in.sysinit")

natoms = lmp.get_natoms()
print(f"running simulation with {natoms} atoms")

lmp.command("run 1000 post no");

for i in range(10):
   lmp.command("run 100 pre no post no")
   pe = lmp.get_thermo("pe")
   ke = lmp.get_thermo("ke")
   print(f"PE = {pe}\nKE = {ke}")

lmp.close()
```

**Methods**:

-   `version() <lammps.lammps.version()>`{.interpreted-text
    role="py:meth"}\_\_: return the numerical version id, e.g. LAMMPS 2
    Sep 2015 -\> 20150902
-   `get_thermo() <lammps.lammps.get_thermo()>`{.interpreted-text
    role="py:meth"}\_\_: return current value of a thermo keyword
-   `last_thermo() <lammps.lammps.last_thermo()>`{.interpreted-text
    role="py:meth"}\_\_: return a dictionary of the last thermodynamic
    output
-   `get_natoms() <lammps.lammps.get_natoms()>`{.interpreted-text
    role="py:meth"}\_\_: total \# of atoms as int
-   `reset_box() <lammps.lammps.reset_box()>`{.interpreted-text
    role="py:meth"}\_\_: reset the simulation box size
-   `extract_setting() <lammps.lammps.extract_setting()>`{.interpreted-text
    role="py:meth"}\_\_: return a global setting
-   `extract_global() <lammps.lammps.extract_global()>`{.interpreted-text
    role="py:meth"}\_\_: extract a global quantity
-   `extract_box() <lammps.lammps.extract_box()>`{.interpreted-text
    role="py:meth"}\_\_: extract box info
-   `create_atoms() <lammps.lammps.create_atoms()>`{.interpreted-text
    role="py:meth"}\_\_: create N atoms with IDs, types, x, v, and image
    flags

**Properties**:

-   `last_thermo_step <lammps.lammps.last_thermo_step>`{.interpreted-text
    role="py:attr"}\_\_: the last timestep thermodynamic output was
    computed
:::

::: tab
PyLammps/IPyLammps API

In addition to the functions provided by
`lammps <lammps.lammps>`{.interpreted-text role="py:class"}\_\_,
`PyLammps <lammps.PyLammps>`{.interpreted-text role="py:class"}\_\_
objects have several properties which allow you to query the system
state:

L.system

:   Is a dictionary describing the system such as the bounding box or
    number of atoms

L.system.xlo, L.system.xhi

:   bounding box limits along x-axis

L.system.ylo, L.system.yhi

:   bounding box limits along y-axis

L.system.zlo, L.system.zhi

:   bounding box limits along z-axis

L.communication

:   configuration of communication subsystem, such as the number of
    threads or processors

L.communication.nthreads

:   number of threads used by each LAMMPS process

L.communication.nprocs

:   number of MPI processes used by LAMMPS

L.fixes

:   List of fixes in the current system

L.computes

:   List of active computes in the current system

L.dump

:   List of active dumps in the current system

L.groups

:   List of groups present in the current system

**Retrieving the value of an arbitrary LAMMPS expressions**

LAMMPS expressions can be immediately evaluated by using the `eval`
method. The passed string parameter can be any expression containing
global [thermo]{.title-ref} values, variables, compute or fix data (see
[Howto_output]{.title-ref}):

``` python
result = L.eval("ke") # kinetic energy
result = L.eval("pe") # potential energy

result = L.eval("v_t/2.0")
```

**Example**

``` python
from lammps import PyLammps

L = PyLammps()
L.file("in.sysinit")

print(f"running simulation with {L.system.natoms} atoms")

L.run(1000, "post no");

for i in range(10):
   L.run(100, "pre no post no")
   pe = L.eval("pe")
   ke = L.eval("ke")
   print(f"PE = {pe}\nKE = {ke}")
```
:::
:::::
