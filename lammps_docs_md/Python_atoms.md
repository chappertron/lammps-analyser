# Per-atom properties

Similar to what is described in [Library_atoms]{.title-ref}, the
instances of `lammps <lammps.lammps>`{.interpreted-text
role="py:class"}\_\_, `PyLammps <lammps.PyLammps>`{.interpreted-text
role="py:class"}\_\_, or
`IPyLammps <lammps.IPyLammps>`{.interpreted-text role="py:class"}\_\_
can be used to extract atom quantities and modify some of them. The main
difference between the interfaces is how the information is exposed.

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

nlocal = lmp.extract_global("nlocal")
x = lmp.extract_atom("x")

for i in range(nlocal):
   print("(x,y,z) = (", x[i][0], x[i][1], x[i][2], ")")

lmp.close()
```

**Methods**:

-   `extract_atom() <lammps.lammps.extract_atom()>`{.interpreted-text
    role="py:meth"}\_\_: extract a per-atom quantity

**Numpy Methods**:

-   `numpy.extract_atom() <lammps.numpy_wrapper.numpy_wrapper.extract_atom()>`{.interpreted-text
    role="py:meth"}\_\_: extract a per-atom quantity as numpy array
:::

::: tab
PyLammps/IPyLammps API

All atoms in the current simulation can be accessed by using the
`PyLammps.atoms <lammps.PyLammps.atoms>`{.interpreted-text
role="py:attr"}\_\_ property. Each element of this list is a
`Atom <lammps.Atom>`{.interpreted-text role="py:class"}\_\_ or
`Atom2D <lammps.Atom2D>`{.interpreted-text role="py:class"}\_\_ object.
The attributes of these objects provide access to their data (id, type,
position, velocity, force, etc.):

``` python
# access first atom
L.atoms[0].id
L.atoms[0].type

# access second atom
L.atoms[1].position
L.atoms[1].velocity
L.atoms[1].force
```

Some attributes can be changed:

``` python
# set position in 2D simulation
L.atoms[0].position = (1.0, 0.0)

# set position in 3D simulation
L.atoms[0].position = (1.0, 0.0, 1.0)
```
:::
:::::
