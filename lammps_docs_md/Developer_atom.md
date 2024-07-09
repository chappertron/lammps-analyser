# Accessing per-atom data

This page discusses how per-atom data is managed in LAMMPS, how it can
be accessed, what communication patters apply, and some of the utility
functions that exist for a variety of purposes.

## Owned and ghost atoms

As described on the [parallel partitioning
algorithms](Developer_par_part) page, LAMMPS uses a domain decomposition
of the simulation domain, either in a *brick* or *tiled* manner. Each
MPI process *owns* exactly one subdomain and the atoms within it. To
compute forces for tuples of atoms that are spread across sub-domain
boundaries, also a \"halo\" of *ghost* atoms are maintained within a the
communication cutoff distance of its subdomain.

The total number of atoms is stored in [Atom::natoms]{.title-ref}
(within any typical class this can be referred to at
[atom-\>natoms]{.title-ref}. The number of *owned* (or \"local\" atoms)
are stored in [Atom::nlocal]{.title-ref}; the number of *ghost* atoms is
stored in [Atom::nghost]{.title-ref}. The sum of
[Atom::nlocal]{.title-ref} over all MPI processes should be
[Atom::natoms]{.title-ref}. This is by default regularly checked by the
Thermo class, and if the sum does not match, LAMMPS stops with a \"lost
atoms\" error. For convenience also the property
[Atom::nmax]{.title-ref} is available, this is the maximum of
[Atom::nlocal + Atom::nghost]{.title-ref} across all MPI processes.

Per-atom properties are either managed by the atom style, or individual
classes. or as custom arrays by the individual classes. If only access
to *owned* atoms is needed, they are usually allocated to be of size
[Atom::nlocal]{.title-ref}, otherwise of size [Atom::nmax]{.title-ref}.
Please note that not all per-atom properties are available or updated on
*ghost* atoms. For example, per-atom velocities are only updated with
[comm_modify vel yes](comm_modify).

## Atom indexing

When referring to individual atoms, they may be indexed by their local
*index*, their index in their [Atom::x]{.title-ref} array. This is
densely populated containing first all *owned* atoms (index \<
[Atom::nlocal]{.title-ref}) and then all *ghost* atoms. The order of
atoms in these arrays can change due to atoms migrating between between
subdomains, atoms being added or deleted, or atoms being sorted for
better cache efficiency. Atoms are globally uniquely identified by their
*atom ID*. There may be multiple atoms with the same atom ID present,
but only one of them may be an *owned* atom.

To find the local *index* of an atom, when the *atom ID* is known, the
[Atom::map()]{.title-ref} function may be used. It will return the local
atom index or -1. If the returned value is between 0 (inclusive) and
[Atom::nlocal]{.title-ref} (exclusive) it is an *owned* or \"local\"
atom; for larger values the atom is present as a ghost atom; for a value
of -1, the atom is not present on the current subdomain at all.

If multiple atoms with the same tag exist in the same subdomain, they
can be found via the [Atom::sametag]{.title-ref} array. It points to the
next atom index with the same tag or -1 if there are no more atoms with
the same tag. The list will be exhaustive when starting with an index of
an *owned* atom, since the atom IDs are unique, so there can only be one
such atom. Example code to count atoms with same atom ID in subdomain:

``` c++
for (int i = 0; i < atom->nlocal; ++i) {
  int count = 0;
  while (sametag[i] >= 0) {
    i = sametag[i];
    ++count;
  }
  printf("Atom ID: %ld is present %d times\n", atom->tag[i], count);
}
```

## Atom class versus AtomVec classes

The [Atom]{.title-ref} class contains all kinds of flags and counters
about atoms in the system and that includes pointers to **all** per-atom
properties available for atoms. However, only a subset of these pointers
are non-NULL and which those are depends on the atom style. For each
atom style there is a corresponding [AtomVecXXX]{.title-ref} class
derived from the [AtomVec]{.title-ref} base class, where the XXX
indicates the atom style. This [AtomVecXXX]{.title-ref} class will
update the counters and per-atom pointers if atoms are added or removed
to the system or migrate between subdomains.
