# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-
# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
#
# MDAnalysis --- https://www.mdanalysis.org
# Copyright (c) 2006-2017 The MDAnalysis Development Team and contributors
# (see the file AUTHORS for the full list of names)
#
# Released under the Lesser GNU Public Licence, v2.1 or any higher version
#
# Please cite your use of MDAnalysis in published work:
#
# R. J. Gowers, M. Linke, J. Barnoud, T. J. E. Reddy, M. N. Melo, S. L. Seyler,
# D. L. Dotson, J. Domanski, S. Buchoux, I. M. Kenney, and O. Beckstein.
# MDAnalysis: A Python package for the rapid analysis of molecular dynamics
# simulations. In S. Benthall and S. Rostrup editors, Proceedings of the 15th
# Python in Science Conference, pages 102-109, Austin, TX, 2016. SciPy.
# doi: 10.25080/majora-629e541a-00e
#
# N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein.
# MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
# J. Comput. Chem. 32 (2011), 2319--2327, doi:10.1002/jcc.21787
#

r"""
ITP topology parser
===================

Writes a GROMACS ITP_ file. The file will
contain atom IDs, segids, residue IDs, residue names, atom names, atom types,
charges, chargegroups, masses, moltypes, and molnums.
Any missing values will be guessed.
Bonds, angles, dihedrals and impropers are also read from the file.

If an ITP file is passed without a ``[ molecules ]`` directive, passing 
``infer_system=True`` (the default option) will create a Universe with 
1 molecule of each defined ``moleculetype``. 
If a ``[ molecules ]`` section is present, ``infer_system`` is ignored.

If files are included with the `#include` directive, they will also be read.
If they are not in the working directory, ITPParser will look for them in the
``include_dir`` directory. By default, this is set to
``include_dir="/usr/local/gromacs/share/gromacs/top/"``.
Variables can be defined with the `#define` directive in files, or by passing
in :ref:`keyword arguments <itp-define-kwargs>`.

Examples
--------

.. code-block:: python

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import ITP_tip5p

    #  override charge of HW2 atom (defined in file as HW2_CHARGE)
    u = mda.Universe(ITP_tip5p, HW2_CHARGE=2, infer_system=True)

.. note::

    AMBER also uses topology files with the .top extension. To use ITPParser
    to read GROMACS top files, pass ``topology_format='ITP'``.

    ::

        import MDAnalysis as mda
        u = mda.Universe('topol.top', topology_format='ITP')


.. _itp-define-kwargs:

Preprocessor variables
----------------------

ITP files are often defined with lines that depend on 
whether a keyword flag is given. For example, this modified TIP5P water file:

.. code-block:: none

    [ moleculetype ]
    ; molname       nrexcl
    SOL             2

    #ifndef HW1_CHARGE
        #define HW1_CHARGE 0.241
    #endif

    #define HW2_CHARGE 0.241

    [ atoms ]
    ; id    at type res nr  residu name     at name         cg nr   charge
    1       opls_118     1       SOL              OW             1       0
    2       opls_119     1       SOL             HW1             1       HW1_CHARGE
    3       opls_119     1       SOL             HW2             1       HW2_CHARGE
    4       opls_120     1       SOL             LP1             1      -0.241
    5       opls_120     1       SOL             LP2             1      -0.241
    #ifdef EXTRA_ATOMS  ; added for keyword tests
    6       opls_120     1       SOL             LP3             1      -0.241
    7       opls_120     1       SOL             LP4             1       0.241
    #endif


Define these preprocessor variables by passing keyword arguments. Any arguments that you 
pass in *override* any variables defined in the file. For example, the universe below 
will have charges of 3 for the HW1 and HW2 atoms::

    import MDAnalysis as mda
    from MDAnalysis.tests.datafiles import ITP_tip5p

    u = mda.Universe(ITP_tip5p, EXTRA_ATOMS=True, HW1_CHARGE=3, HW2_CHARGE=3)

These keyword variables are **case-sensitive**. Note that if you set keywords to
``False`` or ``None``, they will be treated as if they are not defined in #ifdef conditions.

For example, the universe below will only have 5 atoms. ::

    u = mda.Universe(ITP_tip5p, EXTRA_ATOMS=False)


.. _ITP: http://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#molecule-itp-file

.. _TOP: http://manual.gromacs.org/current/reference-manual/file-formats.html#top

Classes
-------

.. autoclass:: ITPParser
   :members:
   :inherited-members:

"""


from .base import TopologyWriterBase

class ITPWriter(TopologyWriterBase):
    """Writes topology information in a GROMACS ITP_ file.

    Creates a Topology with the following Attributes:
    - ids
    - names
    - types
    - masses
    - charges
    - chargegroups
    - resids
    - resnames
    - segids
    - moltypes
    - molnums
    - bonds
    - angles
    - dihedrals
    - impropers

    .. _ITP: http://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html#molecule-itp-file
    """

    format = "ITP"
    units = {"time": None, "length": "nm", "velocity": "nm/ps"}

    #: format strings for the GRO file (all include newline); precision
    #: of 3 decimal places is hard-coded here.
    fmt = {
        "n_atoms": "{0:5d}\n",  # number of atoms
    }


    def __init__(self, filename, **kwargs):
        """Initializes ITPWriter object.

        Parameters
        ----------
        filename : str
            Name of the output ITP file.
        """
        self.filename = filename
        self.kwargs = kwargs

        self.include_directives = self.kwargs.pop("include_directives", None)


    def write_include_directives(self):
        with open(self.filename, 'a') as f:
            f.write("\n; Include directives\n")
            if len(self.include_directives) == 1:
                f.write(f"#include \"{self.include_directives}\"\n")
            else:
                for include in self.include_directives:
                    f.write(f"#include \"{include}\"\n")


    def write_moleculetypes_directive(self,):
        topology = self.topology
        mol_name = topology.moltypes[0] if topology.moltypes else "MOL"
        with open(self.filename, 'w') as f:
            f.write("[ moleculetypes ]")
            f.write(f"""
; molname       nrexcl
{mol_name}            2
            """)

            f.write("[ atoms ]")
            f.write("""
atom type; residue number; residue name; atom name; charge group number; q(e); m(u)
            """)


    def write_atoms_directive(self,):
        atoms = self.topology.atoms
        with open(self.filename, 'a') as f:
            for atom in atoms:
                f.write(f"{atom.id:5d} {atom.type:<10} {atom.resid:5d} "
                        f"{atom.resname:<5} {atom.name:<5} {atom.chargegroup:5d} "
                        f"{atom.charge:8.4f} {atom.mass:8.4f}\n")

    def write(self, universe):
        """Writes the topology to an ITP file.

        Parameters
        ----------
        topology : Topology
            The Topology object to write to the ITP file.
        """
        self.topology = topology

        if self.include_directives is not None:                
            self.write_include_directives()
        self.write_molecules()

        print("Success!")
        # if "include_directives" in self.kwargs:
        #     self.write_include_directives()
        
