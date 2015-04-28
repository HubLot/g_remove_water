#!/usr/bin/env python
# -*- coding:utf-8 -*-

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""
Unit tests for groIO library
"""

from __future__ import division, print_function
import unittest
import os
import sys
import contextlib
import groIO

__author__ = "Hubert Santuz & Jonathan Barnoud"


# Path to the directory containing the material for tests
REFDIR = "test_resources"

# Gro files usually store coordinates with decimals. Let's be precise
# until the fourth one.
PRECISION = 4


class TestReadWrite(unittest.TestCase):
    """
    Test for reading and writing a GRO file.
    """

    def test_read(self):
        """
        Test the reading of a regular gro file
        """

        path = os.path.join(REFDIR, "regular.gro")
        title, atoms, box = groIO.parse_file(path)

        nb_atoms = len(atoms)

        self.assertEqual(title, "Regular      AA\n")
        self.assertEqual(box, "   5.15414   5.15414   7.93645\n")
        self.assertEqual(nb_atoms, 14961)

        # Pick a random atom
        # 24POPC  C220 1235   2.520   4.888   3.113
        atom = atoms[1234]
        keys_tested = ['resid', 'atomid', 'x', 'z']
        values_tested = [24, 1235, 2.520, 3.113]

        for key, value in zip(keys_tested, values_tested):
            self.assertEqual(atom[key], value)

    def test_fail_read(self):
        """
        Test the bad formating of a gro file
        """

        files = ["fail.gro", "fail_100000.gro", "missing_total_atoms.gro"]
        files_desc = [os.path.join(REFDIR, filin) for filin in files]
        print(files_desc)

        for filin in files_desc:
            with self.assertRaises(groIO.FormatError) as context:
                groIO.parse_file(filin)
                self.assertEqual("Something is wrong in the format", context.exception)

    def test_write(self):
        """
        Test the writing of a gro file
        """

        # Create a random file
        title, atoms, box = _generate_file()

        # Write it
        test_write = os.path.join(REFDIR, "test_write.gro")
        with open(test_write, "w") as fout:
            for line in groIO.write_gro(title, atoms, box):
                print(line, end='', file=fout)

        # Reference file
        ref_write = os.path.join(REFDIR, "write.gro")

        with open(ref_write) as f1, open(test_write) as f2:
            ref_readlines = f1.readlines()
            test_readlines = f2.readlines()

            self.assertEqual(ref_readlines, test_readlines)

        os.remove(test_write)


def _generate_file():
    """
    Generate the header and atoms of a random gro file.
    Match the file write.gro to test purposes.

    :Returns:
        - title: the title of the system
        - atoms: a list of atom, each atom is stored as a dictionary
        - box: the box description
    """

    title = "Write\n"
    atoms = []
    atoms.append({'resid': 1, 'resname': "POPC", 'atom_name': "C31",
                  'atomid': 1, 'x': 1.764,  'y': 4.587, 'z': 2.046})
    atoms.append({'resid': 1, 'resname': "POPC", 'atom_name': "N1",
                  'atomid': 2, 'x': 1.824,  'y': 4.555, 'z': 1.916})
    atoms.append({'resid': 1, 'resname': "POPC", 'atom_name': "C32",
                  'atomid': 3, 'x': 1.755,  'y': 4.436, 'z': 1.864})
    atoms.append({'resid': 1, 'resname': "POPC", 'atom_name': "C33",
                  'atomid': 4, 'x': 1.954,  'y': 4.503, 'z': 1.960})

    box = "   1.000   1.000   1.000\n"

    return (title, atoms, box)


class TestGrolib(unittest.TestCase):
    """
    Tests for the other functions in the library.
    """

    def test_renumber(self):
        """
        Test the atom renumbering with the renumber function
        """
        path = os.path.join(REFDIR, "regular.gro")
        title, atoms, box = groIO.parse_file(path)

        removed_res = (10, 50, 60)
        # Remove some residues and renumber atoms and residues
        renumbered = _create_runumbered(atoms, removed_res)
        # Check numbering
        # (number of atom per residue, number of residue)
        topology = ((52, 72 - len(removed_res)), (3, 3739))
        _test_renumber(renumbered, topology)


def residue_numbers(topology, start_res=1):
    """
    Generate residue numbers according to a topology.

    Topology is a list of successive succession of residues described as
    (number of atom per residue, number of residue). For instance, a succession
    of 8 residue of 10 atoms each followed by 5 residues of 3 atoms each is
    described as ((10, 8), (3, 5)).

    :Parameters:
        - topology: the residue succession as described above
        - start_res: the number of the first residue
    """
    resid = start_res - 1
    for natoms, nresidues in topology:
        for residue in range(nresidues):
            resid += 1
            if resid > 99999:
                resid = 1
            for atoms in range(natoms):
                yield resid


def _create_runumbered(atoms, removed_res):
    """
    Remove residues from a structure and renumber the atoms and residues.

    :Parameters:
        - atoms: the list of dictionnary for atoms
        - remove_res: a list of resid to remove from the structure

    :Returns:
        - the new list renumbered
    """
    # Remove some residues
    keep = [atom for atom in atoms if not atom['resid'] in removed_res]
    # Renumber residues and atoms
    renumbered = groIO.renumber(keep)
    return renumbered


def _test_renumber(atoms, topology):
    """
    Test atom renumbering.

    :Parameters:
        - atoms: the list of dictionnary for atoms
        - topology: the residue succession, see :func:`residue_numbers`
    """
    for line_number, (ref_resid, atom) \
            in enumerate(zip(residue_numbers(topology),
                             atoms)):
        resid = atom["resid"]
        atomid = atom["atomid"]
        ref_atomid = line_number + 1
        # Check the residue
        assert resid == ref_resid, \
            ("Residue ID is wrong after renumbering: "
             "{0} instead of {1} at line {2}").format(
                 resid, ref_resid, line_number + 3)
        # Check the atom
        assert atomid == ref_atomid, \
            ("Atom ID is wrong after renumbering: "
             "{0} instead of {1} at line {2}").format(
                 atomid, ref_atomid, line_number + 3)


if __name__ == "__main__":
    unittest.main()
