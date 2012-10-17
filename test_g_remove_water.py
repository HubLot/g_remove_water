#!/usr/bin/env python
#-*- coding:utf-8 -*-

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
Unit tests for g_remove_water.
"""

# This import makes that all divisions are real division even if the two
# members are integers. Use the euclidean division operator // to get euclidean
# divisions.
from __future__ import division
from __future__ import with_statement
from unittest import TestCase, main
import os
import g_remove_water as grw

__author__ = "Marc Gueroult, Hubert Santuz & Jonathan Barnoud"
__version__ = "1.0"

# Path to the directory containing the material for tests
REFDIR = "test_ressources"

class TestGeneral(TestCase) :
    """
    Test genral infrastructure that is shared with every modes.
    """
    def test_find_ref_atom_regular(self) :
        """
        Test find_ref_atom function in the regular case where the reference
        atom is present or absent from the file body.
        """
        path = os.path.join(REFDIR, "regular.gro")
        lines = open(path).readlines()
        # P1 is present in the file
        self.assertTrue(grw.find_ref_atom(lines, "P1"),
                "P1 should be find in {0} body.".format(path))
        # ABS is absent from the file
        self.assertFalse(grw.find_ref_atom(lines, "ABS"),
                "ABS should not be find in {0} body.".format(path))
        # AA is present at the right place in the header but should be ignored
        self.assertFalse(grw.find_ref_atom(lines, "AA"),
                "AA should be ignored in {0} head.".format(path))

    def test_renumber_noarg(self) :
        """
        Test the atom renumbering with the renumber function without the
        start_res argument.
        """
        path = os.path.join(REFDIR, "regular.gro")
        removed_res = (10, 50, 60)
        # Remove some residues and renumber atoms and residues
        with open(path) as infile :
            renumbered = _create_runumbered(infile, removed_res)
        # Check numbering
        # (number of atom per residue, number of residue)
        topology = ((52, 72 - len(removed_res)), (3, 3739))
        _test_renumber(renumbered, topology, 1)

    def test_renumber_start_res(self) :
        """
        Test the atom renumbering with the renumber function with the
        start_res argument. This is testing both the start-res argument and the
        ability to handle residue number bigger than 99999.
        """
        path = os.path.join(REFDIR, "regular.gro")
        removed_res = (10, 50, 60)
        start_res = 98999
        # Remove some residues and renumber atoms and residues
        with open(path) as infile :
            renumbered = _create_runumbered(infile, removed_res, start_res)
        # Check numbering
        # (number of atom per residue, number of residue)
        topology = ((52, 72 - len(removed_res)), (3, 3739))
        _test_renumber(renumbered, topology, start_res)


class TestSlice(TestCase) :
    """
    Test slice mode specific functions.
    """
    pass


class TestSphere(TestCase) :
    """
    Test sphere mode specific functions.
    """
    pass


class TestProgramm(TestCase) :
    """
    Test the program command line.
    """
    pass


def residue_numbers(topology, start_res=1) :
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
    for natoms, nresidues in topology :
        for residue in xrange(nresidues) :
            resid += 1
            if resid > 99999 :
                resid = 1
            for atoms in xrange(natoms) :
                yield resid

def _create_runumbered(infile, removed_res, start_res=None) :
    """
    Remove residues from a structure and renumber the atoms and residues.

    :Parameters:
        - infile: the opened file to read
        - remove_res: a list of resid to remove from the structure
        - start_res: the initial residue number in the renumbered structure

    :Returns:
        - the lines from the renumbered structure
    """
    lines = infile.readlines()
    # Remove some residues
    keep = lines[:2] + [line for line in lines[2:-1]
            if not int(line[0:5]) in removed_res] + [line[-1]]
    # Renumber residues and atoms
    renumbered = grw.renumber(keep, start_res)
    return renumbered

def _test_renumber(lines, topology, start_res) :
    """
    Test atom renumbering in various conditions.

    :Parameters:
        - start_res: the initial residue number in the renumbered structure
    """
    for line_number, (ref_resid, line) \
            in enumerate(zip(residue_numbers(topology, start_res),
                lines[2:-1])) :
        resid = int(line[0:5])
        atomid = int(line[15:20])
        ref_atomid = line_number + 1
        print line[:-1]
        # Check the residue
        assert resid == ref_resid,\
                ("Resisue ID is wrong after renumbering: "
                "{0} instead of {1} at line {2}").format(
                    resid, ref_resid, line_number + 3)
        # Check the atom
        assert atomid == ref_atomid,\
                ("Atom ID is wrong after renumbering: "
                "{0} instead of {1} at line {2}").format(
                    atomid, ref_atomid, line_number + 3)

if __name__ == "__main__" :
    main()
