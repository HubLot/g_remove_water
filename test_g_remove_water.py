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

if __name__ == "__main__" :
    main()
