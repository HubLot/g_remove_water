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
Unit tests for g_remove_water.
"""

# This import makes that all divisions are real division even if the two
# members are integers. Use the euclidean division operator // to get euclidean
# divisions.
from __future__ import division
from __future__ import with_statement
import unittest
import os
import sys
import contextlib
import groIO
import g_remove_water as grw

# Python 2.X/3.X compatibility
try:
    range = xrange
except:
    pass

__author__ = "Marc Gueroult, Hubert Santuz & Jonathan Barnoud"
__version__ = "1.0"

# Path to the directory containing the material for tests
REFDIR = "test_resources"

# Gro files usually store coordinates with decimals. Let's be precise
# until the fourth one.
PRECISION = 4


class TestGeneral(unittest.TestCase):
    """
    Test general infrastructure that is shared with every modes.
    """
    def test_missing_file(self):
        """
        Test if the input file is missing
        """
        with _redirect_stderr(sys.stdout):
            argv = ["-f", "{0}/pouet.gro".format(REFDIR)]
            self.assertRaises(SystemExit, grw.define_options, *[argv])

    def test_find_ref_atom_regular(self):
        """
        Test find_ref_atom function in the regular case where the reference
        atom is present or absent from the file body.
        """
        path = os.path.join(REFDIR, "regular.gro")
        title, atoms, box = groIO.parse_file(path)

        # P1 is present in the file
        self.assertTrue(grw.find_ref_atom(atoms, "P1"),
                        "P1 should be find in {0} body.".format(path))
        # ABS is absent from the file
        self.assertFalse(grw.find_ref_atom(atoms, "ABS"),
                         "ABS should not be find in {0} body.".format(path))
        # AA is present at the right place in the header but should be ignored
        self.assertFalse(grw.find_ref_atom(atoms, "AA"),
                         "AA should be ignored in {0} head.".format(path))

    def test_raise_ref_objects(self):
        """
        Test if the exception is raised when reference objects (atom or residue)
        are missing in the file.
        """

        path = os.path.join(REFDIR, "regular.gro")
        title, atoms, box = groIO.parse_file(path)

        # Lipid atom
        with self.assertRaises(grw.RefObjectError) as context:
            grw.perform_bilayer_removing(atoms, "P12", "OW", verbose=False)

        # Water atom
        with self.assertRaises(grw.RefObjectError) as context:
            grw.perform_bilayer_removing(atoms, "P1", "IEW", verbose=False)

        # Sphere residue
        with self.assertRaises(grw.RefObjectError) as context:
            grw.perform_sphere_removing(atoms, "AZE", 0.4, "SOL", verbose=False)

        # Water residue
        with self.assertRaises(grw.RefObjectError) as context:
            grw.perform_sphere_removing(atoms, "POPC", 0.4, "SO", verbose=False)


class TestSlice(unittest.TestCase):
    """
    Test slice mode specific functions.
    """
    def test_z_mean_values(self):
        """
        Test the mean Z values returned for each leaflet by z_mean_values.

        The result of the function is compared with the result obtained using
        the g_traj gromacs tool with the following command::

            echo 0 1 | g_traj -f regular.gro -s regular.gro -n leaflets.ndx  \
                    -ox regular_com -com -ng 2
        """
        reference = (2.14186, 5.71483)
        # For the record, values for X are (2.23333, 2.5495), values for Y are
        # (2.688, 2.50058)

        path = os.path.join(REFDIR, "regular.gro")
        title, atoms, box = groIO.parse_file(path)

        z_low, z_top = grw.z_mean_values(atoms, "P1")
        self.assertAlmostEqual(reference[0], z_low, PRECISION,
                               ("Z mean value for the lower leaflet do not "
                                "match with {0} places: {1} instead of {2}")
                               .format(PRECISION, reference[0], z_low))
        self.assertAlmostEqual(reference[1], z_top, PRECISION,
                               ("Z mean value for the upper leaflet do not "
                                "match with {0} places: {1} instead of {2}")
                               .format(PRECISION, reference[1], z_top))


class TestSphere(unittest.TestCase):
    """
    Test sphere mode specific functions.
    """
    def test_get_resids(self):
        """
        Test the mapping of resids from resnames
        """
        path = os.path.join(REFDIR, "regular.gro")
        title, atoms, box = groIO.parse_file(path)

        resname = ["POPC"]
        resids = list(range(1,73))

        self.assertEqual(grw.get_resids(atoms, resname), resids)

    def test_sq_distance(self):
        """
        Test square distance calculation without periodic boundary conditions.
        """
        # format is (point A, point B, expected square distance)
        reference = (
                    ((0, 0, 0), (0, 0, 0), 0),
                    ((1, 1, 1), (1, 1, 1), 0),
                    ((0, 0, 0), (1, 0, 0), 1),  # Check X axis
                    ((0, 0, 0), (0, 1, 0), 1),  # Check Y axis
                    ((0, 0, 0), (0, 0, 1), 1),  # Check Z axis
                    ((0.978, 4.892, 3.889), (2.261, 4.973, 2.441), 3.749354),
                    ((-0.044, 3.691, 2.737), (0.227, 4.160, 2.714), 0.2939),
        )
        for point_a, point_b, ref_value in reference:
            value = grw.sq_distance(point_a, point_b)
            self.assertAlmostEqual(value, ref_value, PRECISION,
                                   ("Square distance between {0} and {1} do "
                                    "not match expectations with {2} places: "
                                    "{3} instead of {4}.").format(
                                        point_a, point_b, PRECISION,
                                        value, ref_value))

    def test_geometric_center(self):
        """
        Test the geometric center calculation.

        Reference coordinates were calculated with g_traj. The output of the
        program is available in test_resources/center.xvg.
        """
        reference = {"x": 0.000216667, "y": 0.00045, "z": 5.00003e-05}
        resids = [1]
        path = os.path.join(REFDIR, "center.gro")
        title, atoms, box = groIO.parse_file(path)

        center = grw.geometric_center(atoms, resids)

        for ref, value in zip(reference.values(), center[resids[0]].values()):
            self.assertAlmostEqual(ref, value, PRECISION,
                                   ("Geometric center is wrong: "
                                    "{0} instead of {1}")
                                   .format(center, reference))

    def test_define_options_sphere_radius_type(self):
        """
        Test if --sphere and --radius options return the right types.
        """
        sphere_ref = ["ARG", "LEU", "GLY"]
        argv = ["-f", "{0}/regular.gro".format(REFDIR), "--sphere"] + \
            sphere_ref + ["--radius", "6.4"]
        args = grw.define_options(argv)
        self.assertEqual(args.sphere, sphere_ref,
                         "Option --sphere get {0} instead of {1}.".format(
                             args.sphere, sphere_ref))
        self.assertAlmostEqual(args.radius, 6.4, PRECISION,
                               "Option --radius get {0} instead of {1}."
                                   .format(args.radius, 6.4))

    def test_define_options_radius_non_nul(self):
        """
        Test if define_options complains about nul --radius.
        """
        with _redirect_stderr(sys.stdout):
            argv = ["-f", "{0}/regular.gro".format(REFDIR),
                    "--sphere", "ARG", "--radius", "0.0"]
            self.assertRaises(SystemExit, grw.define_options, *[argv])

    def test_define_options_sphere_radius_missing(self):
        """
        Test if define_options complains about --sphere or --radius missing.
        """
        # Argparse writes in stderr and that screw nosetests output capture.
        # Let's redirect stderr to stdout so argparse output can be
        # captured by nosetests.
        with _redirect_stderr(sys.stdout):
            argv = ["-f", "{0}/regular.gro".format(REFDIR), "--sphere", "ARG"]
            self.assertRaises(SystemExit, grw.define_options, *[argv])
            argv = ["-f", "{0}/regular.gro".format(REFDIR), "--radius", "4.2"]
            self.assertRaises(SystemExit, grw.define_options, *[argv])


class TestProgramm(unittest.TestCase):
    """
    Test the program command line.
    """
    pass


@contextlib.contextmanager
def _redirect_stderr(destination=sys.stdout):
    """
    Redirect sys.stderr to an other file descriptor (sys.stdout by default).

    This function is a context manager and is used like as followed::

        >>> with _redirect_stderr():
        ...     print >> sys.stderr, "Something"

    In this example the print is done in sys.stdout even if sys.stderr is
    explicitly mentioned.

    Stderr can also be redirected to a file descriptor::

        >>> with open("my_file", "wt") as outfile:
        ...     with _redirect_stderr(outfile):
        ...         print >> sys.stderr, "Something"

    """
    old_stderr = sys.stderr
    sys.stderr = destination
    yield
    sys.stderr = old_stderr


if __name__ == "__main__":
    unittest.main()
