#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
g_remove_water.py

Takes a gro file containing a bilayer and remove any water molecules inside the
bilayer.

To remove molecules, it's based on Z mean of the upper and lower leaflet.
At the end, the residues and the atoms are re-numbered and the number of atoms
at the 2nd line is updated.

Title and box vectors are kept.

The atom to determine the upper and lower leaflet (typically the phosphate) can
be set with the option --lipid_atom. Same thing for the reference atom of water
(--water_atom).

For all-atoms files, you have to use the oxygen atom as a reference.

The atom to determine the upper and lower leaflet (typically the phosphate)
can be set with the option --lipid_atom.

Water molecules can also be removed if they are in a sphere centered on the
geometrical center of the atom with a given set of residue name. Use --sphere
to set the list of residue names to use for the center calculation and --radius
to set the radius of the sphere.
Water molecules removed are found based on their residue name (--water_residue).
"""

# This import makes that all divisions are real division even if the two
# members are integers. Use the euclidean division operator // to get euclidean
# divisions.
from __future__ import division

__author__ = "Marc Gueroult, Hubert Santuz & Jonathan Barnoud"
__version__ = "1.0"

import sys
import os
import argparse

# Index of the coordinate slices in gro file atom line
# The gro format is described at http://manual.gromacs.org/online/gro.html
# Columns are organized as :
# * residue number (5 positions, integer)
# * residue name (5 characters)
# * atom name (5 characters)
# * atom number (5 positions, integer)
# * position (in nm, x y z in 3 columns, each 8 positions with 3 decimal
#   places)
# * velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with
#   4 decimal places)
COORDINATE_INDEX = ((20, 28), (28, 36), (36, 44))


def isfile(path):
    """Check if path is an existing file.
    If not, raise an error. Else, return the path."""
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def define_options(argv):
    """Define the script options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", action="store", type=isfile, dest="filin",
                        required=True, metavar="coord.gro",
                        help=("The Structure file newly solvated (.gro) "
                              "REQUIRED"))
    parser.add_argument("-o", action="store", type=str, dest="filout",
                        metavar="out.gro", default="out.gro",
                        help="The Output file (.gro)")
    parser.add_argument("--lipid_atom", action="store", type=str,
                        dest="lipid_atom", metavar="P1", default="P1",
                        help=("The reference atom for the bilayer "
                              "(P1 by default)"))
    parser.add_argument("--water_atom", action="store", type=str,
                        dest="water_atom", metavar="OW", default="OW",
                        help=("The reference atom for the water. "
                              "Use the oxygen. (OW by default)"))
    parser.add_argument("--water_residue", action="store", type=str,
                        dest="water_residue", metavar="SOL", default="SOL",
                        help=("The reference residue name for the water. "
                              "(SOL by default)"))
    parser.add_argument("--sphere", type=str, nargs='+', default=None,
                        metavar="RESNAME",
                        help=("Remove water molecules if they are in a sphere "
                              "centered on the geometric center of atoms with "
                              "the given residue names. You need the "
                              "--radius option to be filled."))
    parser.add_argument("--radius", type=float, default=None,
                        help=("Remove water molecules if they are in a sphere "
                              "of this radius (in nm) centered on a given set of "
                              "residue names. You need the --sphere option to "
                              "be set."))
    args = parser.parse_args(argv)
    # --sphere and --radius have to be used together, let's check this
    if ((args.sphere is None and not args.radius is None) or
            (args.radius is None and not args.sphere is None)):
        raise parser.error(("You need to define both --sphere and --radius in "
                            "order to remove water molecules on a sphere."))
    return args


def z_mean_values(lines, atm_ref):
    """Return the Z mean for the upper and lower leaflet"""

    z_lower = z_upper = 0.0

    #Get the global Z mean
    z_sum = 0.0
    nb_atoms = nb_lower = nb_upper = 0
    for i in lines:
        atom_name = i[10:15].strip()
        if atom_name == atm_ref:
            z_sum += float(i[37:44].strip())  # Z coordinate
            nb_atoms += 1
    z_mean = z_sum / float(nb_atoms)

    #Z mean for the upper and lower leaflet
    for i in lines:
        atom_name = i[10:15].strip()
        if atom_name == atm_ref:
            z = float(i[37:44].strip())  # Z coordinate
            if z < z_mean:
                z_lower += z
                nb_lower += 1
            else:
                z_upper += z
                nb_upper += 1

    z_lower = z_lower / float(nb_lower)
    z_upper = z_upper / float(nb_upper)

    return z_lower, z_upper


def remove_water(lines, z_upper, z_lower, atom_water):
    """Return  a list containing all lines of the gro file to keep"""
    reslines = []  # Lines which will be write in the output file
    #Appends the first 2 lines (i.e. Title and Number of atoms)
    reslines.append(lines[0].strip())
    reslines.append(lines[1].strip())

    resnumber = -1
    nb_atoms = 0  # Total number of atoms
    nb_res_water = 0  # Number of waters molecules removed

    for i in lines[2:-1]:  # Discard for 2 lines and the last one
        to_store = True    # Does we store the line in the result ?
        atom_name = i[10:15].strip()
        ref_resnumber = int(i[0:5])  # Residue number of the current atoms
        if atom_name == atom_water:  # Oxygen of the water molecule
            z = float(i[37:44])
            if z > z_lower and z < z_upper:    # Water in bilayer
                # Store the residue number of the water to removed
                resnumber = int(i[0:5])
                to_store = False  # Don't store the line
                nb_res_water += 1

        # If the ref resnumer and the resnumer are equal means we don't store
        # the hydrogens of the previous water molecule
        if to_store and ref_resnumber != resnumber:
            nb_atoms += 1  # 1 line = 1 atom
            reslines.append(i[:-1])

    #Update the final number of atoms
    reslines[1] = str(nb_atoms)
    #Append the last line of the gro file (Box vectors)
    reslines.append(lines[-1])

    return reslines, nb_res_water


def geometric_center(lines, resids):
    """
    Get the isobarycenter of the atoms having the given residue id.

    The "resids" argument is a list of tuple containing the residue id and
    the numbero of atoms. 
    An atom is taken into account if its residue name is in the list.
    """

    #Get only the resids
    list_resids = [i for i,j in resids]
    #Init the dic
    dic_coords = {}
    for resid in list_resids:
        dic_coords[resid] = [0] * 3

    natoms = 0
    prev_resid = int(lines[3][0:5])

    # Sum up the coordinates of the atoms of interest
    for line in lines[2:-1]:
        resid = int(line[0:5])

        if resid in list_resids:
            for index, (begin, end) in enumerate(COORDINATE_INDEX):
                dic_coords[resid][index] += float(line[begin:end])
            natoms += 1

        if resid != prev_resid:
            prev_resid += 1

    # Do the averaging
    for resid,nb_atoms in resids:
        for i in xrange(3):
            dic_coords[resid][i] /= nb_atoms

    return dic_coords


def sq_distance(point_a, point_b):
    """
    Compute the square of the euclidean distance between two points.
    """
    return sum((a - b) * (a - b) for a, b in zip(point_a, point_b))


def remove_sphere(lines, resnames, center, radius):
    """
    Remove residue when at least one of its atoms is inside the given sphere.
    Return a list of all the gro file lines to keep.
    """

    prev_resid = -1
    nb_res_water = 0
    # We want to include the two first lines anyway
    reslines = []
    reslines.append(lines[0][:-1])
    reslines.append(lines[1][:-1])
    # Comparing squared distances instead of distances avoid to compute a bunch
    # load of square roots
    sq_radius = radius * radius
    # When an atom of interest is inside the sphere we want to remove the whole
    # residue. Therefor we must remove the previous atoms of the residue that
    # we already included and inhibit the inclusion of the following ones.
    inhibit_resid = None
    for line in lines[2:-1]:
        resname = line[5:10].strip()
        resid = int(line[0:5])
        coords = [float(line[begin:end]) for begin, end in COORDINATE_INDEX]
        if (resname in resnames
                and ((not inhibit_resid is None and resid == inhibit_resid)
                     or sq_distance(center, coords) <= sq_radius)):
            #Update the number of water molecules removed
            if (prev_resid != resid): #To avoid adding H
                nb_res_water += 1
            # The atom is inside the sphere
            # Remove previously added atoms from the residue
            while int(reslines[-1][0:5]) == resid:
                del reslines[-1]
            # Set the inhibition
            inhibit_resid = resid

            #Update the prev_resid
            prev_resid = resid
        else:
            reslines.append(line)
            inhibit_resid = None

    # We want to keep the box definition
    #reslines.append(lines[-1][:-1])
    reslines.append(lines[-1])
    # Update the number of atoms
    reslines[1] = str(len(reslines) - 3)

    return reslines, nb_res_water


def find_ref_atom(lines, ref_atom):
    """Return a boolean whether the "atom" was find in the file"""

    for i in lines:
        atom_name = i[10:15].strip()
        if ref_atom == atom_name:
            return True
    return False

def find_ref_residue(lines, ref_residue):
    """Return a boolean whether the "residue" was find in the file"""

    for i in lines:
        residue_name = i[5:10].strip()
        if ref_residue == residue_name:
            return True
    return False

def get_resids(lines, resnames):
    """Return list of tuples containing the resid 
    and the number of atomes from the resnames"""

    res_ids = []
    prev_resid = int(lines[3][0:5])
    nb_atoms = 0

    for line in lines[2:-1]:
        resname = line[5:10].strip()
        resid = int(line[0:5])

        if resname in resnames:
            if resid != prev_resid and nb_atoms != 0:
                res_ids.append((prev_resid, nb_atoms))
                nb_atoms = 1
            else:
                nb_atoms += 1

        if resid != prev_resid:
            prev_resid += 1

    #If the molecule 'resnames' is the last one
    if lines[-2][5:10].strip() in resnames:
        res_ids.append((resid, nb_atoms))

    return res_ids

def renumber(lines, start_res=None):
    """Renum the atoms and the residues in a file"""

    out = []
    #Appends the first 2 lines (ie Title and Number of atoms)
    out.append(lines[0].strip())
    out.append(lines[1].strip())

    start = True
    output_atom = 0

    for i in lines[2:-1]:  # Discard for 2 lines and the last one
        # If first residue, start the counter regarding the number of the
        # starting residue
        if start:
            start = False
            current_res = i[0:5]
            if start_res:
                output_res = start_res
            else:
                output_res = 1

        #New residue
        if i[0:5] != current_res:
            if output_res == 99999:  # cannot go to 100000
                output_res = 1
            else:
                output_res += 1
            current_res = i[0:5]

        if output_atom == 99999:  # cannot go to 100000
            output_atom = 1
        else:
            output_atom += 1

        #Append the line with the correct number of residues and atoms
        out.append("%5i%s%5i%s" % (output_res, i[5:15], output_atom, i[20:]))

    #Save the last line (box vectors)
    out.append(lines[-1])

    return out


if __name__ == '__main__':

    #Command line parsing
    args = define_options(sys.argv[1:])

    filin = args.filin
    filout = args.filout
    lipid_atom = args.lipid_atom
    water_atom = args.water_atom


    print "The input coordinate file is {0}".format(filin)
    print "The output file will be {0}".format(filout)
    print "The reference atom for the lipid bilayer is {0}".format(lipid_atom)
    print "The reference atom for the water is {0}".format(water_atom)

    f = open(filin, 'r')
    data = f.readlines()
    f.close()

    wat = 0
    if args.sphere is None:
        print
        print "Checking the reference atom...",
        if not find_ref_atom(data, lipid_atom):
            print "Oops!"
            print ("The reference atom {0} for the bilayer was not find. "
                   "Exiting...").format(lipid_atom)
            sys.exit()
        if not find_ref_atom(data, water_atom):
            print "Oops!"
            print ("The reference atom {0} for the water was not find. "
                   "Exiting...").format(water_atom)
            sys.exit()
            print "Done!"

        #Get Z mean for  the upper and lower leaflet
        z_lower, z_upper = z_mean_values(data, lipid_atom)

        print "Removing water inside the bilayer...",
        #Remove water molecules inside the bilayer
        temp_lines, wat = remove_water(data, z_upper, z_lower,
                                                 water_atom)
        print "Done!"

        first_res_number = int(data[2][0:5])
        print "Renumber residues and atoms...",
        output = renumber(temp_lines, first_res_number)
        print "Done!"
    else:
        print "The reference residue for the sphere is {0}".format(", ".join(args.sphere))
        print "The reference residue for the water is {0}".format(args.water_residue)
        if not find_ref_residue(data, args.water_residue):
            print "Oops!"
            print ("The reference residue {0} for water  was not find. "
                   "Exiting...").format(args.water_residue)
            sys.exit()
        for res in args.sphere:
            if not find_ref_residue(data, res):
                print "Oops!"
                print ("The reference residue {0} for the sphere  was not find. "
                       "Exiting...").format(res)
                sys.exit()

        print "Get the center of the sphere...",
        resids = get_resids(data, args.sphere)
        print resids
        dic_center = geometric_center(data, resids)
        print "Done! The center is {0}".format(dic_center)
        print "Remove water molecules inside the sphere...",
        output=data
        for resid,center in dic_center.items():
            temp_lines, water_removed = remove_sphere(output, args.water_residue, center, args.radius)
            print "Done!"
            wat += water_removed
            first_res_number = int(data[2][0:5])
            print "Renumber residues and atoms...",
            output = renumber(temp_lines, first_res_number)
            print "Done!"

    print
    print "The old system contained {0} atoms.".format(data[1].strip())
    print "{0} water molecules have been removed.".format(wat)
    print "The new system contains {0} atoms.".format(output[1])

    #Write in the output file
    f = open(filout, 'w')
    f.write('\n'.join(output))
    f.write("\n")
    f.close()
