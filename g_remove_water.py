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


from __future__ import division, print_function

__author__ = "Marc Gueroult, Hubert Santuz & Jonathan Barnoud"
__version__ = "1.0"

import sys
import os
import argparse

# Handle GRO format
import groIO

# Exception to handle the missing references objects (atom or residue)
class RefObjectError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return self.value


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
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Increase verbosity")
    args = parser.parse_args(argv)
    # --sphere and --radius have to be used together, let's check this
    if ((args.sphere is None and args.radius is not None) or
            (args.radius is None and args.sphere is not None)):
        raise parser.error(("You need to define both --sphere and --radius in "
                            "order to remove water molecules on a sphere."))
    return args


def z_mean_values(atoms, atom_ref):
    """Return the Z mean for the upper and lower leaflet"""

    z_lower = z_upper = 0.0

    # Get the global Z mean
    z_sum = 0.0
    nb_atoms = nb_lower = nb_upper = 0
    for atom in atoms:
        if atom['atom_name'] == atom_ref:
            z_sum += atom['z']
            nb_atoms += 1
    z_mean = z_sum / float(nb_atoms)

    # Z mean for the upper and lower leaflet
    for atom in atoms:
        if atom['atom_name'] == atom_ref:
            if atom['z'] < z_mean:
                z_lower += atom['z']
                nb_lower += 1
            else:
                z_upper += atom['z']
                nb_upper += 1

    z_lower = z_lower / float(nb_lower)
    z_upper = z_upper / float(nb_upper)

    return z_lower, z_upper


def remove_water(atoms, z_upper, z_lower, atom_water):
    """Return  a list containing all atoms of the gro file to keep"""

    new_atoms = []  # atoms which will be write in the output file

    atom_resid = -1
    nb_atoms = 0  # Total number of atoms
    nb_res_water = 0  # Number of waters molecules removed

    for atom in atoms: 
        to_store = True    # Does we store the line in the result ?
        ref_resid = atom['resid']
        if atom['atom_name'] == atom_water: 
            if z_lower < atom['z'] < z_upper:    # Water in bilayer
                # Store the residue number of the water to removed
                atom_resid = atom['resid']
                to_store = False
                nb_res_water += 1

        # If the ref resnumer and the resnumer are equal means we don't store
        # the hydrogens of the previous water molecule
        if to_store and ref_resid != atom_resid:
            nb_atoms += 1  # 1 line = 1 atom
            new_atoms.append(atom)


    return new_atoms, nb_res_water


def geometric_center(atoms, resids):
    """
    Get the isobarycenter of the atoms having the given residue id.

    The "resids" argument is a list containing the residue id
    An atom is taken into account if its residue name is in the list.
    """

    dic_coords = {}

    for resid in resids:
        dic_coords[resid] = {}
        for dim in ['x', 'y', 'z']:
            coords = [atom[dim] for atom in atoms if atom['resid'] == resid]
            dic_coords[resid][dim] = sum(coords)/len(coords)

    return dic_coords


def sq_distance(point_a, point_b):
    """
    Compute the square of the euclidean distance between two points.
    """
    return sum((a - b) * (a - b) for a, b in zip(point_a.values(), point_b.values()))


def remove_sphere(atoms, resnames, center, radius):
    """
    Remove residue when at least one of its atoms is inside the given sphere.
    Return a list of all the gro file atoms to keep.
    """

    prev_resid = -1
    nb_res_water = 0

    res_atoms = []
    # Comparing squared distances instead of distances avoid to compute a bunch
    # load of square roots
    sq_radius = radius * radius
    # When an atom of interest is inside the sphere we want to remove the whole
    # residue. Therefor we must remove the previous atoms of the residue that
    # we already included and inhibit the inclusion of the following ones.
    inhibit_resid = None
    for atom in atoms:
        resid = atom['resid']
        coords = {key: atom[key] for key in atom.keys() if key in ['x', 'y', 'z']}
        if (atom['resname'] in resnames
                and ((inhibit_resid is not None and resid == inhibit_resid)
                     or sq_distance(center, coords) <= sq_radius)):
            # Update the number of water molecules removed
            if (prev_resid != resid):  # To avoid adding H
                nb_res_water += 1
            # The atom is inside the sphere
            # Remove previously added atoms from the residue
            while res_atoms[-1]['resid'] == resid:
                del res_atoms[-1]
            # Set the inhibition
            inhibit_resid = resid

            # Update the prev_resid
            prev_resid = resid
        else:
            res_atoms.append(atom)
            inhibit_resid = None

    return res_atoms, nb_res_water


def find_ref_atom(atoms, ref_atom):
    """Return a boolean whether the "atom" was find in the file"""

    for atom in atoms:
        if ref_atom == atom['atom_name']:
            return True
    return False


def find_ref_residue(atoms, ref_residue):
    """Return a boolean whether the "residue" was find in the file"""

    for atom in atoms:
        if ref_residue == atom['resname']:
            return True
    return False


def get_resids(atoms, res_name):
    """Return a list containing the resid from the res_name"""

    res_ids = [atom['resid'] for atom in atoms if atom['resname'] in res_name]

    return list(set(res_ids))


def renumber(atoms):
    """Renum the atoms and the residues in a file"""

    new_atoms = []
    resid = 0
    prev_resid = 0
    for atomid, atom in enumerate(atoms, start=1):
        if atom['resid'] != prev_resid:
            resid += 1
            prev_resid = atom['resid']
        atom['resid'] = resid%100000
        atom['atomid'] = atomid%100000
        new_atoms.append(atom)

    return new_atoms


def perform_bilayer_removing(atoms, ref_lipid_atom, ref_water_atom, verbose):
    """ Remove the water molecules inside a bilayer"""

    if verbose:
        print()
        print("Checking the reference atom...")
    if not find_ref_atom(atoms, ref_lipid_atom):
        raise RefObjectError(("The reference atom {0} for the bilayer was not find. "
                              "Exiting...").format(ref_lipid_atom))

    if not find_ref_atom(atoms, ref_water_atom):
        raise RefObjectError(("The reference atom {0} for the water was not find. "
                              "Exiting...").format(ref_water_atom))

    # Get Z mean for  the upper and lower leaflet
    z_lower, z_upper = z_mean_values(atoms, ref_lipid_atom)

    if verbose:
        print("Removing water inside the bilayer...")

    # Remove water molecules inside the bilayer
    temp_atoms, nb_water_removed = remove_water(atoms, z_upper, z_lower, ref_water_atom)

    if verbose:
        print("Renumber residues and atoms...")
    output = renumber(temp_atoms)

    return (output, nb_water_removed)


def perform_sphere_removing(atoms, sphere_residus, sphere_radius, ref_water_residue, verbose):
    """
    Remove the water molecules inside a sphere centered on the
    geometrical center of the atom with a given set of residue name
    """

    if verbose:
        print()
        print("The reference residue for the sphere is {0}".format(", ".join(sphere_residus)))
        print("The reference residue for the water is {0}".format(ref_water_residue))

    for residue in sphere_residus:
        if not find_ref_residue(atoms, residue):
            raise RefObjectError(("The reference residue {0} for the sphere was not find. "
                                  "Exiting...").format(residue))

    if not find_ref_residue(atoms, ref_water_residue):
        raise RefObjectError(("The reference residue {0} for water was not find. "
                              "Exiting...").format(ref_water_residue))

    resids = get_resids(atoms, sphere_residus)
    dic_center = geometric_center(atoms, resids)

    if not verbose:
        print("Removing water molecules...")

    output_atoms = atoms
    sum_water_removed = 0
    for resid, center in dic_center.items():
        if verbose:
            print("Removing water inside the sphere centered ({x:.1f} {y:.1f} {z:.1f})"
                  " on residue {res}".format(res=resid, **center))

        tmp_atoms, nb_water_removed = remove_sphere(output_atoms, ref_water_residue,
                                                     center, sphere_radius)

        sum_water_removed += nb_water_removed

        output_atoms = renumber(tmp_atoms)

    return (output_atoms, sum_water_removed)


def main():
    """
    Run everythin from the command line
    """

    # Command line parsing
    args = define_options(sys.argv[1:])


    try:
        title, atoms, box = groIO.parse_file(args.filin)
    except groIO.FormatError:
        print(("Something is wrong in the format of {0}").format(args.filin),
               file=sys.stderr)
        return 1

    if args.verbose:
        print("The input coordinate file is {0}".format(args.filin))
        print("The output file will be {0}".format(args.filout))
        print("The reference atom for the lipid bilayer is {0}".format(args.lipid_atom))
        print("The reference atom for the water is {0}".format(args.water_atom))


    try:
        if args.sphere is None:
            output_atoms, nb_water = perform_bilayer_removing(atoms, args.lipid_atom,
                                                        args.water_atom, args.verbose)
        else:
            output_atoms, nb_water = perform_sphere_removing(atoms, args.sphere, args.radius,
                                                       args.water_residue, args.verbose)

        if args.verbose:
            print()
            print("The old system contained {0} atoms.".format(len(atoms)))
            print("{0} water molecules have been removed.".format(nb_water))
            print("The new system contains {0} atoms.".format(len(output_atoms)))

        # Write in the output file
        with open(args.filout, "w") as fout:
            for line in groIO.write_gro(title, output_atoms, box):
                print(line, end='', file=fout)

    except RefObjectError, e:
        print(e)
        return 1

    return 0

if __name__ == '__main__':
    sys.exit(main())
