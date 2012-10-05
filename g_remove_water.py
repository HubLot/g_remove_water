#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Script which removes water molecules inside a membrane (after doing genbox for solvation for example)"""

__author__ = "Marc Gueroult & Hubert Santuz"
__version__ = "1.0"

import sys
import os
import argparse


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


def define_options():
    """Define the script options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", action="store", type=isfile, dest = "filin",required=True,
            help="The Structure file newly solvated (.gro) REQUIRED")
    parser.add_argument("-o", action="store", type=str, dest = "filout",
            help="The Output file (.gro)",default="out.gro")
    parser.add_argument("--refatom", action="store", type=str, dest = "refatom",
            help="The reference atom for the bilayer", default="P1")
    return parser        

def Zmean_values(lines,atm_ref):
    """Return the Z mean for the upper and lower leaflet"""

    Zlower = Zupper = 0.0

    #Get the global Z mean
    Zsum = 0.0
    nb_atoms = nb_lower = nb_upper = 0
    for i in lines:
        atom_name = i[10:15].strip()
        if atom_name == atm_ref:
            Zsum += float(i[37:44].strip()) # Z coordinate
            nb_atoms += 1
    Zmean = Zsum / float(nb_atoms)


    #Z mean for the upper and lower leaflet
    for i in lines:
        atom_name = i[10:15].strip()
        if atom_name == atm_ref:
            z = float(i[37:44].strip()) # Z coordinate
            if z < Zmean:
                Zlower += z
                nb_lower +=1
            else:
                Zupper += z
                nb_upper +=1

    Zlower = Zlower/float(nb_lower)
    Zupper = Zupper/float(nb_upper)

    return Zlower,Zupper


def remove_water(lines, Zupper, Zlower):
    """Return  a list containing all lines of the gro file to keep"""


    reslines = [] #Lines which will be write in the output file
    #Appends the first 2 lines (ie Title and Number of atoms)
    reslines.append(lines[0].strip()) 
    reslines.append(lines[1].strip())

    resnumber = -1
    nb_atoms = 0 #Total number of atoms
    nb_res_water = 0 #Number of waters molecules removed

    for i in lines[2:-1]: #Discard for 2 lines and the last one
        to_store = True  #Does we store the line in the result ?
        atom_name = i[10:15].strip()
        ref_resnumber = int(i[0:5]) #Residue number of the current atoms
        if atom_name == "OW":               #Oxygen of the water molecule
            z = float(i[37:44])
            if z > Zlower and z < Zupper:    # Water in bilayer
                resnumber = int(i[0:5])     # Store the residue number of the water to removed
                to_store = False            #Don't store the line
                nb_res_water +=1

        #If the ref resnumer and the resnumer are equal means we don't store the hydrogens of the previous water molecule
        if to_store and ref_resnumber != resnumber:
            nb_atoms += 1 # 1 line = 1 atom
            reslines.append(i[:-1])


    #Update the final number of atoms
    reslines[1] = str(nb_atoms)
    #Append the last line of the gro file (Box vectors)
    reslines.append(lines[-1])

    return reslines, nb_res_water


def find_ref_atom(lines,ref_atom):
    """Return a boolean wether the "atom" was find in the file"""

    for i in lines:
        atom_name =  i[10:15].strip()
        if ref_atom == atom_name:
            return True
    return False


if __name__ == '__main__' :

    #Command line parsing
    parser = define_options()
    args = parser.parse_args()
    
    filin = args.filin
    filout = args.filout
    refatom = args.refatom


    print "The input coordinate file is {0}".format(filin)
    print "The output file will be {0}".format(filout)
    print "The reference atom for the lipid bilayer is {0}".format(refatom)

    f=open(filin,'r')
    data=f.readlines()
    f.close()

    print
    print "Checking the reference atom...",
    if not find_ref_atom(data,refatom):
        print "Oops!"
        print "The reference atom {0} for the bilayer was not find. Exiting...".format(refatom)
        sys.exit()
    print "Done!"
    
    #Get Z mean for  the upper and lower leaflet
    Zlower, Zupper = Zmean_values(data, refatom)

    print "Removing water inside the bilayer...",
    #Remove water molecules inside the bilayer
    output,water_removed = remove_water(data,Zupper,Zlower)
    print "Done!"

    print
    print "The old system contained {0} atoms.".format(data[1].strip())
    print "{0} water molecules have been removed.".format(water_removed)
    print "The new system contains {0} atoms.".format(output[1])

    #Write in the output file
    f=open(filout,'w')
    f.write('\n'.join(output))
    f.close()



