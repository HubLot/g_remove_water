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


__description__ = \
"""
g_remove_water.py

Takes a gro file containing a bilayer and remove any water molecules inside the bilayer.
To remove molecules, it's based on Z mean of the upper and lower leaflet.
At the end, the residus and the atoms are re-numbered and the number of atoms at the 2nd
line is updated.
Title and box vectors are kept.
The atom to determine the upper and lower leaflet (typically the phosphate) can be set with the option --refatom
"""

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
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument("-f", action="store", type=isfile, dest = "filin",required=True,
            help="The Structure file newly solvated (.gro) REQUIRED")
    parser.add_argument("-o", action="store", type=str, dest = "filout",
            help="The Output file (.gro)",default="out.gro")
    parser.add_argument("--refatom", action="store", type=str, dest = "refatom",
            help="The reference atom for the bilayer (P1 by default)", default="P1")
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


def renumber(lines,start_res=None):
    """Renum the atoms and the residus in a file"""

    out = []
    #Appends the first 2 lines (ie Title and Number of atoms)
    out.append(lines[0].strip()) 
    out.append(lines[1].strip()) 
 
    start = True
    output_atom = 0

    for i in lines[2:-1]: #Discard for 2 lines and the last one
        #If first residue, start the counter regarding the number of the starting residue
        if start:
            start = False
            current_res = i[0:5]
            if start_res:
                output_res = start_res
            else:
                output_res = 1

        #New residue
        if i[0:5] != current_res:
            if output_res == 99999: #cannot go to 100000
                output_res = 1
            else:
                output_res +=1
            current_res = i[0:5]

        if output_atom == 99999: #cannot go to 100000
            output_atom = 1
        else:   
            output_atom +=1

        #Append the line with the correct number of residu and atom
        out.append("%5i%s%5i%s" % (output_res,i[5:15],output_atom,i[20:]))


    #Save the last line (box vectors)
    out.append(lines[-1])

    return out


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
    temp_lines,water_removed = remove_water(data,Zupper,Zlower)
    print "Done!"


    first_res_number = int(data[3][0:5])
    print "Renumber residues and atoms...",
    output = renumber(temp_lines,first_res_number)
    print "Done!"

    print
    print "The old system contained {0} atoms.".format(data[1].strip())
    print "{0} water molecules have been removed.".format(water_removed)
    print "The new system contains {0} atoms.".format(output[1])

    #Write in the output file
    f=open(filout,'w')
    f.write('\n'.join(output))
    f.close()



