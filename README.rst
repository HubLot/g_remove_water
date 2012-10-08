g_remove_water
===========

Takes a gro file containing a bilayer and remove any water
molecules inside the bilayer.

To remove molecules, it's based on Z mean of the upper and lower leaflet.
At the end, the residus and the atoms are re-numbered and the number of atoms at the 2nd line is updated.
Title and box vectors are kept.

The atom to determine the upper and lower leaflet (typically the phosphate) can be set with the option --refatom (default = P1)



Usage
-----
g_remove_water.py [-h] -f FILIN [-o FILOUT] [--refatom REFATOM]

required arguments:
  -f FILIN           The Structure file newly solvated (.gro)

optional arguments:
  -h, --help         show this help message and exit
  -o FILOUT          The Output file (.gro)
  --refatom REFATOM  The reference atom for the bilayer (P1 by default)


Licence
-------

This program is free software: you can redistribute it and/or modify  
it under the terms of the GNU General Public License as published by   
the Free Software Foundation, either version 3 of the License, or      
(at your option) any later version.                                    
                                                                      
This program is distributed in the hope that it will be useful,        
but WITHOUT ANY WARRANTY; without even the implied warranty of         
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
GNU General Public License for more details.                           
                                                                          
A copy of the GNU General Public License is available at
http://www.gnu.org/licenses/gpl-3.0.html.

