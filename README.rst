g_remove_water
===========

.. image:: https://travis-ci.org/HubLot/g_remove_water.svg?branch=master
   :target: https://travis-ci.org/HubLot/g_remove_water
.. image:: https://coveralls.io/repos/HubLot/g_remove_water/badge.svg?branch=master
   :target: https://coveralls.io/r/HubLot/g_remove_water?branch=master 



Takes a gro file containing a bilayer and remove any water
molecules inside the bilayer.

To remove molecules, it's based on Z mean of the upper and lower leaflet.

At the end, the residues and the atoms are re-numbered and the number of atoms at the 2nd line is updated.
Title and box vectors are kept.

The atom to determine the upper and lower leaflet (typically the phosphate) can be set with the lipid_atom option (default = P1)

The atom to determine the water molecule can be set with --water_atom option.
For all-atoms files, you have to use the oxygen atom as a reference.

Water molecules can also be removed if they are in a sphere centered on the
geometrical center of the atom with a given set of residue name. Use --sphere
to set the list of residue names to use for the center calculation and
--radius to set the radius of the sphere.

Usage
-----
g_remove_water.py [-h] -f coord.gro [-o out.gro] [--water_atom OW] [(--lipid_atom P1 | --sphere RES [RES ...] --radius RADIUS)]

required arguments:
    -f FILE             The Structure file newly solvated (.gro)

optional arguments:
    -h, --help                     Show this help message and exit
    -o FILE                        The Output file (.gro)
    --lipid_atom P1                The reference atom for the bilayer (P1 by default)
    --water_atom OW                The reference atom for the water. Use the oxygen.
                                   (OW by default)
    --sphere RESNAME [RESNAME ...] Remove water molecules if they are in a
                                   sphere centered on the geometric center of
                                   atoms with the given residue names. You
                                   need the --radius option to be filled.
    --radius RADIUS                Remove water molecules if they are in a sphere of this
                                   radius centered on a given set of residue names. You
                                   need the --sphere option to be set.

Run tests
---------

Unit tests are available for g_remove_water in test_g_remove_water.py. You can
run them by simply execute test_g_remove_water.py::

    python test_g_remove_water.py

The `nosetests python module <https://nose.readthedocs.org>`_ allows a smarter
display of the report by displaying test function outputs only when the test
fail. If you have nosetests installed you can run the test by typing::

    nostests test_g_remove_water.py

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

