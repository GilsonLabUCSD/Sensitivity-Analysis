#!/usr/bin/env python2
import os as os
import sys as sys
from parmed.tools import summary
from parmed.tools import strip
from parmed.tools import printInfo
from parmed.tools import printBonds
from parmed.tools import printAngles
from parmed.tools import printDihedrals
from parmed.tools import printLJTypes
from parmed.tools import parmout
from parmed.amber import AmberParm
from parmed.tools import writeFrcmod
from parmed.structure import Structure

if len(sys.argv) < 2:
    print("Args: <mol.topo> ")
    print("Aborted. Please specify the name of the topology file.")
    sys.exit()
elif len(sys.argv) > 2:
    print("Args: <mol.topo> ")
    print("Aborted. Too many arguments!")
    sys.exit()

topfile = sys.argv[1] 

if not os.path.isfile(topfile):
    print('The input file %s does not exist.' % topfile)
    sys.exit()

parm = AmberParm(topfile) 

bondedInfoFile = open('moleculeInfo','w')
action = summary(parm)
action.execute()
bondedInfoFile.write('%s\n' %action)  

# Print out force field parameters
action = writeFrcmod(parm, 'FFparameters.dat')
action.execute()

# Generate the "dry" topology file
action = strip(parm, ':WAT')
action.execute()

# Print out the connectivity (bond, angle, torsion) from the dry topology

bondedInfoFile.write('[ Atoms ]\n')
action = printLJTypes( parm );
action.execute()
bondedInfoFile.write('%s\n' %action)

bondedInfoFile.write('[ Bonds ]\n')
action = printBonds(parm)
action.execute()
bondedInfoFile.write('%s\n' %action)

bondedInfoFile.write('[ Angles ]\n')
action = printAngles(parm)
action.execute()
bondedInfoFile.write('%s\n' %action)

bondedInfoFile.write('[ Dihedrals ]\n')
action = printDihedrals(parm)
action.execute()
bondedInfoFile.write('%s\n' %action)

bondedInfoFile.close();
