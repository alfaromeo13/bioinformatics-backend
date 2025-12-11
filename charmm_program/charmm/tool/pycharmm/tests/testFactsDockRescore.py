# pycharmm rigid cdocker test case

## Import module
import numpy as np
from pycharmm import *

################################################################
##
##		Begin of pyCHARMM Rigid CDOCKER
##
################################################################

## Topology and parameter files
settings.set_bomb_level(-1)
read.rtf('toppar/top_all36_prot.rtf')
read.rtf('toppar/top_all36_cgenff.rtf', append = True)
read.prm('toppar/par_all36m_prot.prm', flex = True)
read.prm('toppar/par_all36_cgenff.prm', append = True, flex = True)
settings.set_bomb_level(0)
charmm_script('stream data/benzene.rtf')

## Build system
read.psf_card("data/t4.psf", append = True)
read.pdb("data/t4.pdb", resid = True)
read.sequence_pdb("data/benzene.pdb")
gen.new_segment(seg_name = "LIGA")
read.pdb("data/benzene.pdb", resid = True)

## Prepare for minimization with FACTS
ligand = SelectAtoms().by_seg_id("LIGA")
receptor = ligand.__invert__()

## FACTS rescoring

facts_ener = cdocker.FACTS_rescore(fixAtomSel = receptor, steps = 100)

# Now translate ligand by 500 A away from pocket to get the ligand+receptor
pos = coor.get_positions()
pos.x += 500*np.array(ligand)
coor.set_positions(pos)

# Don't use any minimization
unbound = facts_ener
unbound -= cdocker.FACTS_rescore(fixAtomSel = receptor, steps = 0) 

if abs(facts_ener - (-4419.27931)) <= 0.01 :
    print("Testcase receptor FACTS result: PASS")
else:
    print("Testcase receptor FACTS result: FAIL")

if abs(unbound - (-15.60511)) <= 0.01 :
    print("Testcase binding FACTS result: PASS")
else:
    print("Testcase binding FACTS result: FAIL")
