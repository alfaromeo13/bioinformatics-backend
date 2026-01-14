# pycharmm  FACTS_score test case

## Import module
import numpy as np
from pycharmm import *
from pycharmm.cdocker import FACTS_rescore


## Read in the topology and parameter file
settings.set_bomb_level(-1)
read.rtf('toppar/top_all36_prot.rtf')
read.prm('toppar/par_all36m_prot.prm', flex = True)

## Build protein
read.sequence_string('ALA')

gen.new_segment('ADP', 'ACE', 'CT3', setup_ic=True)

ic.prm_fill(False)
ic.seed(1, 'CAY', 1, 'CY', 1, 'N')
ic.build()

## Loop over trials to ensure FACTS Clear works
xyz_0 = coor.get_positions()
E_trial = []
for i in range(10):
    coor.set_positions(xyz_0)
    E_trial.append(cdocker.FACTS_rescore(fixAtomSel = select_atoms.SelectAtoms(), steps = 100, tolgrd = 0.001))

E_trial = np.array(E_trial)
Error = np.sum(np.abs(E_trial - E_trial[0]))
if (Error < 1e-5): print('>>>>>>>>>>>>>>>>>>>>>>TEST 1 PASSED<<<<<<<<<<<<<<<<<<<')
else: print('>>>>>>>>>>>>>>>>>>>>>>TEST 1 FAILED<<<<<<<<<<<<<<<<<<<')

## Loop over trials with fixed atoms to ensure FACTS Clear works

# Choose random natom tuple of T/F
select = np.random.choice([True,False], size=psf.get_natom())
print(f'Number of fixed atoms = {psf.get_natom() - np.sum(select)}')

xyz_0 = coor.get_positions()
E_trial = []
for i in range(10):
    coor.set_positions(xyz_0)
    E_trial.append(cdocker.FACTS_rescore(fixAtomSel = select_atoms.SelectAtoms().set_selection(select), steps = 100, tolgrd = 0.001))

E_trial = np.array(E_trial)
Error = np.sum(np.abs(E_trial - E_trial[0]))
if (Error < 1e-5): print('>>>>>>>>>>>>>>>>>>>>>>TEST 2 PASSED<<<<<<<<<<<<<<<<<<<')
else: print('>>>>>>>>>>>>>>>>>>>>>>TEST 2 FAILED<<<<<<<<<<<<<<<<<<<')
