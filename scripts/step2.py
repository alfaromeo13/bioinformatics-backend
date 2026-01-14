from multiprocessing.sharedctypes import Value
from collections import defaultdict
import subprocess
import sys

class Step2:
    def __init__(self, mut_chid, pro_mut_chid, resn, wtres, charmm_dir, thread_num):
        self.mut_chid = mut_chid
        self.pro_mut_chid = pro_mut_chid
        self.resn = resn
        self.wtres = wtres
        self.charmm_dir = charmm_dir
        self.thread_num = thread_num
        self.A1 = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                   'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HSD': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',
                   'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                   'TYR': 'Y', 'VAL': 'V'}

    #ne dela nič posebnega, samo pogleda če si pravilno vnesel ime aminokisline, da ni ravno ASZ namesto ASP
    def mutation_check(self):
        if self.wtres in self.A1:
            print("Mutation input listed correctly {} {} {}".format(self.wtres, self.mut_chid, self.resn))
        else:
            print("Mutation input {} {} {} is incorrect - Please check the mutation and the coresponding PDB file".format(self.wtres, self.mut_chid, self.resn))
            exit

    #Pripravi charmm input za mutiranje
    def charmm_inp_maker(self):
        resn = int(self.resn)
        resl = resn - 1
        resd = resn + 1
   
        if self.wtres in self.A1:
            wt = self.A1[self.wtres]
       
        for value in self.A1:
            print("!mutate to {}".format(self.A1[value]))
            print("open read unit 10 card name {}.psf".format(self.pro_mut_chid))
            print("read psf unit 10 card")
            print("open read unit 10 card name {}.crd".format(self.pro_mut_chid))
            print("read coor unit 10 card\n")
            #ukaz patch v charmmu npr. XA2G alanin v glicin, X ne označuje nič, vendar mora biti prisoten zaradi štiričrkovne omejitve komand
            print("patch X{}2{} {} {} {} {} {} {} setup ".format(wt, self.A1[value], self.pro_mut_chid, resn, self.pro_mut_chid, resl, self.pro_mut_chid, resd))
            print("ic build")
            for key, values in self.A1.items():
                if self.A1[value] in values:
                    print("rena resn {} sele segi {} .and. resid {} end".format(key, self.pro_mut_chid, resn))
            print("autogenerate angl dihe\n")
            print("open write unit 11 card name {}_{}_{}2{}.psf".format(self.pro_mut_chid, resn, wt, self.A1[value]))
            print("write psf unit 11 card")
            print("open write unit 11 card name {}_{}_{}2{}.crd".format(self.pro_mut_chid, resn, wt, self.A1[value]))
            print("write coor unit 11 card")
            print("open write unit 11 card name {}_{}_{}2{}.pdb".format(self.pro_mut_chid, resn, wt, self.A1[value]))
            print("write coor pdb unit 11 official\n")

    def generate_charmm_input(self):
        with open('step2_{}{}.inp'.format(self.mut_chid, self.resn), 'w') as f:
            sys.stdout = f

            print("stream charmm_dat/toppar.str")
            print("read rtf card append name charmm_dat/all_pmuta.rtf\n")

            self.charmm_inp_maker()

            print("stop")
            sys.stdout = sys.__stdout__

    def run_step2(self):
        def step2_error_check(mut_chid, resn, error_message):
            print("Error occurred during Step 2 for chain {} residue {}: subprocess returned non-zero exit status 1.".format(mut_chid, resn))
            print("Python error message: {}".format(error_message))
            print("Please check step2_{}{}.out for more information.".format(mut_chid, resn))
        
        print("Mutating the amino-acid {} on residue number {} on chain {}".format(self.wtres, self.resn, self.mut_chid))
        # V ukazu najdemo mpirun in število jeder, kar je lahko v prihodnosti spremenljivka
        command = f'mpirun -n {self.thread_num} {self.charmm_dir} -i step2_{self.mut_chid}{self.resn}.inp > step2_{self.mut_chid}{self.resn}.out'
        try:
            subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            if e.returncode == 1:
                step2_error_check(self.mut_chid, self.resn, str(e))
                sys.exit(1)



    

        

