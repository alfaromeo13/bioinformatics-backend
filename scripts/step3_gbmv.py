import os
import subprocess
import sys
import time
import threading
from collections import defaultdict


class Step3GBMV:
    def __init__(self, cgenff_name, complex_abnr, pro_mut_chid, resn, wtres,protein_chains_list, ligand_chains_list,charmm_dir, pdb_name,thread_num):                 
        self.cgenff_name = cgenff_name
        self.pdb_name = pdb_name
        self.complex_abnr = complex_abnr
        self.pro_mut_chid = pro_mut_chid
        self.resn = resn
        self.wtres = wtres
        self.protein_chains_list = protein_chains_list
        self.ligand_chains_list = ligand_chains_list
        self.thread_num = thread_num
        self.A1 = {'ALA': 'A','ARG': 'R',  'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                   'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HSD': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M',
                   'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                   'TYR': 'Y', 'VAL': 'V'}
        self.all_chains = []
        self.mutpro_chains = []
        self.intpar_chains = []
        self.non_mut_chid = []
        self.charmm_dir = charmm_dir
        self.subprocess_failed = False

    #Za sisteme, ki imajo custom ligande in user vnese parametre iz CGENFF    
    def ligand_rtf_prm(self):
        if self.cgenff_name is not None:
            print("open read card unit 10 name ligand.rtf")
            print("read rtf card unit 10 append")
            print("open read card unit 20 name ligand.prm")
            print("read para flex card unit 20 append")

    #Print za implicitno vodo, kliče se za vsak eksperiment, protein, ligand
    #Identični ukaz kot v članku (Kralj et. al., 2021, Frontiers in Chemistry)
    def implicit_water(self):
            print("GBMV TOL 1E-10 MEM 20 CUTA 20 DN 1.0 BUFR 0.2 -")
            print("EPSILON 80 BETA -12 SHIFT -0.1 SLOPE 0.9 -")
            print("LAMBDA1 0.5 P1 0.45 P2 1.25 P3 0.65 P6 8 -")
            print("ONX 1.9 OFFX 2.1 CORR 1 ALFRQ 1 -")
            print("SON 1.2 SOFF 1.5 -")
            print("FAST 1 SGBFRQ 4 SXD 0.3 -")
            print("WTYP 1 NPHI 5 SA 0.00542 SB 0.9")


############
    #tole je za kompleks GBMV

    def generate_all_chains(self):
        for chain in self.protein_chains_list:
            self.all_chains.append('PRO' + chain)
        for chain in self.ligand_chains_list:
            self.all_chains.append('PRO' + chain)

    def generate_non_mut_chid(self):
        mut_chid = self.pro_mut_chid
        self.non_mut_chid = [chid for chid in self.all_chains if chid != mut_chid]
        for i in range(len(self.non_mut_chid)):
            self.non_mut_chid[i] = 'pro' + self.non_mut_chid[i]

    def process_non_mut(self):
        for chain_name in self.non_mut_chid:
            non_mut_chain_name = chain_name[3:]

            print("read psf card append name {}.psf".format(non_mut_chain_name.lower()))            
            print("read coor card append name {}.crd".format(non_mut_chain_name.lower()))              
 
    #Bere komplekse ki so nastali pri INTE ukazu in pripravi za gbmv 
    def gbmv_comp_read(self, wt):
            for value in self.A1.values():
                #self.process_non_mut()
                mut_chid = self.pro_mut_chid[3:]
                print("read psf card append name joined_{}_{}_{}2{}.psf".format(self.pro_mut_chid, self.resn, wt, value))
                print("read coor card append name joined_{}_{}_{}2{}.crd".format(self.pro_mut_chid, self.resn, wt, value))

                # MINIMIZATION zaenkrat ne rabimo ker se bere fajl od inte ki je minimiziran
                print("set t {}_{}_{}2{}\n".format(self.pro_mut_chid, self.resn, wt, value))
                #dobimo že minimizirano strukturo iz inte ukaza, nujno je da so isti pogoji minimizacije
                #Tukaj izračunamo kompleks oz. comp
                self.implicit_water()
                print("ENERGY atom CDIE eps 1 cutnb 21 ctofnb 18 ctonnb 16 switch vswitch\n")
                print("open write unit 21 form name comp_gbmv_{}{}{}.dat append".format(self.wtres, self.resn, mut_chid).lower())
                print("write title unit 21")
                print("*@t      ?ener    ")
                print("*\n")
                print("close unit 21")
                #Tukaj more biti boljša funkcija če je več chainov tako da iz protein chains lista
                for chain_name in self.protein_chains_list:
                    non_mut_chain_name = chain_name.lower()
                    print("delete atom select segid pro{} end".format(non_mut_chain_name))  
                print("gbmv clear")
                self.implicit_water()
                #ligand
                print("ENERGY atom CDIE eps 1 cutnb 21 ctofnb 18 ctonnb 16 switch vswitch\n")
                print("open write unit 31 form name liga_gbmv_{}{}{}.dat append".format(self.wtres, self.resn, mut_chid).lower())
                print("write title unit 31")
                print("*@t      ?ener    ")
                print("*\n")
                print("close unit 31")
                print("delete atom select all end")
                print("gbmv clear")
                #Protein
                print("read psf card append name joined_{}_{}_{}2{}.psf".format(self.pro_mut_chid, self.resn, wt, value))
                print("read coor card append name joined_{}_{}_{}2{}.crd".format(self.pro_mut_chid, self.resn, wt, value))
                #Zbrišemo ligand chaine
                for chain_name in self.ligand_chains_list:
                    non_mut_chain_name = chain_name.lower()
                    print("delete atom select segid pro{} end".format(non_mut_chain_name))  
                self.implicit_water()
                print("ENERGY atom CDIE eps 1 cutnb 21 ctofnb 18 ctonnb 16 switch vswitch\n")
                print("open write unit 41 form name prot_gbmv_{}{}{}.dat append".format(self.wtres, self.resn, mut_chid).lower())
                print("write title unit 41")
                print("*@t      ?ener    ")
                print("*\n")
                print("close unit 41")
                print("delete atom select all end")
                print("gbmv clear")

    def create_step3_gbmv_comp(self):
            mut_chid = self.pro_mut_chid[3:]
            with open('step3_{}{}{}_gbmv.inp'.format(self.wtres, self.resn, mut_chid).lower(), "w") as f:
                sys.stdout = f

                print("*GENERATED BY SEBALAB")
                print("*SMEHeC: Swift Mutational Energy Heatmap Calculator")
                print("*")
                print("stream charmm_dat/toppar.str")
                self.ligand_rtf_prm()
                mut_chid = self.pro_mut_chid[3:]
                print("open write unit 21 form name comp_gbmv_{}{}{}.dat".format(self.wtres, self.resn, mut_chid).lower())
                print("write title unit 21\n")
                print("open write unit 31 form name liga_gbmv_{}{}{}.dat".format(self.wtres, self.resn, mut_chid).lower())
                print("write title unit 31\n")
                print("open write unit 41 form name prot_gbmv_{}{}{}.dat".format(self.wtres, self.resn, mut_chid).lower())
                print("write title unit 41\n")
                
                wt = self.A1.get(self.wtres)
                self.gbmv_comp_read(wt)

                print("stop")
                sys.stdout = sys.__stdout__

    def run_step3_gbmv_comp(self):
            mut_chid = self.pro_mut_chid[3:].lower()
            mut_chid_caps = self.pro_mut_chid[3:]
            wtres = self.wtres.lower()
            resn = self.resn

            print("Calculating GBMV energy for mutation {} {} {} ".format(self.wtres, mut_chid_caps, self.resn))
            command = f'mpirun -n {self.thread_num} {self.charmm_dir} -i step3_{wtres}{resn}{mut_chid}_gbmv.inp > step3_gbmv_{wtres}{resn}{mut_chid}.out'
            inter_ener_dat_path = "comp_gbmv_{}{}{}.dat".format(self.wtres, self.resn, mut_chid).lower()
            total_lines = 20
            progress_thread = threading.Thread(target=self.print_progress, args=(inter_ener_dat_path, total_lines))
            progress_thread.start()

            try:
                result = subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as e:
                print(f"\nError occurred during Step 3: subprocess returned non-zero exit status {e.returncode}.")
                print(f"Python error message: {str(e)}")
                print(f"Please check step3_{wtres}{resn}{mut_chid}.out for more information.")
                self.subprocess_failed = True 
                progress_thread.join()
                sys.exit(e.returncode)

            self.print_progress_exit = True
            progress_thread.join()

            print("\n")    

    def print_progress(self, inter_ener_dat_path, total_lines):
        self.print_progress_exit = False
        progress = 0
        while progress < total_lines and not self.print_progress_exit and not self.subprocess_failed:
            time.sleep(1)  # Wait 1 second
            line_count = self.count_lines(inter_ener_dat_path)
            progress = min(line_count, total_lines)
            print(f"\rCalculating progress: {progress}/20        ", end='', flush=True)

        if self.subprocess_failed:
            sys.exit(2)

    def count_lines(self, inter_ener_dat_path):
        line_count = 0
        while not os.path.exists(inter_ener_dat_path):
            time.sleep(1)  # čaka eno sekundo za refresh
        with open(inter_ener_dat_path, 'r') as inter_ener_dat:
            for line in inter_ener_dat:
                if line.startswith("PRO"):
                    line_count += 1
                    if line_count == 20:
                        return line_count
        return line_count

#####################
#Za združit rezultate
    def calculate_gbmv(self):
        mut_chid = self.pro_mut_chid[3:].lower()
        mut_chain_id = self.pro_mut_chid[3:]
        complex_file = open('comp_gbmv_{}{}{}.dat'.format(self.wtres, self.resn, mut_chid).lower(), 'r')
        ligand_file = open('liga_gbmv_{}{}{}.dat'.format(self.wtres, self.resn, mut_chid).lower(), 'r')
        protein_file = open("prot_gbmv_{}{}{}.dat".format(self.wtres, self.resn, mut_chid).lower(), 'r')

        complex_lines = [line for line in complex_file if line.startswith("PRO")]
        ligand_lines = [line for line in ligand_file if line.startswith("PRO")]
        protein_lines = [line for line in protein_file if line.startswith("PRO")]
        
        complex_file.close()
        ligand_file.close()
        protein_file.close()

        # Prevermo če mamo res 20 v obeh
        if len(complex_lines) != len(protein_lines):
            print("Error: The number of lines in complex and apo_protein files do not match.")
            return

        output_file_path = 'gbmv_ener_{}{}{}.dat'.format(self.wtres, self.resn, mut_chid).lower()
        output_file = open(output_file_path, 'w')

        # Shranimo rezultate outputa
        for i in range(len(complex_lines)):
            complex_line = complex_lines[i]
            protein_line = protein_lines[i]
            ligand_line = ligand_lines[i]


            # Za dobit energijo iz vrednosti
            complex_energy = float(complex_line.split()[1])
            protein_energy = float(protein_line.split()[1])
            ligand_energy = float(ligand_line.split()[1])

            result = complex_energy - (protein_energy + ligand_energy)
            result_line = f"{complex_line.split()[0]:<20} {result:.2f}"
            output_file.write(result_line + '\n')
        output_file.close()
