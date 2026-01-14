from collections import defaultdict
import os
import re
import subprocess
import sys

class Step1:
    def __init__(self, pdb_name, charmm_dir):
        self.pdb_name = pdb_name
        self.charmm_dir = charmm_dir
        self.gly_1st = False

    def sequence_num(self, chain_pdb):
        with open(chain_pdb, 'r') as file:
            for line in file:
                if line.startswith('ATOM'):
                    seq_num_str = line[22:26].strip()
                    seq_num = int(seq_num_str) if seq_num_str.isdigit() else 0
                    offset = -(seq_num - 1)
                    print('!first sequence number is {}'.format(seq_num))
                    print('read coor pdb unit 3 offset {}'.format(offset))
                    break


    #def check_gly_1st                

    #Spodnje funkcije check preverijo kakšen tip atomov se nahaja v pdb datoteki
    #To nareka kakšen tip CHARMM kode je potrebno uporabiti predvsem je razlika pri gener komandi
    def check_aa(self, line):
        residue_names = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HSD', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        if not self.gly_1st and any(residue_name in line[17:20] for residue_name in residue_names):
            if residue_names.index(line[17:20].strip()) == 7:  # Check if the residue is GLY
                self.gly_1st = True
        return any(residue_name in line for residue_name in residue_names)

    
    def check_ions(self, line):
        ions = ['LIT', 'SOD', 'MG', 'POT', 'CAL', 'RUB', 'CES', 'BAR', 'ZN', 'CAD', 'CLA']
        return any(ion in line for ion in ions)
    
    def check_dna(self, line, chain_letter_up):
        nucleotides_dna = ['ADE', 'CYT', 'GUA', 'THY']
        has_line_code = False

        with open(self.pdb_name, 'r') as file:
            for pdb_line in file:
                chain_id = pdb_line[21:22].strip()
                nuc_name = pdb_line[17:20].strip()
                if chain_id == chain_letter_up and nuc_name in ['DA', 'DC', 'DG', 'DT']:
                    has_line_code = True
                    break

        return has_line_code and any(nuc in line for nuc in nucleotides_dna)
    
    def check_rna(self, line, chain_letter_up):
        nucleotides_rna = ['ADE', 'CYT', 'GUA', 'URA']
        has_line_code = False

        with open(self.pdb_name, 'r') as file:
            for pdb_line in file:
                chain_id = pdb_line[21:22].strip()
                nuc_name = pdb_line[17:20].strip()
                if chain_id == chain_letter_up and nuc_name in ['A', 'C', 'G', 'U']:
                    has_line_code = True
                    break

        return has_line_code and any(nuc in line for nuc in nucleotides_rna)

    #Ta funkcija za vsak posamezni pdb chain generira input datoteko, ki omogoča nastanek rtf in prm 
    def step1_maker(self):
        pdb_files = [file for file in os.listdir() if file.endswith('.pdb') and len(file) == 5]

        for pdb_file in pdb_files:
                
            with open(pdb_file, 'r') as pdb:
                lines = pdb.readlines()
            #To je za special case če je HEM
            with open(pdb_file, 'w') as pdb:
                for line in lines:
                    pdb.write(line.replace('HEM ', 'HEME'))
            
            chain_letter = pdb_file[0]
            chain_letter_up = pdb_file[0].upper()
            with open('step1_' + chain_letter + '.inp', 'w') as f:
                sys.stdout = f

                with open(pdb_file, 'r') as pdb:
                    lines = pdb.readlines()
                    has_amino_acids = any(self.check_aa(line) for line in lines)
                    has_ions = any(self.check_ions(line) for line in lines)
                    has_dna = any(self.check_dna(line,chain_letter_up) for line in lines)
                    has_rna = any(self.check_rna(line,chain_letter_up) for line in lines)

                if has_amino_acids:  #CHARMM ukazi za amino-kisline razen če je GLY prvi se drugače patcha
                    print('bomblev -1')
                    print('stream charmm_dat/toppar.str')
                    print('open unit 3 read card name {}.pdb'.format(chain_letter))
                    print('read sequ pdb unit 3')
                    print('rewind unit 3')
                    if has_amino_acids:
                        if self.gly_1st:
                            print('gener pro{} setup warn first GLYP last CTER'.format(chain_letter))
                        else:
                            print('gener pro{} setup warn first NTER last CTER'.format(chain_letter))
                    if self.gly_1st:
                        self.gly_1st = False        
                    #print('gener pro{} setup warn first NTER last CTER'.format(chain_letter))
                    chain_pdb = chain_letter + '.pdb'
                    self.sequence_num(chain_pdb)
                    print('close unit 3\n')
                    print('!place missing heavy atoms')
                    print('ic purge')
                    print('ic param')
                    print('ic fill preserve')
                    print('ic build')
                    print('define test sele segid pro{} .not. type H* .and. .not. init show end\n'.format(chain_letter))
                    print('!naredimo vse koordinate za vodike')
                    print('coor init sele type H* end')
                    print('hbuild sele type H* end')
                    print('define test sele segid pro{} .not. init show end\n'.format(chain_letter))
                    print('open unit 12 write card name pro{}.psf'.format(chain_letter))
                    print('write psf card unit 12')
                    print('open unit 12 write card name pro{}.crd'.format(chain_letter))
                    print('write coor card unit 12\n')
                    print('stop')
                elif has_ions:  #CHARMM ukazi za ione
                    print('bomblev -1')
                    print('stream charmm_dat/toppar.str')
                    print('open unit 3 read card name {}.pdb'.format(chain_letter))
                    print('read sequ pdb unit 3')
                    print('rewind unit 3')
                    print('gener pro{} setup warn first none last none '.format(chain_letter))
                    chain_pdb = chain_letter + '.pdb'
                    self.sequence_num(chain_pdb)
                    print('close unit 3\n')
                    print('open unit 12 write card name pro{}.psf'.format(chain_letter))
                    print('write psf card unit 12')
                    print('open unit 12 write card name pro{}.crd'.format(chain_letter))
                    print('write coor card unit 12\n')
                    print('stop')
                elif has_dna: #CHARMM ukazi za DNa, ki ima posebno patchanje z deo5
                    print('bomblev -1')
                    print('stream charmm_dat/toppar.str')
                    print('open unit 3 read card name {}.pdb'.format(chain_letter))
                    print('read sequ pdb unit 3')
                    print('gener pro{} setup warn first 5TER last 3TER'.format(chain_letter))
                    chain_pdb = chain_letter + '.pdb'
                    self.sequence_num(chain_pdb)
                    print('close unit 3\n')
                    print('!process DNA chain')
                    print('coor stat sele segid pro{} end'.format(chain_letter))
                    print('calc fres = ?selires')
                    print('calc nres = ?nres')
                    print('set k = @fres')
                    print(f'label DNA{chain_letter}deox')
                    print('coor stat sele segid pro{} .and. ires @k end'.format(chain_letter))
                    print('set resname ?selresn')
                    print('set kthresid = ?selresi')
                    print('set patch deox ! for pyrimidine')
                    print('if k .eq. @fres set patch = deo5')
                    print('patch @patch pro{} @kthresid setup warn'.format(chain_letter))
                    print('incr k by 1')
                    print(f'if k .le. @nres goto DNA{chain_letter}deox')
                    print('autogenerate angles dihedrals\n')
                    chain_pdb = chain_letter + '.pdb'
                    self.sequence_num(chain_pdb)
                    print('close unit 3\n')
                    print('!place missing heavy atoms')
                    print('ic purge')
                    print('ic param')
                    print('ic fill preserve')
                    print('ic build')
                    print('define test sele segid pro{} .not. type H* .and. .not. init show end\n'.format(chain_letter))
                    print('!naredimo vse koordinate za vodike')
                    print('coor init sele type H* end')
                    print('hbuild sele type H* end')
                    print('define test sele segid pro{} .not. init show end\n'.format(chain_letter))
                    print('open unit 12 write card name pro{}.psf'.format(chain_letter))
                    print('write psf card unit 12')
                    print('open unit 12 write card name pro{}.crd'.format(chain_letter))
                    print('write coor card unit 12\n')
                    print('stop')
                elif has_rna:  #RNA ukazi, ki so podobni kot DNA vendar ne potrebujejo dodatnega patchanja   
                    print('bomblev -1')
                    print('stream charmm_dat/toppar.str')
                    print('open unit 3 read card name {}.pdb'.format(chain_letter))
                    print('read sequ pdb unit 3')
                    print('rewind unit 3')
                    print('!process RNA chain')
                    print('gener pro{} setup warn first 5TER last 3TER'.format(chain_letter))
                    chain_pdb = chain_letter + '.pdb'
                    self.sequence_num(chain_pdb)
                    print('close unit 3\n')
                    print('!place missing heavy atoms')
                    print('ic purge')
                    print('ic param')
                    print('ic fill preserve')
                    print('ic build')
                    print('define test sele segid pro{} .not. type H* .and. .not. init show end\n'.format(chain_letter))
                    print('!naredimo vse koordinate za vodike')
                    print('coor init sele type H* end')
                    print('hbuild sele type H* end')
                    print('define test sele segid pro{} .not. init show end\n'.format(chain_letter))
                    print('open unit 12 write card name pro{}.psf'.format(chain_letter))
                    print('write psf card unit 12')
                    print('open unit 12 write card name pro{}.crd'.format(chain_letter))
                    print('write coor card unit 12\n')
                    print('stop')
                else: #Print v primeru da imamo custom male ligande za katere smo podatke dobili z cgenff, nanaša se na ligand_prep 
                    #Kliče namreč ligand.rtf in ligand.prm
                    print('bomblev -1')
                    print('stream charmm_dat/toppar.str')
                    print('!Custom topology and parameter files for lig')
                    print('open read card unit 10 name ligand.rtf')
                    print('read rtf card unit 10 append')
                    print('open read card unit 20 name ligand.prm')
                    print('read para flex card unit 20 append')
                    print('open unit 3 read card name {}.pdb'.format(chain_letter))
                    print('read sequ pdb unit 3')
                    print('rewind unit 3')
                    print('gener pro{} setup warn first none last none'.format(chain_letter))
                    chain_pdb = chain_letter + '.pdb'
                    self.sequence_num(chain_pdb)
                    print('close unit 3\n')
                    print('open unit 12 write card name pro{}.psf'.format(chain_letter))
                    print('write psf card unit 12')
                    print('open unit 12 write card name pro{}.crd'.format(chain_letter))
                    print('write coor card unit 12\n')
                    print('stop')

                sys.stdout = sys.__stdout__

    #Ukaz pokliče step ena in ga zažene
    def run_step1(self):
        def step1_error_check(chain_letter, error_message):
            print(f"Error occurred during Step 1 for chain {chain_letter}: subprocess returned non-zero exit status 1.")
            print(f"Python error message: {error_message}")
            print(f"Please check step1_{chain_letter}.out for more information.")

        pdb_files = [file for file in os.listdir() if file.endswith('.pdb') and len(file) == 5]
        for pdb_file in pdb_files:
            chain_letter = pdb_file[0]
            charmm_directory = self.charmm_dir
            ch_letter = chain_letter.upper()
            print(f"Generating CRD and PSF for chain {ch_letter}")
            command = f'{charmm_directory} -i step1_{chain_letter}.inp > step1_{chain_letter}.out'
            try:
                subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except subprocess.CalledProcessError as e:
                if e.returncode == 1:
                    step1_error_check(chain_letter, str(e))
                    sys.exit(1)
            #Error checker za koordinate, 9999.0 je dodeljeno avtomatsko če manjkajo, posledica PDB napak        
            crd_file = f"pro{chain_letter}.crd"
            with open(crd_file, 'r') as file:
                for line_num, line in enumerate(file, start=1):
                    if re.search(r'\b9999\b', line):
                        print(f"Coordinate issue in {crd_file}, line {line_num}: {line.strip()}")
                        sys.exit(1)

        


    


     

