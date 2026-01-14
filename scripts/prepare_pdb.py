from collections import defaultdict
import os
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

class PreparePdb:
    def __init__(self, pdb_name):
        self.pdb_name = pdb_name

    #Preimenuje vse HETATM v PDBju v ATOM, kar omogoča boljše procesiranje v CHARMMU
    def hetatm2atm(self):
        with open(self.pdb_name, 'r') as file:
            lines = file.readlines()

        with open(self.pdb_name, 'w') as file:
            for line in lines:
                file.write(line.replace("HETATM", "ATOM  "))
   
    #Preimenuje vse HIS v PDB v HSD, poleg tega charmm pozna še HSP in HSE vendar je to stvar protonacije
    #V trenutni verziji tako privzamemo najbolj pogosto protonacijo HSD
    def his2hsd(self):
        with open(self.pdb_name, 'r') as his, open('process_' + self.pdb_name, 'w') as hsd:
            for line in his:
                hsd.write(line.replace("HIS", "HSD"))

    #V kolikor so prisotni ioni je potrebno preimenovati tudi te, saj imajo nekateri drugačna imena kot so navedena v PDB 
    def replace_ion(self):  #te ioni imajo v charmm dokumentaciji drugačno ime
        replace_dict = {' ZN': 'ZN2', 'CAD': 'CD2', ' CA': 'CAL'}

        with open('process_'+self.pdb_name, 'r') as file:
            lines = file.readlines()

        with open('process_'+self.pdb_name, 'w') as file:
            for line in lines:
                if line[17:20] in replace_dict:
                    line = line[:17] + replace_dict[line[17:20]] + line[20:]
                file.write(line)

    #RNA in DNA imata drugačno ime v charmmu, v PDB so ali enočrkovne (RNA) ali dvočrkovne (DNA)
    def replace_dna_rna(self):  
        replace_dict = {' DA': 'ADE', ' DC': 'CYT', ' DG' : 'GUA', ' DT' : 'THY', '  A': 'ADE', '  C': 'CYT', '  G' : 'GUA', '  U' : 'URA' }

        with open('process_'+self.pdb_name, 'r') as file:
            lines = file.readlines()

        with open('process_'+self.pdb_name, 'w') as file:
            for line in lines:
                if line[17:20] in replace_dict:
                    line = line[:17] + replace_dict[line[17:20]] + line[20:]
                file.write(line)

    #Ker je za generacijo rtf in prm potrebno ločiti vse chaine, uporabimo iz biopython funkcijo za chain splitting
    # Če ima PDB chaine a, b,c bodo vsi atomi ki pripadajo določenemu chainu shranjeni potem kot a.pdb, b.pdb, c.pdb
    def chain_spliter(self):
        io = PDBIO()
        parser = PDBParser()
        structure = parser.get_structure('process_' + self.pdb_name, 'process_' + self.pdb_name)
        pdb_chains = structure.get_chains()
        for chain in pdb_chains:
            io.set_structure(chain)
            io.save(chain.get_id() + '.pdb')
        print()
    #Vse A.pdb oz druge chain pdb datoteke preimenuje v small caps za lažje delo
    def small_caps(self):
        directory = os.getcwd()
        mis = '/'
        for file in os.listdir(directory):
            if len(file) == 5 and file.endswith('.pdb'):
                os.rename(directory + mis + file, directory + mis + file.lower())
