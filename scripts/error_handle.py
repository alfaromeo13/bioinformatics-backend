import os
import sys


class Errors:
    def __init__(self, protein_chains_list, ligand_chains_list):                 

        self.protein_chains_list = protein_chains_list
        self.ligand_chains_list = ligand_chains_list
        self.all_chains = []
        
    def generate_all_chains(self):
        for chain in self.protein_chains_list:
            self.all_chains.append('PRO' + chain)
        for chain in self.ligand_chains_list:
            self.all_chains.append('PRO' + chain)

    # za v error checking mogoče
    def check_inte_crd(self):
        all_chains_low = [chain.lower() for chain in self.all_chains]
        for chain in all_chains_low:
            crd_file_path = f"{chain}.crd"
            if not os.path.exists(crd_file_path):
                print(f"Error: Coordinate file '{crd_file_path}' not found. Please check input PDB if you specified the correct chains for interaction.")
                sys.exit(1)