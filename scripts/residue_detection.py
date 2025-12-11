from Bio import PDB

class ResidueDetector:
    def __init__(self, pdb_name, cutoff_distance, protein_chains_list, ligand_chains_list, setup_file):
        self.pdb_file = pdb_name
        self.cutoff_distance = cutoff_distance
        self.protein_chains = protein_chains_list
        self.ligand_chains = ligand_chains_list
        self.setup_file = setup_file
        self.detected_residues = []  

    #Izračun razdalje med izbranimi atomi (Evklidska razdalja)
    def calculate_distance(self, coord1, coord2): 
        return (sum((c1 - c2) ** 2 for c1, c2 in zip(coord1, coord2))) ** 0.5

    def detect_distances(self):
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('pdb', self.pdb_file)
        for model in structure:
            protein_residues = {residue for chain_id in self.protein_chains for residue in model[chain_id]}
            ligand_residues = {residue for chain_id in self.ligand_chains for residue in model[chain_id]}

            for protein_residue in protein_residues:
                for ligand_residue in ligand_residues:
                    if protein_residue != ligand_residue:
                        #Izračuna distanco
                        min_distance = float('inf')
                        for protein_atom in protein_residue:
                            for ligand_atom in ligand_residue:
                                distance = self.calculate_distance(protein_atom.get_coord(), ligand_atom.get_coord())
                                min_distance = min(min_distance, distance)

                        if min_distance <= self.cutoff_distance:
                            res_name_protein = protein_residue.get_resname().strip()
                            chain_id_protein = protein_residue.get_parent().get_id()
                            res_num_protein = protein_residue.get_id()[1]

                            res_name_ligand = ligand_residue.get_resname().strip()
                            chain_id_ligand = ligand_residue.get_parent().get_id()
                            res_num_ligand = ligand_residue.get_id()[1]

                            detected_residue_str = f"{res_name_protein}{chain_id_protein}{res_num_protein} {res_name_ligand}{chain_id_ligand}{res_num_ligand}"
                            self.detected_residues.append(detected_residue_str)

    def print_detected_residues(self):
        for residue_str in self.detected_residues:
            print("Detected residues in contact:")
            print(residue_str)

    def append_detected_residues_to_setup(self):
        with open(self.setup_file, 'r') as file:
            setup_content = file.read()

        #Posodobitev liste mutacij
        mutation_list_line = None
        lines = setup_content.splitlines()
        for idx, line in enumerate(lines):
            if line.startswith('MUTATION_LIST'):
                mutation_list_line = idx
                break

        if mutation_list_line is not None:
            existing_mutations = set(lines[mutation_list_line].split()[1:])
            new_mutations = [mutation for mutation in self.detected_residues if mutation not in existing_mutations]
            if new_mutations:
                lines[mutation_list_line] = 'MUTATION_LIST ' + ' '.join(existing_mutations.union(new_mutations))
                with open(self.setup_file, 'w') as file:
                    file.write('\n'.join(lines))

    #Funkcija preveri mutacije in izbriše duplikate
    def check_mutations(self):
        with open(self.setup_file, 'r') as file:
            lines = file.readlines()

        #Ekstrakcija iz mutation_lista
        mutations = []
        for line in lines:
            if line.startswith('MUTATION_LIST'):
                mutation_entries = line.split()[1:]
                for mutation in mutation_entries:
                    mutations.append(mutation)

        #Odstranitev duplikatov
        unique_mutations = []
        seen_mutations = set()
        for mutation in mutations:
            if mutation not in seen_mutations:
                unique_mutations.append(mutation)
                seen_mutations.add(mutation)

        #Filtriranje če slučajno kakšna AK ni korektna
        valid_amino_acids = set([
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
        ])
        valid_mutations = [mutation for mutation in unique_mutations if mutation[:3] in valid_amino_acids]
        filtered_valid_mutations = []
        for mutation in valid_mutations:
            protein_chain = mutation[3]  
            if protein_chain in self.protein_chains:
                filtered_valid_mutations.append(mutation)

        #Sortiranje mutacij po chainu in zaporedni številki AK
        def sort_key(mutation):
            return (mutation[3], int(mutation[4:].split('.')[0]))

        sorted_mutations = sorted(filtered_valid_mutations, key=sort_key)

        #Generacija nove liste mutacij
        new_mutation_list_line = 'MUTATION_LIST ' + ' '.join(sorted_mutations) + '\n'

        #Menjava stare liste mutacij z novo urejeno listo
        new_lines = []
        for line in lines:
            if line.startswith('MUTATION_LIST'):
                new_lines.append(new_mutation_list_line)
            else:
                new_lines.append(line)
        with open(self.setup_file, 'w') as file:
            file.writelines(new_lines)





