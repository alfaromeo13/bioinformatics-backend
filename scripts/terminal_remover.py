class TerminalRemover:
    def __init__(self, pdb_name, setup_file):
        self.pdb_name = pdb_name
        self.setup_file = setup_file
        self.detected_terminals = {}

    #Najdemo prvo in zadnjo amino-kislino    
    def find_terminals(self):
        chain_terminals = {}
        
        with open(self.pdb_name, 'r') as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    chain_id = line[21]
                    residue_name = line[17:20]
                    residue_number = line[22:26].strip()
                    residue_identifier = f"{residue_name}{chain_id}{residue_number}"
                    
                    if chain_id not in chain_terminals:
                        chain_terminals[chain_id] = {"first": None, "last": None}
                    
                    if chain_terminals[chain_id]["first"] is None:
                        chain_terminals[chain_id]["first"] = residue_identifier
                    chain_terminals[chain_id]["last"] = residue_identifier
                    
        self.detected_terminals = chain_terminals

    #Odstranimo prvo in zadnjo amino-kislino iz setup datoteke    
    def remove_terminals_from_setup(self):
        if not self.detected_terminals:
            return
            
        with open(self.setup_file, 'r+') as setup:
            lines = setup.readlines()
            setup.seek(0)
            for line in lines:
                if line.startswith("MUTATION_LIST"):
                    mutations = line.split()[1:]
                    filtered_mutations = []
                    for mutation in mutations:
                        keep_mutation = True
                        for chain_id, terminals in self.detected_terminals.items():
                            if mutation.lower() == terminals["first"].lower() or mutation.lower() == terminals["last"].lower():
                                keep_mutation = False
                                break
                        if keep_mutation:
                            filtered_mutations.append(mutation)
                    
                    if filtered_mutations:
                        setup.write("MUTATION_LIST " + " ".join(filtered_mutations) + '\n')
                    else:
                        setup.write(line)  
                else:
                    setup.write(line)
            setup.truncate()
