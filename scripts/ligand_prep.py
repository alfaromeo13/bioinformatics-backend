from collections import defaultdict
import sys

class LigandPreparation:
    def __init__(self,cgenff_name):
        self.cgenff_name = cgenff_name

    #Razdeli cgenff output v dve datoteki, ki ju potem kliče CHARMM
    def split_cgen_out(self):
        input_file = self.cgenff_name
        rtf_out = 'ligand.rtf'
        prm_out = 'ligand.prm'

        rtf_lines = []
        prm_lines = []

        with open(input_file, 'r') as file:
            lines = file.readlines()
            start_rtf = False
            start_prm = False

            for line in lines:
                if line.startswith('* Topologies'):
                    start_rtf = True
                    rtf_lines.append(line)
                elif line.startswith('* Parameters'):
                    start_prm = True
                    prm_lines.append(line)

                if start_rtf:
                    rtf_lines.append(line)
                if start_prm:
                    prm_lines.append(line)

                if line.startswith('END') and start_rtf:
                    start_rtf = False
                if line.startswith('END') and start_prm:
                    start_prm = False

        with open(rtf_out, 'w') as rtf_file:
            rtf_file.writelines(rtf_lines)

        with open(prm_out, 'w') as prm_file:
            prm_file.writelines(prm_lines)

    #Vsi vodiki so prerazporejeni pod ustrezne ogljike
    #Avtor funkcije hydrogens_under_carbons Samo Lešnik
    def hydrogens_under_carbons(self):
        bound_to = defaultdict(list)
        atom_line = defaultdict(list)
        beggining = []
        end = []
        element = lambda a: a[0]

        with open('ligand.rtf') as RTF:
            switch = True
            for line in RTF:
                if switch == True and not line.startswith("ATOM"):
                    beggining.append(line)
                spl = line.split()
                if line.startswith("ATOM"):
                    switch = False
                    atom = spl[1].strip()
                    nline = line.strip() + "\n" 
                    atom_line[atom] = nline
                if line.startswith("BOND"):
                    end.append(line)
                    first = spl[1].strip("\n").strip()
                    second = spl[2].strip("\n").strip()
                    if element(first) != "H" and element(second) == "H":
                        bound_to[first].append(second)
                    if element(second) != "H" and element(first) == "H":
                        bound_to[second].append(first)
                if (not line.startswith("ATOM")) and (not line.startswith("BOND")) and switch == False:
                    end.append(line)

        atoms = []

        for atom in atom_line:
            if element(atom) != "H":
                atoms.append(atom_line[atom])
                for h in bound_to[atom]:
                    atoms.append(atom_line[h])
            else:
                if element(atom) != "H":
                    atoms.append(atom_line[atom])

        with open('ligand.rtf', "w") as lig:
            sys.stdout = lig
            for text in beggining:
                print(text, end="")

            for text in atoms:
                print(text, end="")

            for text in end:
                print(text, end="")
            sys.stdout = sys.__stdout__


