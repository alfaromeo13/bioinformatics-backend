import os
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap


class HeatmapGenerator:
    def __init__(self, mut_chid, resn, wtres, method):
        self.method = method
        self.mut_chid = mut_chid.lower()
        self.mut_chid_caps = mut_chid
        self.wtres = wtres.lower()
        self.resn = resn
        self.lines_with_pro = []
        self.ener_dictionary = {}

    def inter_ener_data(self):
        file_pattern = "{}_ener_{}{}{}.dat".format(self.method,self.wtres, self.resn, self.mut_chid)
        if os.path.exists(file_pattern):
            with open(file_pattern, 'r') as file:
                for line in file:
                    if line.startswith("PRO"):
                        self.lines_with_pro.append(line.strip())

    def extract_info(self, line):
        chain_id = line[3]
        resn_dat = re.findall(r'_(\d+)_', line)[0]
        wt = re.findall(r'_(\d+)_([A-Za-z])', line)[0][1]
        mut = re.findall(r'\d([A-Za-z])', line)[0]
        energy = round(float(line.split()[1]), 2)
        if wt == mut:
            ref_ener = energy
        else:
            ref_ener = None
        return chain_id, int(resn_dat), wt, mut, energy, ref_ener

    def get_ref_ener_lines(self):
        ref_ener_lines = {}
        for line in self.lines_with_pro:
            chain_id, resn_dat, wt, mut, energy, ref_ener = self.extract_info(line)
            if ref_ener is not None:
                ref_ener_lines[(chain_id, resn_dat)] = ref_ener
        return ref_ener_lines

    def create_ener_diff_dict(self):
        ref_ener_lines = self.get_ref_ener_lines()
        ener_diff_dict = {}

        for line in self.lines_with_pro:
            chain_id, resn_dat, wt, mut, energy, _ = self.extract_info(line)
            if chain_id == self.mut_chid_caps and resn_dat == self.resn:
                ref_ener = None
                for (ref_chain_id, ref_resn_dat), ref_ener_val in ref_ener_lines.items():
                    if chain_id == ref_chain_id and resn_dat == ref_resn_dat:
                        ref_ener = ref_ener_val
                        break

                if ref_ener is not None:
                    ener_diff = energy - ref_ener
                    ener_diff_dict[mut] = ener_diff

        max_abs_ener_diff = max(abs(val) for val in ener_diff_dict.values())
        max_ener_diff = max_abs_ener_diff if max_abs_ener_diff > 0 else 1
        min_ener_diff = -max_ener_diff

        normalized_ener_diff_dict = {}
        for mut, ener_diff in ener_diff_dict.items():
            normalized_ener_diff = (2 * (ener_diff - min_ener_diff) / (max_ener_diff - min_ener_diff)) - 1
            normalized_ener_diff_dict[mut] = normalized_ener_diff

        with open("{}_ener_diff.dat".format(self.method), "a") as file:
            for mut, ener_diff in normalized_ener_diff_dict.items():
                line = f"Mutation:{mut}{self.mut_chid_caps}{self.resn} Energy Difference: {ener_diff}\n"
                file.write(line)

        return ener_diff_dict, normalized_ener_diff_dict
    
    

    def make_heatmap(self, matrix_norm, matrix_abs, mutation_list, amino_acids, pdb_name, method):
        amino_acid_mapping = {
            'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'E': 'GLU', 'Q': 'GLN',
            'G': 'GLY', 'H': 'HSD', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE',
            'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
        }

        matrix_norm = np.array(matrix_norm)
        matrix_abs = np.array(matrix_abs)
        num_mutations = len(mutation_list)
        fig_height = max(min(0.6 * num_mutations, 4), 1.8)
        font_size = max(min(12 - num_mutations, 9), 7)

        plt.figure(figsize=(12, fig_height))
        plt.rcParams.update({'font.size': font_size})

        # Find any value that is 2 SD above or below the mean in the absolute matrix and set it to dark red
        mean_abs = np.mean(matrix_abs)
        std_abs = np.std(matrix_abs)
        outliers_above = matrix_abs > mean_abs + 5 * std_abs
        outliers_below = matrix_abs < mean_abs - 5 * std_abs

        # Combine outliers in both directions
        outliers = outliers_above | outliers_below

        # Create a colormap with dark red for outliers on the absolute heatmap
        cmap_abs = plt.cm.viridis
        custom_cmap_abs = ListedColormap([cmap_abs(i) for i in range(cmap_abs.N)])
        custom_cmap_abs.set_bad('firebrick')  # Change to dark red

        # Create a copy of the absolute matrix to preserve original values
        matrix_abs_colored = matrix_abs.copy()

        # Replace the outlier values with NaN for plotting purposes on the absolute heatmap
        matrix_abs_colored[outliers] = np.nan

        # Create a colormap for the normalized heatmap
        cmap_norm = plt.cm.viridis

        plt.imshow(matrix_norm, cmap=cmap_norm, aspect='auto', vmin=np.nanmin(matrix_norm), vmax=np.nanmax(matrix_norm))
        for i in range(len(amino_acids)):
            for j in range(len(mutation_list)):
                value = matrix_norm[j][i]
                color = 'white' if value <= 0.7 else 'black'  # Adjust the threshold here (0.7 for black text)
                plt.text(i, j, f'{value:.2f}', ha='center', va='center', color=color, fontsize=min(font_size, 8))

        colorbar = plt.colorbar()
        colorbar.set_label('Interaction Energy')
        plt.yticks(range(len(mutation_list)), mutation_list, fontsize=min(font_size, 8))
        plt.ylabel('Mutations', fontsize=min(font_size, 8))
        plt.xticks(range(len(amino_acids)), [amino_acid_mapping[aa] for aa in amino_acids])  # Use three-letter codes
        plt.xlabel('Amino Acid')
        plt.title('Normalized Energy Differences')
        plt.tight_layout()

        pdb_name_short = pdb_name[0:-4]
        plt.savefig('{}_norm_heatmap_{}.jpg'.format(method, pdb_name_short), format='jpg', dpi=600)
        plt.savefig('{}_norm_heatmap_{}.svg'.format(method, pdb_name_short), format='svg')
        print("Generated the {} heatmap norm_heatmap_{}.svg".format(method, pdb_name_short))
        plt.clf()

        plt.figure(figsize=(12, fig_height))
        plt.rcParams.update({'font.size': font_size})

        # Display the colored absolute matrix with outlier values
        plt.imshow(matrix_abs_colored, cmap=custom_cmap_abs, aspect='auto', vmin=np.nanmin(matrix_abs_colored), vmax=np.nanmax(matrix_abs_colored))

        threshold = np.percentile(matrix_abs, 80)
        for i in range(len(amino_acids)):
            for j in range(len(mutation_list)):
                value = matrix_abs[j][i]
                color = 'black' if value >= threshold else 'white'
                plt.text(i, j, f'{value:.2f}', ha='center', va='center', color=color, fontsize=min(font_size, 8))

        colorbar = plt.colorbar(format='%.0f')
        colorbar.set_label('Interaction Energy')
        plt.yticks(range(len(mutation_list)), mutation_list, fontsize=min(font_size, 8))
        plt.ylabel('Mutations', fontsize=min(font_size, 8))
        plt.xticks(range(len(amino_acids)), [amino_acid_mapping[aa] for aa in amino_acids])  # Use three-letter codes
        plt.xlabel('Amino Acid')
        plt.title('Absolute Energy Differences')
        plt.tight_layout()

        plt.savefig('{}_abs_heatmap_{}.jpg'.format(method, pdb_name_short), format='jpg', dpi=600)
        plt.savefig('{}_abs_heatmap_{}.svg'.format(method, pdb_name_short), format='svg')
        print("Generated the {} heatmap abs_heatmap_{}.svg\n".format(method, pdb_name_short))
