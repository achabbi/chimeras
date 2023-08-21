import Bio.PDB
from Bio import AlignIO
import random
import argparse



amino_acids = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     		'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     		'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     		'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}



def create_matrix_intra(pdbfile, main_chain):
	'''
	Creates contact matrix represented as a list of lists for specific chain of PDB (intramolecular contacts only)

	Input: --pdbfile (string): path for PDB file
			--main_chain (string): chain of PDB file

	Output: list of lists, where each row/column corresponds to a residue pair. 
			Tuple entries in matrix have the following format:
			(1 if residues contacting or 0 if residues not contacting, 
			residue 1 (row) position actual, residue 2 (column) position actual)
	'''

	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbfile.split(".")[0] + "_chain" + main_chain, pdbfile)

	matrix = []
	principle_chain = structure[0][main_chain]


	for residue in principle_chain:
		if residue.get_id()[0] == ' ' and residue.get_resname() in amino_acids.keys():
			previous_position1 = residue.get_id()[1]
			previous_position2 = residue.get_id()[1]
			break


	for residue1 in principle_chain:
		if residue1.get_id()[0] != ' ' and residue1.get_resname() not in amino_acids.keys():
			continue

		difference1 = residue1.get_id()[1] - previous_position1
		if difference1 > 1:
			for idx1 in range(difference1 - 1): # -1 for indexing purposes
				matrix.append([(0, previous_position1 + idx1 + 1, pos) for pos in range(matrix[0][0][2], matrix[0][-1][2] + 1)])

		contacts = []

		for residue2 in principle_chain:
			if residue2.get_id()[0] != ' ' and residue2.get_resname() not in amino_acids.keys():
				continue

			difference2 = residue2.get_id()[1] - previous_position2
			if difference2 > 1:
				for idx2 in range(difference2 - 1): # - 1 for indexing purposes
					contacts.append((0, residue1.get_id()[1], previous_position2 + idx2 + 1))

			contacting = False

			for atom1 in residue1:
				for atom2 in residue2:

					if atom1 - atom2 <= 4.5:
						contacts.append((1, residue1.get_id()[1], residue2.get_id()[1]))
						contacting = True
						break

				if contacting:
					break

			if not contacting:
				contacts.append((0, residue1.get_id()[1], residue2.get_id()[1]))
				
			previous_position2 = residue2.get_id()[1]

		matrix.append(contacts)
		previous_position1 = residue1.get_id()[1]

	return matrix






def modify_alignment(pdbfile, main_chain, alignment):
	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbfile.split(".")[0] + "_chain", pdbfile)

	principle_chain = structure[0][main_chain]

	sequence = ''
	for residue in principle_chain:
		if residue.get_id()[0] == ' ' and residue.get_resname() in amino_acids.keys():
			sequence += amino_acids[residue.get_resname()]


	temp_alignment = []

	for seq in alignment:
		temp = ''

		for residue in seq:
			temp += str(residue)

		temp_alignment.append(temp)

	missed_residues_front = temp_alignment[0].find(sequence[:3])

	new_alignment = []

	for seq in temp_alignment:
		new_alignment.append(seq[missed_residues_front:])

	missed_residues_back = new_alignment[0].find(sequence[-3:])
	missed_addon = missed_residues_back + 3

	new_alignment_2 = []

	for idx in range(len(new_alignment)):
		new_alignment_2.append(new_alignment[idx][:missed_addon])
		
	return new_alignment_2



def check_chain(pdbfile, chain):
	parser = Bio.PDB.PDBParser()
	structure = parser.get_structure(pdbfile.split(".")[0], pdbfile)

	for residue in structure[0][chain]:
		if residue.get_id()[0] == ' ':
			if residue.get_resname() not in amino_acids.keys():
				return False

	return True



def merge_alignment(pdbfile, alignment):
	parser = Bio.PDB.PDBParser()
	structure = parser.get_structure(pdbfile.split(".")[0], pdbfile)
	alignment_list = []
	final_alignment = []

	for chain in structure[0]:
		chain_str = str(chain).split('=')[1][0]

		if check_chain(pdbfile, chain_str):
			alignment_list.append(modify_alignment(pdbfile, chain_str, alignment))

	for alignment_idx in range(len(alignment)):
		alignment = ''
		for chain in alignment_list:
			alignment += chain[alignment_idx]

		final_alignment.append(alignment)

	return final_alignment
		




def create_combination(num_crossovers, alignment, contact_matrix):

	start_positions = []

	for sequence in alignment:
		for idx in range(len(sequence)):
			if sequence[idx] != '-':
				start_positions.append(idx)


	start = min([position for position in start_positions if start_positions.count(position) > 1])
	end_positions = []

	for sequence in alignment:
		for idx in range(start, len(sequence)):
			if sequence[idx] != '-':
				end_positions.append(idx)

	end = max([position for position in end_positions if end_positions.count(position) > 1])
	print(start, end)

	combination_unsorted = random.sample(range(start, end), num_crossovers)
	combination_unsorted.sort()

	actual_crossover_points = []
	for num in combination_unsorted:
		actual_crossover_points.append(contact_matrix[0][num][2])

	print("Actual Crossover Points: " + str(actual_crossover_points))
	print("Modified Crossover Points: " + str(combination_unsorted))

	return combination_unsorted, actual_crossover_points



def new_algo(contact_matrix, alignment):

	modified_matrix = []

	for row in contact_matrix:
		new_row = []
		for tple in row:
			new_row.append(tple[0])
		modified_matrix.append(new_row)

	contact_array = [row.count(1) for row in modified_matrix]
	sum_array = []

	for idx in range(len(contact_array)):
		residues = []

		for parent in alignment:
			residues.append(parent[idx])

		val = contact_array[idx]*(len(set(residues))-1)
		sum_array.append(val)

	return sum_array





def main(num_crossovers, pdbfile, alignment_file, specified_points, chain):

	if chain != "all":
		contact_matrix = create_matrix_intra(pdbfile, chain)
		temp_alignment = AlignIO.read(alignment_file, "clustal")
		alignment = modify_alignment(pdbfile, chain, temp_alignment)
	else:
		contact_matrix = create_matrix_inter(pdbfile)
		temp_alignment = AlignIO.read(alignment_file, "clustal")
		alignment = merge_alignment(pdbfile, temp_alignment)

	array = new_algo(contact_matrix, alignment)
	print(array)




def save_txt(crossover_points, actual_crossover_points):
	f = open("crossover_points.txt", "a")

	for idx in range(len(crossover_points) - 1):
		f.write(str(crossover_points[idx]))
		f.write(',')

	f.write(str(crossover_points[-1]))
	f.write('\n')

	for idx in range(len(actual_crossover_points) - 1):
		f.write(str(actual_crossover_points[idx]))
		f.write(',')

	f.write(str(actual_crossover_points[-1]))
	f.close()




if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('--pdb', action='store', required=True, help='PDB file name')
	parser.add_argument('--align', action='store', required=True, help='Alignment file name')
	parser.add_argument('-c', '--crossovers', type=int, required=True, help='Number of crossovers')
	parser.add_argument('-nr', '--notrandom', action='store_false', default=True, help='Randomize crossover points')
	parser.add_argument('--chain', action='store', default='A', help='Chain to make chimeras')
	args = parser.parse_args()



	pdb_file = "../PDB_Files/{}.pdb".format(args.pdb)
	alignment_file = "../Alignment_Files/{}.clustal_num".format(args.align)

	crossover_points, actual_crossover_points = main(args.crossovers, pdb_file, alignment_file, args.notrandom, args.chain)
	save_txt(crossover_points, actual_crossover_points)

	main(args.crossovers, pdb_file, alignment_file, args.notrandom, args.chain)

