import Bio.PDB
from Bio import AlignIO
import itertools
import argparse
import csv
from tqdm import tqdm



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
			for idx1 in range(difference1 - 1): # -1 for indexing
				matrix.append([(0, previous_position1 + idx1 + 1, pos) for pos in range(matrix[0][0][2], matrix[0][-1][2] + 1)])

		contacts = []

		for residue2 in principle_chain:
			if residue2.get_id()[0] != ' ' and residue2.get_resname() not in amino_acids.keys():
				continue

			difference2 = residue2.get_id()[1] - previous_position2
			if difference2 > 1:
				for idx2 in range(difference2 - 1): # - 1 for indexing
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
	parser = Bio.PDB.PDBParser(QUIET=True)
	structure = parser.get_structure(pdbfile.split(".")[0], pdbfile)

	for residue in structure[0][chain]:
		if residue.get_id()[0] == ' ':
			if residue.get_resname() not in amino_acids.keys():
				return False

	return True



def merge_alignment(pdbfile, alignment):
	parser = Bio.PDB.PDBParser(QUIET=True)
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
		



def create_proteins(combination, num_proteins, contact_matrix, alignment, batch_size, myjob_id): # this is the function that makes the dictionaries
	'''
	Protein representation as dictionary: {(crossover_position1, crossover_position2): parent in that range}
	
	Input: single tuple from create_combinations()
	Output: returns list of dictionaries of proteins from one crossover combination
	'''


	crossover_ranges = [(0, combination[0])]

	for idx in range(1, len(combination)):
		crossover = (combination[idx - 1] + 1, combination[idx])
		crossover_ranges.append(crossover)

	crossover_ranges.append((combination[-1] + 1, len(contact_matrix) - 1))
	protein_identities = list(itertools.product(range(num_proteins), repeat=len(combination) + 1)) # this generates all combinations

	proteins = []
	sequences = []
	disruptions = []
	distances = []

	start_idx = myjob_id*batch_size

	if (start_idx + batch_size) >= (len(protein_identities) - 1):
		stop_idx = len(protein_identities)
	else:
		stop_idx = start_idx + batch_size
	
	print('START IDX: {}'.format(start_idx))
	print('STOP IDX: {}'.format(stop_idx))


	for i in tqdm(range(start_idx, stop_idx)): # this goes through each of the tuples generated above to make the dictionaries
		
		protein_dict = {}

		for idx in range(len(protein_identities[i])):
			protein_dict[crossover_ranges[idx]] = protein_identities[i][idx] # makes each individual dictionary

		sequence = get_sequence(alignment, protein_dict)

		if sequence not in sequences:
			disruption = calculate_disruptions(contact_matrix, alignment, protein_dict)
			distance = calculate_distance(protein_dict, alignment)

			sequences.append(sequence)
			distances.append(distance)
			disruptions.append(disruption)
			proteins.append(protein_dict)

	return proteins, disruptions, distances



def A(alignment, parent, k):
	return alignment[parent][k]


def probability(alignment, i, j, s_i, s_j):
	for parent in range(len(alignment)):
		if A(alignment, s_i, i) == A(alignment, parent, i) and A(alignment, s_j, j) == A(alignment, parent, j):
			if parent == s_i or parent == s_j: # not sure about this
				return 0
	return 1


def get_parent(protein, residue):
	'''
	Returns parent for residue position of a recombined protein
	{(crossover_position1, crossover_position2): parent in that range}
	'''

	for key in protein.keys():
		if residue >= key[0] and residue <= key[1]:
			return protein[key]



def calculate_disruptions(contact_matrix, alignment, protein):
	'''
	Calculates number of contacts broken for a protein chimera
	'''

	contacts_broken = 0
	residue_pairs = list(itertools.combinations(range(len(contact_matrix)), 2))

	for pair in residue_pairs:
		i = pair[0]
		j = pair[1]

		s_i = get_parent(protein, i)
		s_j = get_parent(protein, j)
		contact_value = contact_matrix[i][j][0]*probability(alignment, i, j, s_i, s_j)
		contacts_broken += contact_value

	return contacts_broken



def get_sequence(alignment, protein):
	sequence = ''
	parents = list(range(len(alignment)))

	for i in range(len(alignment[0])):

		s_i = get_parent(protein, i)
		residue = A(alignment, s_i, i)

		if residue != "-":
			sequence += residue
		else:

			for parent in parents:
				if parent != s_i:
					residue = A(alignment, parent, i)

					if residue != '-':
						sequence += residue
						break

				if parent == len(parents) - 1:
					residue = A(alignment, s_i, i)
					sequence += residue


	return sequence




def get_real_crossover_points(points, contact_matrix):
	real_points = []

	for point in points:
		real_points.append(contact_matrix[point][0][1])

	return real_points




def calculate_distance(protein, alignment):
	distance_counts = []
	
	for parent in alignment:
		count = 0

		regions = list(protein.keys())
		regions.sort()

		for idx in range(regions[-1][1]):
			residue_chimera = A(alignment, get_parent(protein, idx), idx)
			residue_parent = parent[idx]

			if residue_chimera != residue_parent:
				count += 1

		distance_counts.append(count)

	return min(distance_counts)





def get_jobid(txt_filename, jobid):

	d = open('{}.txt'.format(txt_filename,'r'))
	t = d.read()
	d.close()

	jobs_temp = [i for i in t.split('\n')]

	if jobs_temp[-1] == '' or jobs_temp[-1] == ' ':
		jobs_temp.pop(-1)

	jobs = [int(i) for i in jobs_temp]
	myjob_id = jobs[jobid]

	print("ARRAY TASK ID: " + str(jobid))
	print("MYJOB ID: " + str(myjob_id))

	return myjob_id




def read_crossover_points():
	f = open("crossover_points.txt", "r")
	t = f.read()
	f.close()
	
	temp_points = t.split('\n')[0]
	points = [int(i) for i in temp_points.split(',')]

	return points





def main(num_proteins, pdbfile, alignment_file, specified_points, chain, batch_size, jobid):
	if chain != "all":
		contact_matrix = create_matrix_intra(pdbfile, chain)
		temp_alignment = AlignIO.read(alignment_file, "clustal")
		alignment = modify_alignment(pdbfile, chain, temp_alignment)

	if not specified_points:
		points = input("Enter crossover points:\n")
		crossover_positions = [int(num) for num in points.split(" ")]
	else:
		crossover_positions = read_crossover_points()

	myjob_id = get_jobid('jobs', jobid)
	proteins, disruptions, distances = create_proteins(crossover_positions, num_proteins, contact_matrix, alignment, batch_size, myjob_id)

	return proteins, disruptions, distances, contact_matrix, myjob_id




def get_crossover_position_identity(protein, contact_matrix):
	protein_lst = [key[1] for key in protein.keys()]
	end_point = max(protein_lst)
	protein_lst.remove(end_point)
	protein_lst.sort()

	actual_points = get_real_crossover_points(protein_lst, contact_matrix)

	for position in protein_lst.copy():
		actual_points.append(get_parent(protein, position))

	actual_points.append(get_parent(protein, end_point))

	return actual_points





def to_csv(proteins, disruptions, distances, contact_matrix, myjob_id, name):
	fieldnames_crossovers = []
	fieldnames_regions = []

	num_regions = len(proteins[0].keys())
	num_crossovers = num_regions - 1

	fieldnames_crossovers = ['crossover_{}'.format(str(idx + 1)) for idx in range(num_crossovers)]
	fieldnames_regions = ['region_{}'.format(str(idx + 1)) for idx in range(num_regions)]

	fieldnames = fieldnames_crossovers + fieldnames_regions
	fieldnames.extend(['contacts_broken', 'distance'])


	with open(('{}_{}.csv'.format(name, myjob_id)), 'w') as csvfile:

		writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
		writer.writeheader()


		for idx in range(len(proteins)):

			format_protein = {}
			values = get_crossover_position_identity(proteins[idx], contact_matrix)

			for value_idx in range(len(values)):
				format_protein[fieldnames[value_idx]] = values[value_idx]

			format_protein['contacts_broken'] = disruptions[idx]
			format_protein['distance'] = distances[idx]

			writer.writerow(format_protein)




if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument('--pdb', action='store', required=True, help='PDB file name')
	parser.add_argument('--align', action='store', required=True, help='Alignment file name')
	parser.add_argument('-p', '--proteins', type=int, required=True, help='Number of proteins')
	parser.add_argument('-nr', '--notrandom', action='store_false', default=True, help='Randomize crossover points')
	parser.add_argument('--chain', action='store', default='A', help='Chain to make chimeras')
	parser.add_argument('-sc', '--savecsv', action='store_true', default=False, help='Save results to csv')
	parser.add_argument('-n', '--name', action='store', required=True, help='Name for saved files')
	parser.add_argument('-b', '--batch_size', type=int, required=True, help='Number of entries per output file',default=40000)
	parser.add_argument('-j', '--jobid', type=int, required=True, help='Job array id')
	args = parser.parse_args()


	pdb_file = "../PDB_Files/{}.pdb".format(args.pdb)
	alignment_file = "../Alignment_Files/{}.clustal_num".format(args.align)

	proteins, contacts_broken, distances, matrix, myjob_id = main(args.proteins, pdb_file, alignment_file, args.notrandom, args.chain, args.batch_size, args.jobid)

	if args.savecsv:
		to_csv(proteins, contacts_broken, distances, matrix, myjob_id, args.name)
		print("FINISHED BATCH %d SUCCESSFULLY" %myjob_id)



