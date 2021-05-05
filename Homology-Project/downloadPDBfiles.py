import Bio
from Bio.PDB import PDBList
'''Selecting structures from PDB'''

pdbList_file = open(pdbList_filename)
pdbl = PDBList()


def downloadPDBs(pdbList_filename):

	pdbList=['4B97','4IPH','5W9F']
	filenames = []
	for i in pdbList:
		filename= pdbl.retrieve_pdb_file(i,pdir='PDB', file_format='pdb')
		filenames.append(filename)
	return filenames