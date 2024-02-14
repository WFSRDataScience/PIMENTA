"""


usage: get_taxids_mzg.py [-h] -i INPUT -p DB -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT    the fasta retrieved from MZGdb
  -p DB, --dump DB folder including nodes.dmp and names.dmp retrieved from NCBI, script will automatically download these files if they do not exist.
  -o OUTPUT, --output OUTPUT output file containing 

Outputs:
  <OUTPUT>  taxid map, there are some sequences still without taxid, these have to be manually added
  <INPUT>.short.fasta  This file can be used to create a BLASTdb for PIMENTA. Place the taxonomy BLAST database (available at NCBI) in the same folder before creating BLAST database. 



"""
import sys
import os
import taxopy
import pandas as pd
import argparse

#import warnings
#warnings.filterwarnings("ignore")

def get_header(filename,taxdb,output):
	"""Translate the header from inputfile to BLAST format"""
	f = open(output, 'w')
	with open(filename, 'r') as fd:
		for line in fd.readlines():
			if '>' in line and "BOLD:" in line:
				header=line.strip(">")
				key = header.split("\t")[0]
				seqid = key.split("_")[0]
				species = ' '.join(key.split("_")[-2:])	
				#Translates species names to names that are included in names/nodes.dmp
				if species.startswith('sp'):
#					print(header)
					speciestmp=key.split("_")[-3:]
					species=speciestmp[0]+" "+speciestmp[1]+". "+speciestmp[2]
				elif " sp" in species and species.endswith("sp"):
					species=species+"."
				elif len(key.split("_")[-2]) < 2:
				#	print(header)
					species=""
				elif species == "Acanthacartia tonsa":
					species= "Acartia tonsa"
				elif species == "Acanthacartia californiensis":
					species= "Acartia californiensis"
				elif species == "Acanthacartia clausi":
					species= "Acartia clausii"
				elif "coryphaenae" in species and "Caligus" in species:
					species="Caligus elongatus"
				elif "coryphaenae" in species:
					print(species)
				elif species =="Anomalocera patersonii":
					species="Anomalocera patersoni"
				elif species == "Acanthacartia longiremis":
					species= "Acartia longiremis"
				elif species == "Acartiura longiremis":
					species= "Acartia longiremis"
				elif species == "Acartiura clausi":
					species= "Acartia clausii"
				if species != "":
					try:
						taxid=get_taxid(taxdb,species)
						#print(key+"\t"+str(taxid))
						f.write(seqid+"\t"+str(taxid)+"\n")
					except IndexError:
						speciestmp=key.split("_")[-4:]
						species=speciestmp[0]+" "+speciestmp[1]+". "+speciestmp[2]+" "+speciestmp[3]
						try:
							taxid=get_taxid(taxdb,species)
							f.write(seqid+"\t"+str(taxid)+"\n")
						except IndexError:
							speciestmp=key.split("_")[-3:]
							species=speciestmp[0]+" "+speciestmp[1]+"."
							try:
								taxid=get_taxid(taxdb,species)
								f.write(seqid+"\t"+str(taxid)+"\n")
							except IndexError: 			
								try:
									species=speciestmp[0]+" "+speciestmp[1]+" "+speciestmp[2]			
									taxid=get_taxid(taxdb,species)
									f.write(seqid+"\t"+str(taxid)+"\n")
								except IndexError:
									try:
										species=speciestmp[0]+" "+speciestmp[2]
										taxid=get_taxid(taxdb,species)
										f.write(seqid+"\t"+str(taxid)+"\n")
									except IndexError:
										print(key, species)	
										f.write(seqid+"\n")
#				print(species)
	f.close()


def get_bold(filename,taxdb,output):
	"""WIP""" 
	f = open(output, 'w')
#	atlanta=["Atlanta gaudichaudi","Atlanta inclinata","Atlanta peronii", "Atlanta selvagensis", "Atlanta tokiokai"]

	with open(filename, 'r') as fd:
		for line in fd.readlines():
			key, species =line.rstrip().split("\t")
			speciestmp=species.split(" ")
			if "Neocallichirus guassutinga" in species:
				species="Sergio guassutinga"
			elif "Omalacantha bicornuta" in species:
				species="Microphrys bicornutus"
			#elif "Atlanta" in species and species not in atlanta:
		#		species="Atlanta"
			elif species == "Acanthacartia clausi":
				species= "Acartia clausii"
			if species != "":
				try:
					taxid=get_taxid(taxdb,species)
					#print(key+"\t"+str(taxid))
					f.write(key+"\t"+str(taxid)+"\n")
				except IndexError:
					try:
						species=speciestmp[0]+" "+speciestmp[2]
						taxid=get_taxid(taxdb,species)
						#print(key+"\t"+str(taxid)+"\n")
						f.write(key+"\t"+str(taxid)+"\n")
					except IndexError:
						try:
							taxid=get_taxid(taxdb,speciestmp[0])
							f.write(key+"\t"+str(taxid)+"\n")
						except IndexError:
							print(key+"\t"+species)
							#f.write(key+"\n")
	f.close()

def get_taxid(taxdb,common):
	"""Get taxid from given species/genus/family etc."""
	taxid=taxopy.taxid_from_name(common, taxdb)
	tax=taxopy.Taxon(taxid[0], taxdb)
	level=tax.rank
	return taxid[0]


def rewrite_headers(InFile):
	""" Rewrites input fasta so it can be used with taxonomy db from NCBI and PIMENTA, this still needs to be transformed into a BLAST db"""
	OutFile=InFile+".short.fasta"
	with open(OutFile, 'w') as outfile:
		with open(InFile, 'r') as infile:
			for line in infile:
				if ">NC" in line:
					line = line.split("_")[0]+ "_" + line.split("_")[1]  + " ".join(line.split("_")[2:]) + '\n'
				elif line.startswith('>'):
					line = line.split("_")[0]  + " ".join(line.split("_")[1:]) + '\n'
				print(line, end='', file=outfile)
def main(args):
	taxdb = taxopy.TaxDb(nodes_dmp=args.db+"/nodes.dmp", names_dmp=args.db+"/names.dmp", keep_files=True)
#	taxid=get_taxid(taxdb,"Philine")
	#print(taxid)
#	get_bold(args.input,taxdb,args.output)
	rewrite_headers(args.input)
	get_header(args.input,taxdb,args.output)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", dest="input", required=True)
	parser.add_argument("-p", "--dump", dest="db", required=True)
	parser.add_argument("-o", "--output", dest="output", required=True)       
	args = parser.parse_args()
	main(args)

