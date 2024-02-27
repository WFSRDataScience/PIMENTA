#!/usr/bin/env python3
import sys
import os
import taxopy
import pandas as pd
import argparse

print(sys.path)
#replaces Taxonomy.R, retrieves taxonomy from blast output
#Does this for whole dir instead for each marker/target
#Retrieves cluster size (from cluster name) and seq length from fasta
#
def get_excluded_seqids(file):
	df = pd.read_csv(file, sep=";")
	return df["sacc"].tolist()

def Taxonomy(file, filename, taxdb, excluded_seqids, dir):
	df = pd.read_csv(file, sep="\t")
	df_tax=[]
	marker=filename.split(".")[-1]
	print(excluded_seqids)
	for qseqid, df_qseqid in df.groupby('qseqid'):
		#Expand df sallacc
		df_qseqid.sallacc=df_qseqid.sallacc.str.split(';')
		df_qseqid_expanded=df_qseqid.explode('sallacc')
		#Filtering out seqids
		#print(df_qseqid_expanded)
		filtered_df = df_qseqid_expanded[~df_qseqid_expanded["sallacc"].isin(excluded_seqids)]
#		print(df_qseqid_expanded[df_qseqid_expanded["sallacc"].isin(excluded_seqids)])
		#retrieve top3 bitscores
		tophit=filtered_df.nlargest(1, 'bitscore')
		#Retrieve all results wilt top1 bitscores
		tophit_expanded=filtered_df.loc[filtered_df['bitscore'].isin(tophit["bitscore"].tolist())]
		#Get number of Hits
		hits=len(tophit_expanded["sallacc"].tolist())
		#Expand df taxids
		try:
			tophit_expanded.staxids=tophit_expanded.staxids.str.split(';')
			tophit_expanded=tophit_expanded.explode('staxids')
		except AttributeError:
			print(tophit_expanded["staxids"])

	#	print(tophit_expanded)
		#Get clusterID and size
		try:
			clusterID=qseqid.split(";")[0]
			size=qseqid.split(";")[1].lstrip("size=")
		except AttributeError:
			try:
				size=get_size(qseqid,dir)
			except FileNotFoundError:
				size=1
		#Get result
		level="species"
		taxid="1"
		try:
			lineage,tax,taxid,level=lineage_from_taxopy(tophit_expanded["staxids"].tolist(), taxdb,hits)
		except (IndexError, taxopy.exceptions.TaxidError) as error:
			print(tophit_expanded)
			print(tophit_expanded["staxids"].tolist())
			lineage='NA'
			tax='NA'
			taxid='NA'
			level='NA'
		result=[qseqid,taxid,level,tax,hits,lineage,size,marker]
		print(result)
		df_tax.append(result)
	df = pd.DataFrame(df_tax, columns = ["#qseqid","TAXID","Taxonomic Level","Taxonomic Name","Number of Top Hits","Lineage","Reads", "Marker"])
	df=df.sort_values(by="Lineage")
	return df


def lineage_from_taxopy(taxids, taxdb, hits):
	lineages=[]
	keys=["species", "genus", "family", "order","class", "phylum", "kingdom", "superkingdom"]
	set_taxids = set(taxids)
	taxids = list(set_taxids)
	if len(set_taxids) == 1: #and hits >= 3:
		#print("only one taxid")
		taxid_final=taxids[0]
		tax=taxopy.Taxon(int(taxid_final), taxdb)
		lineage=get_lineage(tax)
		lineage.reverse()
		final_lineage=",".join(lineage)
		taxresult=tax.name
		level=tax.rank
	else:
		for taxid in taxids:
			try:
				tax=taxopy.Taxon(int(taxid), taxdb)
				#Get lineage
				lineage=get_lineage(tax)
				lineages.append(lineage)
				#print(lineage)
				if lineage.count("None") > 2:
					print("WARNING: Possible bad TaxID:")
					print(taxid)
					print(lineage)

			except taxopy.exceptions.TaxidError:
				print("WARNING: The input integer is not a valid NCBI taxonomic identifier. taxid:")
				print(taxid)
		lineages_first=lineages[0]
		if len(lineages) == 1:
			common_links=lineages[0]
		else:
			#get list with only common ranks
			common_links = sorted(set(lineages_first).intersection(*lineages[1:]), key=lambda x:lineages_first.index(x))
		if common_links[0] == 'None':
			#print(common_links)
			common_links=common_links[1:]
#			print("None is last common rank")
		if common_links[0] == "Metazoa" or common_links[0] == "Eukaryote":
			print("WARNING: Cluster identified as Eukaryote or Metazoa")
			print(lineages)	
		common_links.reverse()
#		print(lineages_first)
		difference=len(lineages_first) - len(common_links)
		final_lineage=','.join(common_links)
		taxresult=common_links[-1]
		taxid_final, level=get_rank_and_id_common(taxdb,taxresult)

	return final_lineage, taxresult, taxid_final, level

def get_lineage(tax):
	keys=["species", "genus", "family", "order", "class", "subphylum", "phylum", "kingdom", "superkingdom"]
	dic=tax.rank_name_dictionary
	lineage=[dic.get(key) for key in keys]
	lineage_final=['None' if v is None else v for v in lineage]
	return lineage_final

def get_rank_and_id_common(taxdb,common):
	taxid=taxopy.taxid_from_name(common, taxdb)
	print(taxid)
	tax=taxopy.Taxon(taxid[0], taxdb)
	level=tax.rank
	return taxid[0], level


def create_platform_table(taxdf):
	print(taxdf)

def get_size(ClusterID,dir):
	filename=dir+"/../ClustCons/multi-seq/"+str(ClusterID)+".fasta"
	count = 0
	with open(filename) as f:
		for i in f:
			if i.startswith('>'):
				count += 1
	return count

def main(args):
	excluded_seqids=get_excluded_seqids(args.excludedseqids)
	taxdb = taxopy.TaxDb(nodes_dmp=args.db+"/nodes.dmp", names_dmp=args.db+"/names.dmp", keep_files=True)
	for file in os.listdir(args.blast):
		print(file)
		if file.endswith(".blast.filtered.tsv"):
			filename = file.replace('.blast.filtered.tsv','')
			taxdf=Taxonomy(os.path.join(args.blast, file),filename,taxdb,excluded_seqids,args.blast)
			#save original Taxonomy.R output
			Taxonomy_filename=args.outdir+"/"+filename+".Taxonomy.summary"
			print(Taxonomy_filename)
			taxdf.to_csv(Taxonomy_filename, sep ="\t" , index=False)



if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", dest="blast", required=True)
	parser.add_argument("-o", "--outdir", dest="outdir", required=True)
	parser.add_argument("-p", "--dump", dest="db", required=True)
	parser.add_argument("-e", "--excluded", dest="excludedseqids", required=True)

	args = parser.parse_args()

	main(args)

