import sys
import os 
def main():
	file=sys.argv[1] 
	bestand = open(file, 'r')   
	size = os.path.getsize(file)
	y = file.split('/')
	x = y[-1].split('.')
	samplename = x[0]
	sample = samplename[-2:]
	barcode=x[2]
	Taxonomy_count(bestand, sample, barcode)
	
def Taxonomy_count(bestand, sample, barcode): 
	org1='' 
	level1=''
	taxid1=''
	lineage1=''
	count=0
	count_clusters=0
	result=""
	print("Taxonomic Name\tTaxonomic Level\tTAXID\tClusters\tReads\tLineage")
	for regel in bestand:
		split=regel.split("\t")
		org2=split[3]
		level2=split[2]
		taxid2=split[1]
		lineage2=split[5]
		if len(org2) > 2:
			try: 
				if org2==org1:
					count2 = float(split[6].strip('\n'))
					count_clusters+=1
					count=count+count2
				else:
					if org1!='' and count != 0:
						result=org1+'\t'+str(level1)+'\t'+str(taxid1)+'\t'+str(count_clusters)+'\t'+str(count)+'\t'+lineage1
						print(org1,'\t',level1,'\t',taxid1,'\t',count_clusters,'\t',int(count),'\t',lineage1)
					org1=org2
					count = float(split[6].strip('\n'))
					count_clusters = 1
					level1=level2
					taxid1=taxid2
					lineage1=lineage2
			except ValueError:
				x=0		
	if count !=0 and org1!='' and org1 != result.split("\t")[0]:
		print(org1,'\t',level1,'\t',taxid1,'\t',count_clusters,'\t',int(count),'\t',lineage1)
main()
