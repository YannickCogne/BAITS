import csv
import os
import sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Alphabet import generic_protein
from Bio.Emboss.Applications import NeedleCommandline
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord





def readtablebioM(contig_nom,ref_db):
##### Step1: Preparation of the nucleotide bait #####
	dico={}
	dico_biom=SeqIO.to_dict(SeqIO.parse(ref_db, "fasta"))
	
	with open('prot_pept_bioM.tab', 'rb') as csvfile:
		reader = csv.DictReader(csvfile,delimiter='\t')
		for row in reader:
			if not (dico.has_key(row["Nom_file"])):
				dico[row["Nom_file"]]={}
				dico[row["Nom_file"]]["pept"]=set()
				dico[row["Nom_file"]]["Blast"]={}
				dico[row["Nom_file"]]["Blastseq"]={}
				dico[row["Nom_file"]]["peptind"]={}
				dico[row["Nom_file"]]["pept"].add(Seq(row["BioM"], generic_protein))
				dico[row["Nom_file"]]["peptind"][row["BioM"]]={}
				if dico_biom.has_key(contig_nom+str(row["ID_contig"])):
					dico[row["Nom_file"]]["nucl"]=dico_biom[contig_nom+str(row["ID_contig"])]
				else:
					print "error seq "+contig_nom+str(row["ID_contig"])+" not found"
			else:
				
				dico[row["Nom_file"]]["pept"].add(Seq(row["BioM"], generic_protein))
				dico[row["Nom_file"]]["peptind"][row["BioM"]]={}
	return dico



def blast_nucl(setofind,dico):
#### Step2: Fishing in the ocean of species-assemblies ###
	cwd=os.getcwd()
	with open("./tmp/bioM.fasta", "w") as output_handle:
		for elem in dico.keys():
			seq=dico[elem]["nucl"]
			seq.id=elem
			seq.description=elem
			SeqIO.write(seq, output_handle, "fasta")
		
	for ind in setofind:

		os.system("blastn.exe  -db "+cwd+"\\blast_index\\"+ind+" -query "+cwd+"\\tmp\\bioM.fasta -outfmt 5 -out "+cwd+"\\tmp\\"+ind+"XbioM.xml" )
		
def makedbblast(setofind):
	cwd=os.getcwd()
	
	for elem in setofind:
		presfile=os.listdir(cwd+"\\blast_index\\")
		if not ((elem+".nhr") in presfile):
			os.system("makeblastdb.exe  -in "+cwd+"\\nucl_bdd\\CDS-"+elem+".fasta -out "+cwd+"\\blast_index\\"+elem+" -dbtype \"nucl\"")
		else:
			print "blast database for "+elem+" already done"

def makedbblastprot(setofind):
	cwd=os.getcwd()
	
	for elem in setofind:
		presfile=os.listdir(cwd+"\\blast_index\\")
		if not ((elem+"_prot.phr") in presfile):
			os.system("makeblastdb.exe  -in "+cwd+"\\prot_bdd\\T-"+elem+".fasta -out "+cwd+"\\blast_index\\"+elem+"_prot -dbtype \"prot\"")
		else:
			print "blast database for "+elem+"_prot already done"

						


			
def parsexmlblast(setofind,dico):

	for elem in setofind:
		result_handle = open("./tmp/"+elem+"XbioM.xml")
		blast_records = NCBIXML.parse(result_handle)

		for record in list(blast_records):
			dico[record.query]["Blast"][elem]=list()
			
			for alignment in record.alignments:
				dico[record.query]["Blast"][elem].append(alignment.hit_def.split(" ")[0])

	return dico

def blasttoseq(setofind,dico):
	dicotemp={}
	for elem in setofind:
		dicotemp[elem]=SeqIO.to_dict(SeqIO.parse("./prot_bdd/T-"+elem+".fasta", "fasta"))
	for bioM in dico.keys():
		for elem in dico[bioM]["Blast"].keys():
			dico[bioM]["Blastseq"][elem]=list()
			for seq in dico[bioM]["Blast"][elem]:
				dico[bioM]["Blastseq"][elem].append(dicotemp[elem][seq])
		
	return dico




	
def pairwise_aln(setofind,dico):
	###Step 3: Selection of the best fishes###
	matrix = matlist.blosum62
	for elem in dico.keys():
		for pept in dico[elem]["pept"]:
			for ind in dico[elem]["Blastseq"]:
				trouve=0
				dico[elem]["peptind"][str(pept)][ind]={}
				dico[elem]["peptind"][str(pept)][ind]["prot"]=""
				dico[elem]["peptind"][str(pept)][ind]["pepsub"]=""
				dico[elem]["peptind"][str(pept)][ind]["score"]=np.nan
				dico[elem]["peptind"][str(pept)][ind]["iseq"]=0
				dico[elem]["peptind"][str(pept)][ind]["trypok"]=0
				dico[elem]["peptind"][str(pept)][ind]["ident"]=0
				dico[elem]["peptind"][str(pept)][ind]["align"]=""
				dico[elem]["peptind"][str(pept)][ind]["round"]=0
				dico[elem]["peptind"][str(pept)][ind]["tuple"]=set()
				dico[elem]["peptind"][str(pept)][ind]["Rpept"]=""
				for hit in dico[elem]["Blastseq"][ind]:
					if not trouve:
						align=pairwise2.align.globalds( pept,hit.seq, matrix,-15,-0.5)
					
						for al in align:
						
							

							trypoks=0
							trypoke=0
							start=0
							end=0
							start2=0
							end2=0
							nbmatches=0
							trouvestart=0
							i=0
							for letter in al[0]:
								if not trouvestart:
									if letter !="-":
										trouvestart=1
										start=i
										nbmatches+=1
								else:
									if letter !="-":
										nbmatches+=1
										if nbmatches==len(str(pept)):
											end=i
							
								i+=1
							start2=start
							end2=end
							if not (str(al[1])[start-1]=="K" or str(al[1])[start-1]=="R"):
								while str(al[1])[start2-1]!="K" and str(al[1])[start2-1]!="R" and start2>1:
									start2=start2-1
							else:
								trypoks=1
					
					

							if not (str(al[1])[end]=="K" or str(al[1])[end]=="R"):
								first=1
								while str(al[1])[end2]!="K" and str(al[1])[end2]!="R" and end2<(len(al[1])-1):
									if first:
										first=0
										end2=start2+1
									else:
										end2+=1
							else:
								trypoke=1
							if start2 >0:
								peptsub=str(al[1])[start2-1:end2+1]
							else: 
								peptsub=str(al[1])[start2:end2+1]

							#print al
							tmppept=str(al[1])[start:end+1]
							align2=pairwise2.align.globalxs( pept,tmppept.replace("-",""),-15,-0.5)
							
							if align2:
								score=align2[0][2]
							else:
								score=np.nan
							
							identity=score/len(str(pept))*100
							###Step 4: Validation of the compliance between the best fishes and the biomarker###
							if ((score > dico[elem]["peptind"][str(pept)][ind]["score"]) or(np.isnan(dico[elem]["peptind"][str(pept)][ind]["score"]))) and not (np.isnan(score)) :
								if identity==100 and (trypoke & trypoks):
									trouve=0
								dico[elem]["peptind"][str(pept)][ind]["prot"]=hit.id
								dico[elem]["peptind"][str(pept)][ind]["pepsub"]=peptsub
								dico[elem]["peptind"][str(pept)][ind]["score"]=score
								dico[elem]["peptind"][str(pept)][ind]["iseq"]=(peptsub[1:]==str(pept))
								dico[elem]["peptind"][str(pept)][ind]["trypok"]=(trypoke & trypoks)
								dico[elem]["peptind"][str(pept)][ind]["ident"]=identity
								dico[elem]["peptind"][str(pept)][ind]["align"]=al
								dico[elem]["peptind"][str(pept)][ind]["round"]=1
								toadd=(peptsub[1:],hit.id)
								dico[elem]["peptind"][str(pept)][ind]["tuple"]=set()
								dico[elem]["peptind"][str(pept)][ind]["tuple"].add(toadd)
								dico[elem]["peptind"][str(pept)][ind]["Rpept"]=tmppept
							if ((score == dico[elem]["peptind"][str(pept)][ind]["score"]) and ((trypoke & trypoks) >dico[elem]["peptind"][str(pept)][ind]["trypok"])):
								if identity==100 and (trypoke & trypoks):
									trouve=0
								dico[elem]["peptind"][str(pept)][ind]["prot"]=hit.id
								dico[elem]["peptind"][str(pept)][ind]["pepsub"]=peptsub
								dico[elem]["peptind"][str(pept)][ind]["score"]=score
								dico[elem]["peptind"][str(pept)][ind]["iseq"]=(peptsub[1:]==str(pept))
								dico[elem]["peptind"][str(pept)][ind]["trypok"]=(trypoke & trypoks)
								dico[elem]["peptind"][str(pept)][ind]["ident"]=identity
								dico[elem]["peptind"][str(pept)][ind]["align"]=al
								dico[elem]["peptind"][str(pept)][ind]["round"]=1
								toadd=(peptsub[1:],hit.id)
								dico[elem]["peptind"][str(pept)][ind]["tuple"]=set()
								dico[elem]["peptind"][str(pept)][ind]["tuple"].add(toadd)
								dico[elem]["peptind"][str(pept)][ind]["Rpept"]=tmppept
							if ((score == dico[elem]["peptind"][str(pept)][ind]["score"]) and ((trypoke & trypoks) ==dico[elem]["peptind"][str(pept)][ind]["trypok"])):

								toadd=(peptsub[1:],hit.id)
								dico[elem]["peptind"][str(pept)][ind]["tuple"].add(toadd)

	return dico	









def blast_second(setofind,dico):
	matrix = matlist.blosum62
	dicotemp={}
	###Step 5: Preparation of the protein bait### 
	for elem in setofind:
		dicotemp[elem]=SeqIO.to_dict(SeqIO.parse("./prot_bdd/T-"+elem+".fasta", "fasta"))
	for elem in dico.keys() :
		 for elem2 in dico[elem]["peptind"].keys():
			 for ind in dico[elem]["peptind"][elem2]:
				if  dico[elem]["peptind"][elem2][ind]["ident"]<100 or not dico[elem]["peptind"][elem2][ind]["trypok"] :
					hist=dico[elem]["peptind"][elem2][ind]["pepsub"]
					
					with open("./tmp/"+ind+"_"+elem2+"_bioM.fasta", "w") as output_handle:
						seq=SeqRecord(Seq(elem2,generic_protein),id=elem2,description=elem2)
						SeqIO.write(seq, output_handle, "fasta")								
						for ind2 in dico[elem]["peptind"][elem2].keys():
							prot=dico[elem]["peptind"][elem2][ind2]["prot"]
							if prot:
								seq=dicotemp[ind2][prot]
								SeqIO.write(SeqRecord(seq.seq,id=prot,description=prot), output_handle, "fasta")
					cwd=os.getcwd()
					###Step 6: Re-fishing for uncompliant fishes###
					os.system("blastp.exe  -db "+cwd+"\\blast_index\\"+ind+"_prot -query "+cwd+"\\tmp\\"+ind+"_"+elem2+"_bioM.fasta -outfmt 5 -out "+cwd+"\\tmp\\"+ind+"_"+elem2+"_bioM.xml" )
					listetmp=set()
					result_handle = open("./tmp/"+ind+"_"+elem2+"_bioM.xml")
					blast_records = NCBIXML.parse(result_handle)
					trouve=0
					###Step 7: Selection of potential new fishes###
					for record in list(blast_records):			
						for alignment in record.alignments:
							listetmp.add(alignment.hit_def.split(" ")[0])
					for hitprot in  listetmp:
						hit=dicotemp[ind][hitprot]
						pept=Seq(elem2,generic_protein)
						if not trouve:
							align=pairwise2.align.globalds( pept,hit.seq, matrix,-15,-0.5)

							for al in align:
						
							

								trypoks=0
								trypoke=0
								start=0
								end=0
								start2=0
								end2=0
								nbmatches=0
								trouvestart=0
								i=0
								for letter in al[0]:
									if not trouvestart:
										if letter !="-":
											trouvestart=1
											start=i
											nbmatches+=1
									else:
										if letter !="-":
											nbmatches+=1
											if nbmatches==len(str(pept)):
												end=i
							
									i+=1
								start2=start
								end2=end
								
								if not (str(al[1])[start-1]=="K" or str(al[1])[start-1]=="R"):
									while str(al[1])[start2-1]!="K" and str(al[1])[start2-1]!="R" and start2>1:
										start2=start2-1
								else:
									trypoks=1
					
					

								if not (str(al[1])[end]=="K" or str(al[1])[end]=="R"):
									first=1
									while str(al[1])[end2]!="K" and str(al[1])[end2]!="R" and end2<(len(al[1])-1):
										if first:
											first=0
											end2=start2+1
										else:
											end2+=1
								else:
									trypoke=1
								if start2 >0:
									peptsub=str(al[1])[start2-1:end2+1]
								else: 
									peptsub=str(al[1])[start2:end2+1]

								
								tmppept=str(al[1])[start:end+1]
								align2=pairwise2.align.globalxs( pept,tmppept.replace("-",""),-15,-0.5)
							
								if align2:
									score=align2[0][2]
								else:
									score=np.nan
								
								identity=score/len(str(pept))*100
								###Step 8: Validation of the compliance between the new fishes and the biomarker###
								###Step 9: Suggestion of a substitute biomarker###
								if ((score > dico[elem]["peptind"][str(pept)][ind]["score"]) or(np.isnan(dico[elem]["peptind"][str(pept)][ind]["score"]))) and not (np.isnan(score)) :
									if identity==100 and (trypoke & trypoks):
										trouve=0
									dico[elem]["peptind"][str(pept)][ind]["prot"]=hit.id
									dico[elem]["peptind"][str(pept)][ind]["pepsub"]=peptsub
									dico[elem]["peptind"][str(pept)][ind]["score"]=score
									dico[elem]["peptind"][str(pept)][ind]["iseq"]=(peptsub[1:]==str(pept))
									dico[elem]["peptind"][str(pept)][ind]["trypok"]=(trypoke & trypoks)
									dico[elem]["peptind"][str(pept)][ind]["ident"]=identity
									dico[elem]["peptind"][str(pept)][ind]["align"]=al
									dico[elem]["peptind"][str(pept)][ind]["round"]=2
									toadd=(peptsub[1:],hit.id)
									dico[elem]["peptind"][str(pept)][ind]["tuple"]=set()
									dico[elem]["peptind"][str(pept)][ind]["tuple"].add(toadd)
									dico[elem]["peptind"][str(pept)][ind]["Rpept"]=tmppept
								if ((score == dico[elem]["peptind"][str(pept)][ind]["score"]) and ((trypoke & trypoks) >dico[elem]["peptind"][str(pept)][ind]["trypok"])):
									if identity==100 and (trypoke & trypoks):
										trouve=0
									dico[elem]["peptind"][str(pept)][ind]["prot"]=hit.id
									dico[elem]["peptind"][str(pept)][ind]["pepsub"]=peptsub
									dico[elem]["peptind"][str(pept)][ind]["score"]=score
									dico[elem]["peptind"][str(pept)][ind]["iseq"]=(peptsub[1:]==str(pept))
									dico[elem]["peptind"][str(pept)][ind]["trypok"]=(trypoke & trypoks)
									dico[elem]["peptind"][str(pept)][ind]["ident"]=identity
									dico[elem]["peptind"][str(pept)][ind]["align"]=al
									dico[elem]["peptind"][str(pept)][ind]["round"]=2
									toadd=(peptsub[1:],hit.id)
									dico[elem]["peptind"][str(pept)][ind]["tuple"]=set()
									dico[elem]["peptind"][str(pept)][ind]["tuple"].add(toadd)
									dico[elem]["peptind"][str(pept)][ind]["Rpept"]=tmppept

								if ((score == dico[elem]["peptind"][str(pept)][ind]["score"]) and ((trypoke & trypoks) ==dico[elem]["peptind"][str(pept)][ind]["trypok"])):

									toadd=(peptsub[1:],hit.id)
									dico[elem]["peptind"][str(pept)][ind]["tuple"].add(toadd)

	return dico
def Tab1(dico,setofind):
	toprint=list()
	for i in range(0,len(setofind)+2,1):
		toprint.append("")
	
	listind=list(setofind)
	for i in range(2,len(setofind)+2,1):
		toprint[i]+=listind[i-2]
	for elem in dico.keys():
		for pept in dico[elem]["pept"]:
			toprint[0]+="\t"+elem
			toprint[1]+="\t"+str(pept)
			for i in range(2,len(setofind)+2,1):
				if  (dico[elem]["peptind"][str(pept)][listind[i-2]]["ident"]==100) and (dico[elem]["peptind"][str(pept)][listind[i-2]]["trypok"]) :
					toprint[i]+="\tOK"
				else:
					tmpstr=""
					tmpset=set()
					for item in dico[elem]["peptind"][str(pept)][listind[i-2]]["tuple"]:
						tmpset.add(item[0])
					for item in tmpset:
						tmpstr+=";"+item
					toprint[i]+="\t"+tmpstr[1:]
	with open("Table_pept.tab", 'w') as outfile	:			
		for ligne in toprint:
			outfile.write(ligne+"\n")
	
def Tab2(dico,setofind):
	toprint=list()
	for i in range(0,len(setofind)+2,1):
		toprint.append("")
	
	listind=list(setofind)
	for i in range(2,len(setofind)+2,1):
		toprint[i]+=listind[i-2]
	for elem in dico.keys():
		for pept in dico[elem]["pept"]:
			toprint[0]+="\t"+elem
			toprint[1]+="\t"+str(pept)
			for i in range(2,len(setofind)+2,1):
				toprint[i]+="\t"+str(dico[elem]["peptind"][str(pept)][listind[i-2]]["ident"])
					
	with open("Table_pept_ident.tab", 'w') as outfile	:			
		for ligne in toprint:
			outfile.write(ligne+"\n")
	
		
def Tab3(dico,setofind):
	toprint=list()
	for i in range(0,len(setofind)+2,1):
		toprint.append("")
	
	listind=list(setofind)
	for i in range(2,len(setofind)+2,1):
		toprint[i]+=listind[i-2]
	for elem in dico.keys():
		for pept in dico[elem]["pept"]:
			toprint[0]+="\t"+elem
			toprint[1]+="\t"+str(pept)
			for i in range(2,len(setofind)+2,1):
				toprint[i]+="\t"+str(dico[elem]["peptind"][str(pept)][listind[i-2]]["trypok"])
					
	with open("Table_pept_tryp.tab", 'w') as outfile:				
		for ligne in toprint:
			outfile.write(ligne+"\n")
		
def Tab4(dico,setofind):
	toprint=list()
	for i in range(0,len(setofind)+2,1):
		toprint.append("")
	
	listind=list(setofind)
	for i in range(2,len(setofind)+2,1):
		toprint[i]+=listind[i-2]
	for elem in dico.keys():
		for pept in dico[elem]["pept"]:
			toprint[0]+="\t"+elem
			toprint[1]+="\t"+str(pept)
			for i in range(2,len(setofind)+2,1):

				tmpstr=""
				tmpset=set()
				for item in dico[elem]["peptind"][str(pept)][listind[i-2]]["tuple"]:
					tmpset.add(item[1])
				for item in tmpset:
					tmpstr+=";"+item
				toprint[i]+="\t"+tmpstr[1:]
					
	with open("Table_prot.tab", 'w') as outfile	:			
		for ligne in toprint:
			outfile.write(ligne+"\n")	
			
def list1(dico,setofind):
	dicotmp={}
	line=""
	for elem in dico.keys() :
		 for elem2 in dico[elem]["peptind"].keys():
			 for ind in dico[elem]["peptind"][elem2]:
				for res in dico[elem]["peptind"][elem2][ind]["tuple"]:
					if dicotmp.has_key(res[0]):
						dicotmp[res[0]]["prot"].add(res[1])
						dicotmp[res[0]]["orgn"].add(ind)
						dicotmp[res[0]]["ref"].add(elem2)
						dicotmp[res[0]]["protref"].add(elem)
						dicotmp[res[0]]["round"].add(dico[elem]["peptind"][elem2][ind]["round"])				
					else:
						dicotmp[res[0]]={}
						dicotmp[res[0]]["prot"]=set()
						dicotmp[res[0]]["protref"]=set()
						dicotmp[res[0]]["orgn"]=set()
						dicotmp[res[0]]["ref"]=set()
						dicotmp[res[0]]["round"]=set()
						dicotmp[res[0]]["TrypOK"]=0
						dicotmp[res[0]]["ident"]=0
						dicotmp[res[0]]["same"]=0
						dicotmp[res[0]]["prot"].add(res[1])
						dicotmp[res[0]]["orgn"].add(ind)
						dicotmp[res[0]]["ref"].add(elem2)
						dicotmp[res[0]]["protref"].add(elem)
						dicotmp[res[0]]["TrypOK"]=dico[elem]["peptind"][elem2][ind]["trypok"]
						dicotmp[res[0]]["ident"]=dico[elem]["peptind"][elem2][ind]["ident"]
						dicotmp[res[0]]["round"].add(dico[elem]["peptind"][elem2][ind]["round"])
						dicotmp[res[0]]["same"]=dico[elem]["peptind"][elem2][ind]["iseq"]
	i=0
	with open("liste_pept.tab", 'w') as outfile	:
		for elem in dicotmp:
			
			line="pept"
			if i==0:
				for elem2 in dicotmp[elem].keys():
					line+="\t"+str(elem2)
				outfile.write(line+"\n")
			line=elem
			
			for elem2 in dicotmp[elem].keys():
				if not(type(dicotmp[elem][elem2])==type(set())):
					line+="\t"+str(dicotmp[elem][elem2])
				else:
					temp=""
					for item in dicotmp[elem][elem2]:
						temp+=";"+str(item)
					line+="\t"+temp[1:]
			outfile.write(line+"\n")
			i+=1	
		
def main():
	if (len(sys.argv)!=3):
		##To modify with the common name of contigs"
		contig_nom="Contig_Gammarus_90_"
		##To modify with the file localisation of fasta reference files
		ref_db="./bioM_nucl/contigs_rna_gfoss.fasta"
	else :
		contig_nom=sys.argv[1]
		ref_db=sys.argv[2]
	dico={}
	##### Step1: Preparation of the nucleotide bait #####
	dico=readtablebioM(contig_nom,ref_db)
	setofind=set()
	
	for elem in os.listdir("./prot_bdd/"):
		setofind.add( elem.split("-")[1].split(".")[0])

	makedbblast(setofind)
	makedbblastprot(setofind)
	blast_nucl(setofind,dico)
	dico=parsexmlblast(setofind,dico)
	#### Step2: Fishing in the ocean of species-assemblies ###
	dico=blasttoseq(setofind,dico)
	###Step 3: Selection of the best fishes###
	###Step 4: Validation of the compliance between the best fishes and the biomarker###
	dico=pairwise_aln(setofind,dico)
	###Step 5: Preparation of the protein bait### 
	###Step 6: Re-fishing for uncompliant fishes###
	###Step 7: Selection of potential new fishes###
	###Step 8: Validation of the compliance between the new fishes and the biomarker###
	###Step 9: Suggestion of a substitute biomarker###
	
	dico=blast_second(setofind,dico)

	Tab1(dico,setofind)
	Tab2(dico,setofind)
	Tab3(dico,setofind)
	Tab4(dico,setofind)
	list1(dico,setofind)
main()