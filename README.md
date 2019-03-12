# BAITS

## Recommandation ##
Install blast and set up to the PATH.
Install Biopython
Use Python v2.7


##Usage##
1- Delete files in repository blast_index; nucl_bdd and prot_bdd.
2- Copy CDS databases in nucl_bdd and name it "CDS_nameofsample"
3- Copy protein databases in nucl_bdd and name it "T_nameofsample"
4- Accession must be the same between CDS and protein databases
5- Change prot_pep_BioM.tab with peptides you search. ID_contig can correspond to the variable part of accession from reference.
For windows:
BAITS_v1.0vwindows.py "common part of contigs accession from reference" "reference database accession"
For MAC:
Python BAITS_v1.0vmaclinux.py "common part of contigs accession from reference" "reference database accession"

Using no parameters in input will use :
BAITS_v1.0v??.py "Contig_Gammarus_90_" "./bioM_nucl/contigs_rna_gfoss.fasta"
