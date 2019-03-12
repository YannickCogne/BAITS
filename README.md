# BAITS

## Recommandation ##
-Install blast and add it to the PATH.<br>
-Install Biopython<br>
-Use Python v2.7<br>


## Usage ##
<ol>
<li>Delete files in repository blast_index, nucl_bdd and prot_bdd.</li>
<li>Copy CDS fasta files in nucl_bdd and name it "CDS_nameofsample"</li>
<li>Copy protein fasta files in prot_bdd and name it "T_nameofsample"</li>
<li>Sequence accession must be the same between CDS and protein fasta files for a given sample</li>
<li>Change prot_pep_BioM.tab with your own peptides. Prot and Nom_files has to be set by the user. ID_contig corresponds to the variable part of the sequence accession from the reference fasta file. <br>

<li> To execute: <br>
For windows: <br>
Python BAITS_v1.0vwindows.py "common part of sequence accession from reference" "reference database file location"<br>
For MAC/LINUX:<br>
Python BAITS_v1.0vmaclinux.py "common part of contigs accession from reference" "reference database file location"<br>
Using no parameters in input will use :<br>
BAITS_v1.0v??.py "Contig_Gammarus_90_" "./bioM_nucl/contigs_rna_gfoss.fasta"<br>
</li>
</ol>

## Output ##
All output table are in TSV format <br>
- liste_pept.tab : List of all peptide found by BAITS <br>
- Table_pept.tab: Table with presence and substitute for reference biomarkers <br>
- Table_pept_ident.tab: Table with the identity score between pairwise alignement area and biomarkers <br>
- Table_pept_tryp.tab: Verification of Tryptic site<br>
- Table_prot.tab: Table with proteins corresponding to peptides found or substitutes <br>
