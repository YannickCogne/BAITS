# BAITS

## Recommandation ##
-Install blast and set up to the PATH.<br>
-Install Biopython<br>
-Use Python v2.7<br>


## Usage ##
<ol>
<li>Delete files in repository blast_index; nucl_bdd and prot_bdd.</li>
<li>Copy CDS databases in nucl_bdd and name it "CDS_nameofsample"</li>
<li>Copy protein databases in prot_bdd and name it "T_nameofsample"</li>
 <li> Accession must be the same between CDS and protein databases</li>
<li> Change prot_pep_BioM.tab with peptides you search. ID_contig can correspond to the variable part of accession from reference. <br>
<li>For windows: <br>
BAITS_v1.0vwindows.py "common part of contigs accession from reference" "reference database accession"<br>
For MAC:<br>
Python BAITS_v1.0vmaclinux.py "common part of contigs accession from reference" "reference database accession"<br>

Using no parameters in input will use :<br>
BAITS_v1.0v??.py "Contig_Gammarus_90_" "./bioM_nucl/contigs_rna_gfoss.fasta"<br>
   
  
  </li>
</ol>

## Output ##
All output table are in TSV format <br>
-liste_pept.tab : List of all peptide found by BAITS <br>
-Table_pept.tab: Table with presence and substitute for reference biomarkers <br>
-Table_pept_ident.tab: Table with the identity score between pairwise alignement area and biomarkers <br>
-Table_pept_tryp.tab: Verification of Tryptic site<br>
-Table_prot.tab: Table with proteins corresponding to peptides found or substitutes <br>
