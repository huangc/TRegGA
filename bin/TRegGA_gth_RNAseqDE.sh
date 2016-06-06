#! /usr/bin/bash
# This is to look up genes found by gth with the differentially expressed (DE) genes from RNA-Seq results

## General parameters
WORK_DIR=$(pwd)
TRegGA_DIR=/projects/huangcy/MYGIT/TRegGA

RNASEQ_DIR=${TRegGA_DIR}/doc
RNASEQ_TABLE=131108_CY_RSEM_EBSeq_Genes_DE_Deliverable.txt

# SAMPLE="IRBB7 DV86 IRBB62"
TARGET_SEQNAME=OsaXA7
SAMPLE="IRBB7 DV86 IRBB62 XA7 XA7A XA7B"


##--------------------------------------------------------------
## 
cd ${TRegGA_DIR}/run
mkdir -p gth_RNAseqDE
cd gth_RNAseqDE
\cp ${RNASEQ_DIR}/${RNASEQ_TABLE} .

# RNA-Seq DE gene analysis result table has the following tab-delimited features:
# GENE_ID Ensembl_ID      IRGSP1_MSU_ID   Transcript_ID   OsjChr  Start   End     Strand  Length  
# L11-a_count     L11-a_FPKM      L11-b_count     L11-b_FPKM      L12_count       L12_FPKM        L22_count        L22_FPKM
# L13_count       L13_FPKM        L23_count       L23_FPKM        Lx2_vs_Lx1_PostFC       Lx3_vs_Lx1_PostFC     Lx3_vs_Lx2_PostFC
# Gene_Note       R_Struc GO_Biological_Process   GO_Cellular_Component   GO_Molecular_Function   InterPro

# Noted that lines in the tab-delimited file exported from Excel are ended with "\r", and that needs to be changed to "\n"
tr "\r" "\n" < ${RNASEQ_TABLE} > RNASEQ.table

for SYNONYM in $SAMPLE
do
\cp ${TRegGA_DIR}/assembly/rfguided/${SYNONYM}-on-${TARGET_SEQNAME}/EVALUATION/gth.summary gth.summary.${SYNONYM}-${TARGET_SEQNAME}
\cp ${TRegGA_DIR}/assembly/rfguided/${SYNONYM}-on-${TARGET_SEQNAME}/EVALUATION/gth.OsjPRT-on-${TARGET_SEQNAME} .
\cp ${TRegGA_DIR}/assembly/rfguided/${SYNONYM}-on-${TARGET_SEQNAME}/EVALUATION/gth.OsjPRT-on-${SYNONYM}-${TARGET_SEQNAME} .

## Note that gth found genes in gth.OsjPRT-on-${SYNONYM}-${TARGET_SEQNAME} has the following annotation:
# Protein Sequence: file=/projects/huangcy/MYGIT/TRegGA/reference/rice_japonica/OsjPRT, description=OS06T0672700-01 
# pep:known chromosome:IRGSP-1.0:6:27923541:27926280:-1 gene:OS06G0672700 transcript:OS06T0672700-01 
# description:"Note\x3dConserved hypothetical protein., Transcript_evidence\x3dAK101470 (DDBJ, Best hit), 
# ORF_evidence\x3dN\P_001058337.1 (RefSeq), NIAS_FLcDNA\x3dJ033040K01,"

# Extract GENE_ID in gth, such as "OS06G0672700", which should match the GENE_ID feature in the ${RNASEQ_TABLE}
grep "^Protein Sequence" gth.OsjPRT-on-${SYNONYM}-${TARGET_SEQNAME} |\
    awk -F"gene:" '{ print $2 }' |\
    awk -F" " '{ print $1 }' | uniq > gth.genelist.${SYNONYM}-${TARGET_SEQNAME}

# Find matching GENE_ID feature in the ${RNASEQ_TABLE}
genelist=`tr "\n" " " < gth.genelist.${SYNONYM}-${TARGET_SEQNAME}`
\rm -f DEgene_${SYNONYM}
for k in $genelist
do
grep "^$k" RNASEQ.table >> DEgene_${SYNONYM}
done
head -1 RNASEQ.table | cat - DEgene_${SYNONYM} > gth.RNAseqDE.${SYNONYM}-${TARGET_SEQNAME}
done

# Remove temp files
\rm -f RNASEQ.table
\rm -f DEgene_*

