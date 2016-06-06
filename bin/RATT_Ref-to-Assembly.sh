#!/usr/bin/bash
# This is to transfer gene annotation from reference to an assembly sequence.

## General parameters

TRegGA_DIR=/projects/huangcy/MYGIT/TRegGA
PAGIT_DIR=/usr/local/src/NGS-DIR//PAGIT/PAGIT
WORK_DIR=${TRegGA_DIR}/run

# SYNONYM=LABEL # Replace LABEL with an appropriate cultivar name; e.g. IRBB7
# SCAFFOLDS=${DeNOVO_DIR}/${SYNONYM}-soap.scafSeq # the regular SOAP scaffolds

## REFSEQ can be the Embl file of the TRegGA rfguided-assembled sequence and annotation that are to be transferred to the QUERYSEQ
TARGET_DIR=${TRegGA_DIR}/targets # TARGET_SEQs are in this folder
TARGET_SEQNAME=OsaXA7 # Replace GUIDESEQ with an appropriate target name; e.g. Nb58k or OsaXA7
TARGET_SEQ=${TARGET_SEQNAME}.fasta
TARGET_EMBL=${TARGET_SEQNAME}.embl
TARGET_GFF3=${TARGET_SEQNAME}.gff3
REFSEQ=${TARGET_DIR}/${TARGET_EMBL}
REFSEQNAME=${TARGET_SEQNAME}
REFFASTA=${TARGET_DIR}/${TARGET_SEQ}

# Alternatively, REFSEQ can be the whole reference chromosome  
# REFSEQ=${TRegGA_DIR}/targets/OsaCHR6.embl
# REFSEQNAME=OsaCHR6
# In case whole reference chromosome is used, the fasta file needs to be generated from the embl file.
cd ${TRegGA_DIR}/targets
# python ${TRegGA_DIR}/targets/getTarget.py ${REFSEQ} ${REFSEQNAME}
seqret -sequence ${REFSEQNAME}.embl -feature -snucleotide T -supper1 -sformat1 embl -osformat2 fasta -outseq ${REFSEQNAME}.fasta
# REFFASTA=${TRegGA_DIR}/targets/${REFSEQNAME}.fasta

# QUERYSEQNAME in regular TRegGA is RATT_QUERY_SEQNAME=${SYNONYM}-${TARGET_SEQNAME}
# QUERYSEQ=${TRegGA_DIR}/assembly/rfguided/${SYNONYM}-on-${TARGET_SEQNAME}/EVALUATION/${SYNONYM}-${TARGET_SEQNAME}.fa
# QUERYSEQNAME=${SYNONYM}-${TARGET_SEQNAME}
QUERYSEQ=/projects/huangcy/MYGIT/TRegGA/run/DJ123/OsdDNA-scaffold_10.fa
QUERYSEQNAME=DJ123_scaffold10

## PAGIT
# ------------------------------------------------------------------------------
# The PAGIT package includes a script that must be placed in the user's .bashrc
# file for PAGIT to work properly.
source ${PAGIT_DIR}/sourceme.pagit

###
################################################################################
### ! Typically there would be no need for further editing below this line ! ###
##
################################################################################
#
## Setup parameters for RATT
#
RATT_DIR=${PAGIT_DIR}/RATT
RATT_OUTPUT=RATT
# RATT_QUERY_SEQNAME=${SYNONYM}-${TARGET_SEQNAME}
RATT_QUERY_SEQNAME=${QUERYSEQNAME}
RATT_QUERY_SEQ=${RATT_QUERY_SEQNAME}.fa
RATT_QUERY_EMBL=${RATT_QUERY_SEQNAME}.embl
#
RTO_d=embl # Directory name with embl-annotation files.
RTO_q=${RATT_QUERY_SEQ} # A multifasta file to which the annotation will be mapped.
RTO_r=RATT # The prefix you wish to give to each result file
RTO_t=Strain #  Transfer type, one of Assembly, Strain, Species, Multiple, Free.
RATT_OPTIONS="${RTO_d} ${RTO_q} ${RTO_r} ${RTO_t}"
#
RATT_REPORT=${RATT_OUTPUT}/${RTO_r}.${RATT_QUERY_SEQNAME}.Report.txt
#
EVALUATION_DIR=EVALUATION
GTH_PRTDIR=${TRegGA_DIR}/reference/rice_japonica
GTH_PRTDB=OsjPRT

##------------------------------------------------
## RATT to transfer the annotation from reference genome to the new assembly using RATT (Rapid Annotation Transfer Tool) from PAGIT
#
cd ${WORK_DIR}
mkdir -p ratt_${QUERYSEQNAME}-to-${REFSEQNAME}
cd ratt_${QUERYSEQNAME}-to-${REFSEQNAME}

# make embl directory and copy the embl annotation file to here.
#
mkdir -p ${RATT_OUTPUT}
mkdir -p ${RATT_OUTPUT}/${RTO_d}
\cp ${REFSEQ} ${RATT_OUTPUT}/${RTO_d}/${REFSEQNAME}.embl
\cp ${QUERYSEQ} ${RATT_OUTPUT}/${RATT_QUERY_SEQ}

# Modify ${QUERYSEQ} such that the seqname reflects the ${QUERYSEQNAME}
sed -i -r "s/>.*$/>${QUERYSEQNAME}/;" ${RATT_OUTPUT}/${RATT_QUERY_SEQ}

# Run RATT
#
cd ${RATT_OUTPUT}
${RATT_DIR}/start.ratt.sh ${RATT_OPTIONS}

## Evaluation
#
cd ${WORK_DIR}/ratt_${QUERYSEQNAME}-to-${REFSEQNAME}
mkdir -p ${EVALUATION_DIR}

echo "This directory contains final results and evaluations of the workflow." > ${EVALUATION_DIR}/0README
\cp ${RATT_OUTPUT}/${QUERYSEQNAME}.fa ${EVALUATION_DIR}
\cp ${RATT_OUTPUT}/${RTO_r}.${QUERYSEQNAME}.final.embl ${EVALUATION_DIR}/${QUERYSEQNAME}.embl
\cp ${RATT_OUTPUT}/${RTO_r}.${QUERYSEQNAME}.Report.gff ${EVALUATION_DIR}/${QUERYSEQNAME}.gff3
\cp ${REFSEQ} ${EVALUATION_DIR}/${REFSEQNAME}.embl
\cp ${RATT_OUTPUT}/nucmer.${RTO_r}.snp ${EVALUATION_DIR}/${QUERYSEQNAME}.snp

#
cd ${EVALUATION_DIR}
gth -genomic ${REFFASTA} -protein ${GTH_PRTDIR}/${GTH_PRTDB} > gth.${GTH_PRTDB}-on-${REFSEQNAME}
gth -genomic ${QUERYSEQNAME}.fa -protein ${GTH_PRTDIR}/${GTH_PRTDB} > gth.${GTH_PRTDB}-on-${QUERYSEQNAME}

#
echo "CDS annotated in ${REFSEQNAME}.embl:" > gth.summary
echo "" >> gth.summary
egrep "CDS" ${REFSEQNAME}.embl >> gth.summary
echo "" >> gth.summary
echo "CDS transferred to ${QUERYSEQNAME}.embl:" >> gth.summary
echo "" >> gth.summary
egrep "CDS" ${QUERYSEQNAME}.embl >> gth.summary
echo "" >> gth.summary
echo "PGS annotated in gth.${GTH_PRTDB}-on-${REFSEQNAME} (might contain partial matches):" >> gth.summary
echo "" >> gth.summary
egrep "^  PGS" gth.${GTH_PRTDB}-on-${REFSEQNAME} >> gth.summary
echo "" >> gth.summary
echo "PGS annotated in gth.${GTH_PRTDB}-on-${QUERYSEQNAME}" >> gth.summary
echo "(might contain partial matches and matches outside the target-transferred gene models):" >> gth.summary
echo "" >> gth.summary
egrep "^  PGS" gth.${GTH_PRTDB}-on-${QUERYSEQNAME} >> gth.summary
echo "" >> gth.summary

#
makeblastdb -in ${REFFASTA} -dbtype nucl -out TRGT -parse_seqids
blastn -db TRGT -query ${QUERYSEQNAME}.fa -outfmt 6 > blastn.out

echo "" >> gth.summary
echo "To visually examine the reference-guided assembly relative to the target sequence," >> gth.summary
echo "you could run:" >> gth.summary
echo "" >> gth.summary
echo "cd ${EVALUATION_DIR}; act ${REFSEQNAME}.embl blastn.out ${QUERYSEQNAME}.embl" >> gth.summary
echo "" >> gth.summary

