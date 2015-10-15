KBaseRNASeq
===========================

This module provides functionality to perform  RNASeq analysis to enable users to quantify gene expression, identify splice junctions and measure differential expression using the Tuxedo suite of tools.



Examples of tools that will be developed in this module:

1) Build Bowtie2Index

2) Align Reads to Tophat

3) Align Reads to Bowtie2

4) Compute Gene Expression -Cufflinks

5) Merge transcripts to Transcriptome - Cuffmerge

6) Identify Differential Expression  - Cuffdiff 

7) Create ExpressionSeries 

8) Create ExpressionMatrix

9) View Expression Histogram

10) View Alignment Statistics - Pie Chart


Notes for development:
=====================

1) Delete the PATHONPATH in the 'bin/run_KBaseRNASeq.sh'

2) . /kb/dev_container/user-env.sh

3) export PYTHONPATH='kb/deployment/lib'

4) insert 'token=token.rstrip()' in line526 within lib/biokbase/RNASeq/KBaseRNASeq.py

5) make deploy

6) testing the command:  ./bin/run_KBaseRNASeq.sh test/script_test/input.json output.txt /mnt/project/mytoken.txt

References:
============

Trapnell C,et al. 2009 TopHat: discovering splice junctions with RNA-Seq.Bioinformatics

Trapnell C,et al. 2010 Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell
differentiation.Nature Biotechnology

Kim D, et al. 2011 TopHat2: accurate alignment of transcriptomes in the presence of insertions, deletions and gene fusions.
Genome Biology
