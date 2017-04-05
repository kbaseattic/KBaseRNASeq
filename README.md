KBaseRNASeq
===========================

This module provides functionality to perform  RNASeq analysis to enable users to quantify gene expression, identify splice junctions and measure differential expression using the Tuxedo suite of tools.

[![Build Status](https://travis-ci.org/kbase/KBaseRNASeq.svg?branch=master)](https://travis-ci.org/kbase/KBaseRNASeq)
Code coverage: (develop branch)
[![Coverage Status](https://coveralls.io/repos/github/kbase/KBaseRNASeq/badge.svg?branch=master)](https://coveralls.io/github/kbase/KBaseRNASeq?branch=master)

https://travis-ci.org/arfathpasha/KBaseRNASeq.svg?branch=master

Build status:</br>
master:  [![Build Status](https://travis-ci.org/arfathpasha/KBaseRNASeq.svg?branch=master)](https://travis-ci.org/arfathpasha/KBaseRNASeq)</br>
staging: [![Build Status](https://travis-ci.org/arfathpasha/KBaseRNASeq.svg?branch=staging)](https://travis-ci.org/arfathpasha/KBaseRNASeq)</br>
develop: [![Build Status](https://travis-ci.org/arfathpasha/KBaseRNASeq.svg?branch=develop)](https://travis-ci.org/arfathpasha/KBaseRNASeq)</br>

Code coverage: (master branch)
[![Coverage Status](https://coveralls.io/repos/github/arfathpasha/KBaseRNASeq/badge.svg?branch=master)](https://coveralls.io/github/arfathpasha/KBaseRNASeq?branch=master)</br>

Examples of tools that will be developed in this module:

1) Associate Reads to RNASeq Sample

2) Build Bowtie2Index

3) Align Reads to Tophat

4) Align Reads to Bowtie2

5) Compute Gene Expression -Cufflinks

6) Merge transcripts to Transcriptome - Cuffmerge

7) Identify Differential Expression  - Cuffdiff 

8) Create ExpressionSeries 

9) Create ExpressionMatrix

10) View Expression Histogram

11) View Alignment Statistics - Pie Chart


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
