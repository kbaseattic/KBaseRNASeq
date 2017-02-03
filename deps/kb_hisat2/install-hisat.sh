#!/bin/bash

dest=${TARGET-/usr/}/bin
echo "using $dest as installation directory";
mkdir -p $dest

# downlownload version
VERSION='2.0.4'
rm -rf hisat-${VERSION}*
wget "ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-${VERSION}-Linux_x86_64.zip"
unzip hisat2-${VERSION}-Linux_x86_64.zip
rm hisat2-${VERSION}-Linux_x86_64.zip
# compile and copy binaries
cd hisat2-${VERSION}
#make
#mkdir -p dest/bin
cp hisat2-inspect-s hisat2_test_HLA_genotyping.py hisat2_build_genotype_genome.py hisat2_simulate_reads.py hisat2_extract_snps_haplotypes_UCSC.py hisat2_extract_HLA_vars.py hisat2-build-l-debug hisat2-inspect-s-debug hisat2_extract_exons.py hisat2_test_BRCA_genotyping.py hisat2_extract_splice_sites.py hisat2-build-l hisat2-inspect-l-debug hisat2_genotype.py hisat2-inspect hisat2-align-l hisat2_extract_snps_haplotypes_VCF.py hisat2 extract_exons.py hisat2-align-s-debug hisat2-build-s hisat2-align-s extract_splice_sites.py hisat2-inspect-l hisat2-build-s-debug hisat2-build hisat2-align-l-debug $dest
cd ..
rm -rf hisat-${VERSION}
ver=`hisat2 --version`
echo "${ver} is successfully installed"
