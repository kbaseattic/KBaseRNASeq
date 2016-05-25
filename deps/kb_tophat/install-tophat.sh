#!/bin/bash

dest=${TARGET-/usr/}/bin
echo "using $dest as installation directory";
mkdir -p $dest

# downlownload version
VERSION='2.1.1'
rm -rf tophat-${VERSION}*
wget "https://ccb.jhu.edu/software/tophat/downloads/tophat-${VERSION}.Linux_x86_64.tar.gz"
tar -xzvf tophat-${VERSION}.Linux_x86_64.tar.gz
rm tophat-${VERSION}.Linux_x86_64.tar.gz
# compile and copy binaries
cd tophat-${VERSION}.Linux_x86_64
#make
#mkdir -p dest/bin
cp bam2fastx bam_merge bed_to_juncs contig_to_chr_coords fix_map_ordering gtf_juncs gtf_to_fasta juncs_db long_spanning_reads map2gtf prep_reads sam_juncs samtools_0.1.18 segment_juncs sra_to_solid tophat tophat2 tophat-fusion-post tophat_reports $dest
cd ..
rm -rf tophat-${VERSION}.Linux_x86_64
ver=`tophat --version`
echo "${ver} is successfully installed"
