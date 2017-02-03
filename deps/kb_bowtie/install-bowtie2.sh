#!/bin/bash

dest=${TARGET-/usr/}/bin
echo "using $dest as installation directory";
mkdir -p $dest

# downlownload version
VERSION='2.2.9'
rm -rf bowtie2-${VERSION}*
wget https://github.com/BenLangmead/bowtie2/releases/download/v${VERSION}/bowtie2-${VERSION}-linux-x86_64.zip
unzip bowtie2-${VERSION}-linux-x86_64.zip
# compile and copy binaries
cd bowtie2-${VERSION}
#make
#mkdir -p dest/bin
cp bowtie2 bowtie2-align-l bowtie2-align-l-debug bowtie2-align-s bowtie2-align-s-debug bowtie2-build bowtie2-build-l bowtie2-build-l-debug bowtie2-build-s bowtie2-build-s-debug bowtie2-inspect bowtie2-inspect-l bowtie2-inspect-l-debug bowtie2-inspect-s bowtie2-inspect-s-debug $dest
cd ..
rm -rf bowtie2-${VERSION}
ver=`bowtie2 --version`
echo "${ver} is successfully installed"
