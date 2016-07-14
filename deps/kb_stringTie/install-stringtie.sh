#!/bin/bash
dest=${TARGET-/usr/}/bin
echo "using $dest as installation directory";
mkdir -p $dest

# downlownload version
VERSION='1.2.3'
rm -rf stringtie-${VERSION}*
wget "http://ccb.jhu.edu/software/stringtie/dl/stringtie-${VERSION}.Linux_x86_64.tar.gz"
tar -xzvf stringtie-${VERSION}.Linux_x86_64.tar.gz
rm stringtie-${VERSION}.Linux_x86_64.tar.gz
# compile and copy binaries
cd stringtie-${VERSION}.Linux_x86_64
#make
cp `find . -maxdepth 1 -perm -111 -type f` ${dest}
cd ../
rm -rf stringtie-${VERSION}.Linux_x86_64
ver=`stringtie --version`
echo "${ver} is successfully installed"
