#!/bin/bash
dest=${TARGET-/usr/}/bin
echo "using $dest as installation directory";
mkdir -p $dest

# downlownload version
VERSION='2.1.1'
rm -rf tablemaker${VERSION}*
wget -O tablemaker${VERSION}.tar.gz https://ndownloader.figshare.com/files/3193031
tar -xzvf tablemaker${VERSION}.tar.gz
rm tablemaker${VERSION}.tar.gz
# compile and copy binaries
cd tablemaker-${VERSION}.Linux_x86_64
#make
cp `find . -maxdepth 1 -perm -111 -type f` ${dest}
cd ../
rm -rf tablemaker-${VERSION}.Linux_x86_64
#ver=`tablemaker --version`
echo "tablemaker is successfully installed"
