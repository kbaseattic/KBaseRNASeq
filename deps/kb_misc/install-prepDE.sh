#!/bin/bash
# get a copy of prepDE.py for expression matrix construction
#
dest=${TARGET-/usr/}/bin
echo "using $dest as installation directory";
mkdir -p $dest

cd $dest
curl 'https://ccb.jhu.edu/software/stringtie/dl/prepDE.py' -O
chmod a+rx prepDE.py
echo "prepDE.py is successfully installed"

