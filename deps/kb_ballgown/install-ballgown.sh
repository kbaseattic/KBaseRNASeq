#!bin/bash

echo Starting R and Bioconductor update, followed by Ballgown installation

echo "deb http://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
apt-get update
rm -rf /usr/lib/R
rm -rf /usr/local/lib/R

echo "##########" About to install a new R "###########"

DEBIAN_FRONTEND=noninteractive apt-get -y -o Dpkg::Options::="--force-confdef" -o Dpkg::Options::="--force-confnew" install r-base

echo "##########" About to load bioconductor and ballgown "###########"

R -q -e 'chooseCRANmirror(ind=52); install.packages(c("getopt")); source("http://bioconductor.org/biocLite.R"); biocLite("ballgown")' 

echo "##########" Finished loading bioconductor and ballgown "###########"

echo Completed: R and Bioconductor update, followed by Ballgown installation

