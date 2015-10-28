if [ -z "$(which Rscript)" ]; then
  apt-get install r-base
fi
R -q -e 'source("http://bioconductor.org/biocLite.R")'
R -q -e 'biocLite("BiocUpgrade")'
R -q -e 'if (!require("cummeRbund")) biocLite("cummeRbund")'
