#!/usr/bin/env Rscript

# ballgown_fpkmgenematrix.R - first shot at using ballgown to compute 
# gene-level differential expression matrix

# example call:
# ./ballgown_fpkmgenematrix.R --ballgown_dir . --sample_dir_pat "_rep[12]_stringtie" --experiment_groups "1100" --out_csvfile "testout.csv"

cat( "############## Hello this ballgown_fpkmgenematrix.R ##############\n" )

cat( "############## getwd() is ", getwd(), " ###################\n" )

suppressMessages( require("getopt") )
suppressMessages( require( ballgown ) )
#suppressMessages(require (rjson))

options(showWarnCalls=FALSE)
options(showErrorCalls=FALSE)

spec = matrix( c(
                 'ballgown_dir',      'd', 1, "character",
                 'sample_dir_pat',    'p', 1, "character",  # re pattern for selecting subdirectories
                 'experiment_groups', 'e', 1, "character",  # character string of 0s and 1s indicating group membership of each column
                 'out_csvfile',       'o', 1, "character",  # output differential gene expression CSV file name
                 'help',              'h', 0, "logical"
                ), 
              byrow=TRUE, ncol=4 );
opt = getopt( spec )

if ( ! is.null( opt$help ) ) {
    cat( getopt( spec, usage=TRUE ) );
    q( status = 1 );
}
cat( "############## ls -lR ################\n" )
system( paste( "ls -lR", opt$ballgown_dir ) ) 

cat( "############## here is the spec matrix of arguments ##############\n" )
print( spec )
cat( "############## about to create ballgown ", opt$ballgown_dir, ", ", opt$sample_dir_pat, " #############\n" )

bg <- ballgown( dataDir = opt$ballgown_dir, samplePattern = opt$sample_dir_pat, meas="all" )

cat( "############## here is bg str ##############\n" )
print( str( bg ) )

# convert string "111000" integer into vector c(1,1,1,0,0,0)
experiment_group_vector <- as.integer( strsplit( opt$experiment_groups, "")[[1]] )   

cat( "############## here is experiment_group_vector ##############\n" )
print( experiment_group_vector )

cat( "############## about to pData ##############\n" )

# map group membership to column names
pData(bg) <- data.frame( id = sampleNames(bg), group = experiment_group_vector )

cat( "############## about to stattest ##############\n" )

# create gene-level differential expression table
gene_diff_ex_tab <- stattest( bg, feature='gene', meas='FPKM', covariate='group' )

cat( "############## about to write.table ##############\n" )

write.table( gene_diff_ex_tab[,-1], file=opt$out_csvfile, sep="\t", row.names=F )  # remove first column before saving

cat( "############## kthxbye ##############\n" )

q( save="no" )

