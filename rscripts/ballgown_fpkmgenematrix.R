#!/usr/bin/env Rscript
#
# ballgown_fpkmgenematrix.R - first shot at using ballgown to compute gene-level differential expression matrix
# 
#
# example call:
#
# ./ballgown_fpkmgenematrix.R  --sample_dir_group_table samples_groups.txt --output_dir ballgown_out --output_csvfile testout.csv
#
#          samples_groups.txt contains lines of two columns 
#                                - sample subdirectory name (used as sample name
#                                - corresponding group (0, or 1)
#
#                example
#                         hy5_rep1_stringtie 0
#                         hy5_rep2_stringtie 0
#                         WT_rep1_stringtie  1
#                         WT_rep2_stringtie  1
#                (sub directory column can also have the full pathname specified)
#
#          testout.csv file is gene-level FPKM output 
#
#                example
#
#                         "id"	"fc"	"pval"	"qval"
#                         "AT1G01010"	1	NA	NA
#                         "AT1G01020"	1	NA	NA
#                         "AT1G01030"	1334.97985720067	0.440550490118355	0.787884152938117
#                         "AT1G01040"	0.190911435650667	0.181401258304825	0.787884152938117
#                         "AT1G01046"	1			NA			NA
#                         "AT1G01050"	1			NA			NA
#                         "AT1G01060"	0.293968217480655	0.212709415508708	0.787884152938117
#               


                         #############################
                         # Subroutines and functions #
                         #############################


dmesg <- function( msg, ... ) 
  { cat( "##################", msg, ..., "##################\n" ) }

get_and_check_samples_groups <- function( sgfile )     # load file which contains sample dir names and corresponding groups
  {
   # put in checks here:  directory existance, group validity, etc
   return( read.table( sgfile, col.names = c( "sample_dir", "group" ) ) )
  }



                              ################
                              # Main program #
                              ################

dmesg( "Hello this ballgown_fpkmgenematrix.R" )

dmesg( "getwd() is", getwd() )

suppressMessages( require("getopt") )
suppressMessages( require( ballgown ) )
#suxppressMessages(require (rjson))

options( showWarnCalls = FALSE )
options( showErrorCalls = FALSE )

# parse command line arguments/options

option_tab = matrix( c(
                       #'input_dir',             'd', 1, "character",  # directory which holds input subdirectories
                       'sample_dir_group_table', 's', 1, "character",  # file containing table of sample_dir, group id (integer)
                       #'sample_dir_pat',        'p', 1, "character",  # re pattern for selecting subdirectories
                       #'experiment_groups',     'e', 1, "character",  # character string of 0s and 1s indicating group membership of each column
                       'output_dir',             'O', 1, "character",  # where output file(s) go
                       'output_csvfile',            'o', 1, "character",  # output differential gene expression CSV file name
                       'help',                   'h', 0, "logical"
                      ), 
                    byrow=TRUE, ncol=4 );
opt = getopt( option_tab )

if ( ! is.null( opt$help ) ) {
    dmesg( getopt( option_tab, usage=TRUE ) );
    q( status = 1 );
}
#dmesg( "ls -lR" )
#system( paste( "ls -lR", opt$input_dir ) ) 

dmesg( "Here is the option_tab matrix of arguments" )
print( option_tab )


dmesg( "about to get samples_groups table" )
smg <- get_and_check_samples_groups( opt$sample_dir_group_table )
print( smg )

dmesg( "creating phenotype data frame for ballgown() constructor" )
pg <- data.frame( id    = basename( as.character( smg$sample_dir) ),
                  group = smg$group )
print( pg )

dmesg( "about to create ballgown" )
#bg <- ballgown( dataDir = opt$input_dir, samplePattern = opt$sample_dir_pat, meas="all" )
bg <- ballgown( samples = as.character( smg$sample_dir ),
                pData = pg,           
                verbose = TRUE,
                meas = "all" 
              )

dmesg( "here is bg str" )
print( str( bg ) )

dmesg( "here are sampleNames(bg)" )
print( sampleNames(bg) )

dmesg( "here is pData(bg)" )
print( pData(bg) )


## convert string "111000" integer into vector c(1,1,1,0,0,0)
#experiment_group_vector <- as.integer( strsplit( opt$experiment_groups, "")[[1]] )   
#dmesg( "here is experiment_group_vector" )
#print( experiment_group_vector )
# dmesg( "about to pData")
# map group membership to column names
#pData(bg) <- data.frame( id = sampleNames(bg), group = experiment_group_vector )

dmesg( "about to stattest" )

# create gene-level differential expression table
#gene_diff_ex_tab <- stattest( bg, feature='gene', meas='FPKM', covariate='group', getFC=TRUE )
gene_diff_ex_tab <- stattest( bg, feature='gene', meas='FPKM', covariate='group', getFC=TRUE )

dmesg( "about to write.table" )

write.table( gene_diff_ex_tab[,-1], 
             file=paste( opt$output_dir, opt$output_csvfile, sep="/"), 
             sep="\t", row.names=F )  # remove first column before saving

dmesg( "kthxbye" )

q( save="no" )     # don't bother saving the workspace

