#include <KBaseAssembly.spec>
#include <KBaseExpression.spec>

 module KBaseRNASeq{

   /* Importing datatype objects from other modules */
   
   /*
      Create an analysis id RNASeq analysis object
      @id KBaseRNASeq.RNASeqAnalysis 
   */
      
      typedef string ws_rnaseq_analysis_id;

   /*
      reference genome id for mapping the RNA-Seq fastq file
      @id ws KBaseGenomes.Genome
   */
      
	typedef string ws_genome_id;
 
   /*
      Id for expression sample
      @id ws KBaseExpression.ExpressionSample

   */
        typedef string ws_expression_sample_id;

   /* 
    	List of Expression sample ids

   */
	
   	typedef list<ws_expression_sample_id> ws_expression_sample_ids;
  
   /* 
      Id for KBaseAssembly.SingleEndLibrary
      @id ws KBaseAssembly.SingleEndLibrary
   */
 
      	typedef string ws_singleEndLibrary_id;
   /*
      Id for KBaseGenomes.ContigSet
      @id ws KBaseGenomes.ContigSet
   */

        typedef string ws_reference_assembly_id;
 
   /* 
      Id for KBaseAssembly.PairedEndLibrary
      @id ws KBaseAssembly.PairedEndLibrary
   */
 
      	typedef string ws_pairedEndLibrary_id;
   
   /* 
     Id for the handle object
   */
	typedef string HandleId;
 
   /*
      @optional hid file_name type url remote_md5 remote_sha1
   */
	
   	typedef structure {
       		HandleId hid;
       		string file_name;
       		string id;
       		string type;
       		string url;
       		string remote_md5;
       		string remote_sha1;
   	} Handle;


  /*
      @optional genome_id genome_scientific_name
      @metadata ws handle.file_name
      @metadata ws handle.type
      @metadata ws genome_scientific_name
      @metadata ws genome_id
   */

   	typedef structure {
       		Handle handle;
                int size;	
       		ws_genome_id genome_id;
		string genome_scientific_name;
   	} ReferenceAnnotation;

  /*
      Id for KBaseRNASeq.ReferenceAnnotation
      @id ws KBaseRNASeq.ReferenceAnnotation
   */

        typedef string ws_referenceAnnotation_id;

  /*
      @optional genome_id genome_scientific_name handle ftp_url
      @metadata ws handle.file_name
      @metadata ws handle.type
      @metadata ws genome_id
      @metadata ws genome_scientific_name
   */

        typedef structure {
                Handle handle;
		int size;
                ws_genome_id genome_id;
		string ftp_url;
		string genome_scientific_name;
        }Bowtie2Indexes;

  /*
      Id for KBaseRNASeq.Bowtie2Indexes
      @id ws KBaseRNASeq.Bowtie2Indexes
   */

        typedef string ws_bowtieIndex_id;

  /* Kbase SampleAnnotation ID */
    typedef string sample_annotation_id;

  /* Kbase OntologyID  */
    typedef string ontology_id;

  /* list of Kbase Ontology IDs */
    typedef list<ontology_id> ontology_ids;

  /* Kbase OntologyName */
    typedef string ontology_name;

  /* Kbase OntologyDefinition */
    typedef string ontology_definition;

  /* Data structure for top level information for sample annotation and ontology */
    typedef structure {
        sample_annotation_id sample_annotation_id;
        ontology_id ontology_id;
        ontology_name ontology_name;
        ontology_definition ontology_definition;
    } SampleAnnotation;

  /*
	List of KBaseExpression.SampleAnnotation
  */
	typedef list<SampleAnnotation> sample_annotations;

  /* Specification for the RNASeqFastq Metadata
   
    Object for the RNASeq Metadata
    @optional library_type platform source tissue condition source_id ext_source_date sample_desc title sample_annotations genome_id genome_scientific_name custom*/
   	
    typedef structure {
		string sample_id;
       		string library_type;
		string replicate_id;
       		string platform;
       		string sample_desc;
       		string title;
       		string source;
       		string  source_id;
       		string ext_source_date;
       		string domain; 
       		ws_genome_id genome_id;
		string genome_scientific_name;
		sample_annotations sample_annotations;
       		string tissue;
       		string condition;
		string custom;
    	}RNASeqSampleMetaData;

/*
      Complete List of RNASeq MetaData
    
*/
     typedef list<RNASeqSampleMetaData> RNASeqSamplesMetaData;
 

/*
     RNASeq fastq  object
     @optional  singleend_sample pairedend_sample metadata analysis_id analysis_desc
     @metadata ws analysis_id
     @metadata ws analysis_desc
     @metadata ws metadata.sample_id
     @metadata ws metadata.replicate_id
     @metadata ws metadata.library_type
     @metadata ws metadata.platform
     @metadata ws metadata.title
     @metadata ws metadata.source
     @metadata ws metadata.source_id
     @metadata ws metadata.sample_desc
     @metadata ws length(metadata.tissue)
     @metadata ws length(metadata.condition)
     @metadata ws metadata.genome_id
     @metadata ws metadata.ext_source_date
     @metadata ws length(metadata.sample_annotations)
*/

     typedef structure {
	 KBaseAssembly.SingleEndLibrary singleend_sample;
	 KBaseAssembly.PairedEndLibrary pairedend_sample;
	 ws_rnaseq_analysis_id  analysis_id;
	 string analysis_desc;
	 RNASeqSampleMetaData metadata;  
     }RNASeqSample;

  
/* 
    list of RNASeqSamples
*/
    
    typedef list<RNASeqSample> RNASeqSamplesSet;

/*
  The workspace id of a RNASeqSample
  @id ws KBaseRNASeq.RNASeqSample   
*/

    typedef string ws_rnaseqSample_id;

    
 /*
    Object for the RNASeq Alignment bam file
    @optional aligner_opts aligner_version aligned_using metadata
    @metadata ws metadata.sample_id
    @metadata ws metadata.replicate_id
    @metadata ws metadata.platform
    @metadata ws metadata.title
    @metadata ws metadata.source
    @metadata ws metadata.source_id
    @metadata ws metadata.ext_source_date
    @metadata ws metadata.sample_desc
    @metadata ws metadata.genome_id
    @metadata ws length(metadata.tissue)
    @metadata ws length(metadata.condition)
    @metadata ws aligned_using
    @metadata ws aligner_version

    */

    typedef structure{
	string aligned_using;
	string aligner_version;
	list<mapping<string opt_name, string opt_value>> aligner_opts;
	Handle file;
	int size;
	RNASeqSampleMetaData metadata;
    }RNASeqSampleAlignment;

/*
      list of RNASeqSampleAlignment
*/

 typedef list<RNASeqSampleAlignment> RNASeqSampleAlignmentSet;


/* 
  The workspace id for a RNASeqSampleAlignment object
  @id ws KBaseRNASeq.RNASeqSampleAlignment
*/

 typedef string ws_samplealignment_id;

/* Structure Read_mapping_sections
   @optional introns exons splice_junctions intergenic_regions
*/

   typedef structure{
        int introns;
        int exons;
        int splice_junctions;
        int intergenic_regions;
        }Read_mapping_sections;
	
/*
    Object - getAlignmentStats method
    @optional singletons multiple_alignments properly_paired alignment_rate unmapped_reads mapped_sections total_reads mapped_reads
*/
    typedef structure{
        /* later change this to ws_samplealignment_id */
	string alignment_id;
        int properly_paired;
        int multiple_alignments;
        int singletons;
        float alignment_rate;
        Read_mapping_sections mapped_sections;
        int unmapped_reads;
        int mapped_reads;
        int total_reads;
        }AlignmentStatsResults;

/*
  Object type to define replicate group
/*

  typedef structure{
	ws_rnaseq_analysis_id analysis_id;
        string sample_name;
        ws_rnaseqSample_ids sample_replicate_group;
    }RNASeqSampleReplicateGroup;

/*
  Object type to define id to RNASeqSampleReplicateGroup
  @id ws KBaseRNASeq.RNASeqSampleReplicateGroup
*/
   
  typedef string ws_RNASeqSampleReplicateGroup_id;

/*
  The mapping object for the ws_rnaseqSample_id with the ws_samplealignment_id
*/

  typedef mapping<ws_rnaseqSample_id sampleID,ws_samplealignment_id sampleAlignmentID> mapped_sample_alignment;

/*
  The mapping object for the ws_rnaseqSample_id with the ws_expression_sample_id
*/

  typedef mapping<ws_rnaseqSample_id sampleID,ws_expression_sample_id exp_sampleID> mapped_sample_expression;


/*
  Object to Describe the RNASeq analysis
  @optional experiment_desc num_replicates title platform genome_id source tissue condition annotation_id publication_id source_ids external_source_date sample_ids alignments expression_values sample_annotations_map genome_scientific_name sample_rep_groups alignments expression_values transcriptome_id cuffdiff_diff_exp_id experiment_design
  @metadata ws experiment_id
  @metadata ws title
  @metadata ws experiment_desc
  @metadata ws experiment_design
  @metadata ws platform
  @metadata ws genome_id
  @metadata ws genome_scientific_name
  @metadata ws num_samples
  @metadata ws num_replicates 
  @metadata ws annotation_id 
*/

 typedef structure {
        string experiment_id;
        string title;
	string experiment_desc;
	string experiment_design;
	string domain;
        string platform;
        ws_genome_id genome_id;
	string genome_scientific_name;
        int num_samples;
        int num_replicates;
        list<ws_rnaseqSample_id> sample_ids;
	list<ws_RNASeqSampleReplicateGroup_id>  sample_rep_groups;
	list<mapped_sample_alignment> alignments;
	list<mapped_sample_expression> expression_values;
	string transcriptome_id;
	string cuffdiff_diff_exp_id;
        list<string> tissue;
        list<string> condition;
	list<mapping<string sample_name,sample_annotations>> sample_annotations_map;
        ws_referenceAnnotation_id annotation_id;
        string source;
        string Library_type;
        string publication_id;
        list<string> source_ids;
        string external_source_date;
        }RNASeqAnalysis;


/*
    Object for RNASeq Merged transcriptome assembly
 	@metadata ws analysis.experiment_id
	@metadata ws analysis.title
	@metadata ws analysis.Library_type
	@metadata ws analysis.platform
	@metadata ws analysis.num_samples
	@metadata ws analysis.num_replicates
	@metadata ws length(analysis.sample_ids)
	@metadata ws length(analysis.tissue)
	@metadata ws length(analysis.condition)
/*

   typedef structure{
	Handle file;
	RNASeqAnalysis analysis;
	}RNASeqCuffmergetranscriptome;

/*
   Object RNASeqDifferentialExpression file structure
*/
   typedef structure {
      	Handle file;
	RNASeqAnalysis analysis;
      	}RNASeqCuffdiffdifferentialExpression;
	
/* FUNCTIONS used in the service */

/* Function parameters to call fastqc */

   typedef structure{
       ws_rnaseqSample_id sample_id;
       string output_obj_name;
       list<mapping<string parameter ,string values>> results;
       }fastqcParams;

async funcdef fastqcCall(fastqcParams params)
     returns(string job_id) authentication required;
	
   typedef structure{
	string ws_id;
	ws_singleEndLibrary_id singleend_sample;
	ws_pairedEndLibrary_id pairedend_sample;
	ws_rnaseq_analysis_id analysis_id;
	string sample_id;
	string library_type;
	string replicate_id;
	string platform;
	string sample_desc;
	string  title;
	string source;
	string ext_source_date;
	string domain;
	ws_genome_id genome_id;
	list<string> tissue;
	list<string> condition;
   }associateReadsParams;
	
async funcdef associateReads(associateReadsParams params)
     returns(UnspecifiedObject) authentication required;
 	
   typedef structure{
	string ws_id;
   	string experiment_id;
   	string title;
   	string experiment_desc;
	string experiment_design;
	string domain;
   	string platform;
   	ws_genome_id genome_id;
   	int num_samples;
   	int num_replicates;
	list<ws_rnaseqSample_id> sample_ids;
   	list<string> tissue;
   	list<string> condition;
   	ws_referenceAnnotation_id annotation_id;
   	string source;
   	string Library_type;
   	string publication_id;
   	string external_source_date;
   	}SetupRNASeqAnalysisParams;
   	
async funcdef SetupRNASeqAnalysis(SetupRNASeqAnalysisParams params)
	returns(UnspecifiedObject) authentication required;
	
   typedef structure{
	string ws_id;
	ws_reference_assembly_id reference;
	string output_obj_name;
	}Bowtie2IndexParams;

async funcdef BuildBowtie2Index(Bowtie2IndexParams params)
     returns(UnspecifiedObject) authentication required;
   
   typedef structure{
	int skip;
	int upto;
	int trim5;
	int trim3;
	string phred33;
	string phred64;
	string int-quals;
	string local;
	string end-to-end;
	string very-fast;
	string fast;
	string very-sensitive;
	string sensitive;
	string very-fast-local;
	string very-sensitive-local;
	string fast-local;
	string fast-sensitive;
	int N;
	int L;
	int dpad;
	int gbar;
	int ma;
	int mp;
	int np;
	int threads;
	int offrate;
	string qc-filter;
	string seed;
	string non-deterministic;
	}b_opts;

   typedef mapping<string Bowtie2_opts,b_opts opts_bowtie2> b_opts_str;
	
   typedef structure{
	string ws_id;
	ws_rnaseqSample_id sample_id;
        ws_bowtieIndex_id bowtie2_index;
	string output_obj_name;
       /* b_opts_str opts_dict;	*/
        string phred33;
        string phred64;
        string local;
        string very-fast;
        string fast;
        string very-sensitive;
        string sensitive;
        string very-fast-local;
        string very-sensitive-local;
        string fast-local;
        string fast-sensitive;
	}Bowtie2Params;

async funcdef Bowtie2Call(Bowtie2Params params) 
     returns(string job_id) authentication required;

typedef structure{
     string read-mismatches;
     string read-gap-length;
     int read-edit-dist;
     int read-realign-edit-dist;
     string bowtie1;
     string output-dir;
     int mate-inner-dist;
     int mate-std-dev;
     int min-anchor-length;
     int splice-mismatches;
     int min-intron-length;
     int max-intron-length;
     int max-insertion-length;
     int num-threads;
     int max-multihits;
     string report-secondary-alignments;
     string no-discordant;
     string no-mixed;
     string no-coverage-search;
     string coverage-search;
     string microexon-search;
     string library-type;
     int segment-mismatches;
     int segment-length;
     int min-segment-intron;
     int max-segment-intron;
     int min-coverage-intron;
     int max-coverage-intron;
     string b2-very-fast;
     string b2-fast;
     string b2-sensitive;
     string b2-very-sensitive;
     string fusion-search;
     int fusion-anchor-length;
     int fusion-min-dist;
     int fusion-read-mismatches;
     int fusion-multireads;
     int fusion-multipairs;
     string fusion-ignore-chromosomes;
     }t_opts;

typedef mapping<string Tophat_opts,t_opts opts_tophat> t_opts_str;

 typedef structure{
     string ws_id;
     ws_rnaseqSample_id sample_id;
     string output_obj_name;
     ws_reference_assembly_id reference;
     ws_bowtieIndex_id bowtie_index;
     /* t_opts_str opts_dict; */ 
     int read_mismatches;
     int read_gap_length;
     int read_edit_dist;
     int min_intron_length;
     int max_intron_length;
     int num_threads;
     string report_secondary_alignments;
     string no_coverage_search;
     string library_type; 
     ws_referenceAnnotation_id annotation_gtf;
     }TophatParams;

async funcdef TophatCall(TophatParams params)
     returns (UnspecifiedObject) authentication required;

 typedef structure{
        int num_threads;
        string library-type;
        string library-norm-method;
        int frag-len-mean;
        int frag-len-std-dev;
        string upper-quartile-norm;
        string total-hits-norm;
        string compatible-hits-norm;
        int max-mle-iterations;
        int max-bundle-frags;
        string no-effective-length-correction;
        string no-length-correction;
        float min-isoform-fraction;
        float pre-mrna-fraction;
        int max-intron-length;
        float junc-alpha;
        string small-anchor-fraction;
        int min-frags-per-transfrag;
        int overhang-tolerance;
        int max-bundle-length;
        int min-intron-length;
        int trim-3-avgcov-thresh;
        int trim-3-dropoff-frac;
        float max-multiread-fraction;
        }opts_cufflinks;

typedef mapping <string Cufflinks_opts,opts_cufflinks> cuff_opts;

typedef structure{
	string ws_id;
        ws_samplealignment_id alignment_sample_id;
        string output_obj_name;
        ws_referenceAnnotation_id annotation_gtf;
        /*cuff_opts opts_dict;*/
	int num_threads;
        /*string library-type; */
        /*string library-norm-method; */
	int min-intron-length;
	int max-intron-length;
	int overhang-tolerance;
        }CufflinksParams;

async funcdef CufflinksCall(CufflinksParams params)
    returns (UnspecifiedObject) authentication required;

typedef structure{
	string ws_id;
        RNASeqAnalysis analysis;
        string output_obj_name;
        /*mapping <string Cuffmerge_opts, int num_threads> opts_dict; *
        }CuffmergeParams;

 
async funcdef CuffmergeCall(CuffmergeParams params)
    returns (UnspecifiedObject) authentication required;

typedef structure{
        int num-threads;
        string time-series;
        string total-hits-norm;
        string compatible-hits-norm;
        string multi-read-correct;
        string min-alignment-count;
        float FDR;
        string library-type;
        string library-norm-method;
        string dispersion-method;
        int frag-len-mean;
        int frag-len-std-dev;
        int max-mle-iterations;
        string poisson-dispersion;
        }opts_cuffdiff;

typedef mapping <string diff_opts,opts_cuffdiff> cuffdiff_opts;

	typedef structure{
	string ws_id;
        RNASeqAnalysis rnaseq_exp_details;
        string output_obj_name;
        ws_referenceAnnotation_id annotation_gtf;
        /*cuffdiff_opts  opts_dict;*/
	int num-threads;
	list<string> labels;
        string time-series;
	string library-type;
        string library-norm-method;
	string multi-read-correct;
        int  min-alignment-count;
        }CuffdiffParams;

async funcdef CuffdiffCall(CuffdiffParams params)
   returns (UnspecifiedObject) authentication required;


typedef structure{
	string ws_id;
        /*ws_samplealignment_id alignment_sample_id;*/
	string alignment_sample_id;
	string output_obj_name;
        }AlignmentStatsParams;

async funcdef getAlignmentStats(AlignmentStatsParams params)
   returns (AlignmentStatsResults) authentication required;

typedef structure{
	string ws_id;
        ws_expression_sample_id expression_sample;
        int number_of_bins;
	string output_obj_name;
        }ExpressionHistogramParams;
        
async funcdef createExpressionHistogram(ExpressionHistogramParams params)
   returns (UnspecifiedObject) authentication required;

typedef structure{
	string ws_id;
        RNASeqAnalysis analysis;
	string out_obj_name;
        }ExpressionSeriesParams;

typedef structure{
	string ws_id;
	RNASeqAnalysis analysis;
	string out_obj_name;
	}CummeRbundParams;
	
async funcdef cummeRbundCall(CummeRbundParams params)
   returns (UnspecifiedObject) authentication required;

async funcdef createExpressionSeries(ExpressionSeriesParams params)
   returns (UnspecifiedObject) authentication required;

typedef structure{
	string ws_id;
        RNASeqAnalysis rnaseq_exp_details;
	string out_obj_name; /* final expression matrix name */
        }ExpressionMatrixParams;

funcdef createExpressionMatrix(ExpressionMatrixParams params)
   returns (UnspecifiedObject) authentication required;

};
