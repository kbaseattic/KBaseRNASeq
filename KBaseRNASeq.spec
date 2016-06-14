#include <KBaseAssembly.spec>
#include <KBaseExpression.spec>
#include <MAK.spec>

 module KBaseRNASeq{

   /* Importing datatype objects from other modules */
   
   /*
      Create an analysis id RNASeq analysis object
      @id KBaseRNASeq.RNASeqSampleSet 
   */
      
      typedef string ws_rnaseq_sampleset_id;
   
   /*
      Id for KBaseRNASeq.RNASeqDifferentialExpression
      @id ws KBaseRNASeq.RNASeqDifferentialExpression
   */

        typedef string ws_cuffdiff_diff_exp_id;

   /* 
      Id for KBaseAssembly.SingleEndLibrary
      @name ws KBaseAssembly.SingleEndLibrary
   */
 
      	typedef string ws_SingleEndLibrary;
   /*
      Id for KBaseAssembly.PairedEndLibrary
      @name ws KBaseAssembly.PairedEndLibrary
   */

        typedef string ws_PairedEndLibrary;

   /*
      reference genome annotation id for mapping the RNA-Seq fastq file
      @id ws KBaseGenomeAnnotations.GenomeAnnotation
   */

        typedef string ws_genome_annotation_id;
   
   /* 
     Id for the handle object
     @id handle
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
      @optional genome_scientific_name
      @metadata ws handle.file_name
      @metadata ws handle.type
      @metadata ws genome_scientific_name
      @metadata ws genome_id
   */

   	typedef structure {
       		Handle handle;
                int size;	
       		ws_genome_annotation_id genome_id;
		string genome_scientific_name;
   	} GFFAnnotation;

  /*
      Id for KBaseRNASeq.GFFAnnotation
      @id ws KBaseRNASeq.GFFAnnotation
   */

        typedef string ws_referenceAnnotation_id;

  /*
      @optional genome_scientific_name handle ftp_url
      @metadata ws handle.file_name
      @metadata ws handle.type
      @metadata ws genome_id
      @metadata ws genome_scientific_name
   */

        typedef structure {
                Handle handle;
		int size;
                ws_genome_annotation_id genome_id;
		string ftp_url;
		string genome_scientific_name;
        }Bowtie2Indexes;

  /*
      Id for KBaseRNASeq.Bowtie2Indexes
      @id ws KBaseRNASeq.Bowtie2Indexes
   */

        typedef string ws_bowtieIndex_id;

/*
  Object to Describe the RNASeq SampleSet
  @optional platform num_replicates source publication_Id external_source_date sample_ids
  @metadata ws sampleset_id
  @metadata ws platform
  @metadata ws num_samples
  @metadata ws num_replicates
  @metadata ws length(condition)
*/

 typedef structure {
        string sampleset_id;
        string sampleset_desc;
	string domain;
        string platform;
        int num_samples;
        int num_replicates;
        list<string> sample_ids;
        list<string> condition;
        string source;
        string Library_type;
        string publication_Id;
        string external_source_date;
        }RNASeqSampleSet;
/*
      Id for KBaseRNASeq.RNASeqSampleSet
      @id ws KBaseRNASeq.RNASeqSampleSet
   */

        typedef string ws_Sampleset_id;



  /* Specification for the RNASeqFastq Metadata
   
    Object for the RNASeq Metadata
    @optional library_type replicate_id platform custom*/
   	
    typedef structure {
		string sample_id;
       		string library_type;
		string replicate_id;
       		string platform;
       		string condition;
		string custom;
    	}RNASeqSampleMetaData;

/*
     RNASeq fastq  object
     @optional  singleend_sample pairedend_sample metadata 
     @metadata ws singleend_sample.handle.file_name
     @metadata ws pairedend_sample.handle_1.file_name
     @metadata ws pairedend_sample.handle_2.file_name
     @metadata ws metadata.replicate_id
     @metadata ws metadata.library_type
     @metadata ws metadata.platform
     @metadata ws metadata.condition
*/

     typedef structure {
	 KBaseAssembly.SingleEndLibrary singleend_sample;
	 KBaseAssembly.PairedEndLibrary pairedend_sample;
	 ws_rnaseq_sampleset_id  sampleset_id;
	 RNASeqSampleMetaData metadata;  
     }RNASeqSample;

/*
  The workspace id of a RNASeqSample
  @id ws KBaseRNASeq.RNASeqSample   
*/

    typedef string ws_rnaseqSample_id;

/* Structure Read_mapping_sections
   @optional five_UTR three_UTR exons TSS TES introns intergenic_regions
*/

   typedef structure{
        float five_UTR;
        float three_UTR;
        float TSS;
        float TES;
        float exons;
        float introns;
        float intergenic_regions;
        }Read_mapping_sections;

/*
    Object - getAlignmentStats method
    @optional singletons multiple_alignments properly_paired alignment_rate unmapped_reads mapped_sections total_reads mapped_reads
*/
    typedef structure{
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
    Object for the RNASeq Alignment bam file
    @optional aligner_opts aligner_version aligned_using replicate_id platform size mapped_sample_id sampleset_id alignment_stats 
    @metadata ws aligned_using
    @metadata ws aligner_version
    @metadata ws genome_id
    @metadata ws size
    @metadata ws alignment_stats.total_reads
    @metadata ws alignment_stats.mapped_reads
    @metadata ws alignment_stats.alignment_rate
    @metadata ws read_sample_id
    @metadata ws library_type
    @metadata ws replicate_id
    @metadata ws condition
    @metadata ws platform
    */

    typedef structure{
	string aligned_using;
	string aligner_version;
	string library_type;
	string read_sample_id;
	string replicate_id;
	string condition;
	string platform;
        ws_genome_annotation_id genome_id;
        ws_bowtieIndex_id bowtie2_index;
	list<mapping<string opt_name, string opt_value>> aligner_opts;
        mapping<string condition,mapping<string sample_id , string replicate_id>> mapped_sample_id;
	ws_Sampleset_id sampleset_id;
	Handle file;
	int size;
        AlignmentStatsResults alignment_stats;	
    }RNASeqAlignment;


/* 
  The workspace id for a RNASeqAlignment object
  @id ws KBaseRNASeq.RNASeqAlignment
*/

 typedef string ws_samplealignment_id;

/*
  Set object for RNASeqAlignment objects
  @optional condition sample_alignments
*/

  typedef structure {
	ws_Sampleset_id sampleset_id;
	ws_genome_annotation_id genome_id;
	ws_bowtieIndex_id bowtie2_index;
	list<string> read_sample_ids;
	list<string> condition;
	list<ws_samplealignment_id> sample_alignments;
        list<mapping<string read_sample_name , string  alignment_name>> mapped_rnaseq_alignments;
        list<mapping<string read_sample_id , ws_samplealignment_id alignment_id>> mapped_alignments_ids;
	}RNASeqAlignmentSet;

/*
  The workspace id for a RNASeqAlignmentSet object
  @id ws KBaseRNASeq.RNASeqAlignmentSet
*/

 typedef string ws_alignmentSet_id;


/*
  The workspace object for a RNASeqExpression
  @optional description data_quality_level original_median external_source_date source file processing_comments mapped_sample_id
  @metadata ws type
  @metadata ws numerical_interpretation
  @metadata ws description
  @metadata ws genome_id
  @metadata ws platform
*/
  	
   typedef structure {
        string id;
        string type;
        string numerical_interpretation;
        string description;
        int data_quality_level;
        float original_median;
        string external_source_date;
        list<mapping<string feature_id,float feature_value>> expression_levels; 
        ws_genome_annotation_id genome_id; 
        ws_referenceAnnotation_id annotation_id;
	string condition;
	mapping<string sample_id,ws_samplealignment_id alignment_id> mapped_rnaseq_alignment;
        mapping<string condition,mapping<string sample_id , string replicate_id>> mapped_sample_id;
        string  platform; 
        string source; 
        Handle file;
        string processing_comments;
        string tool_used;
        string tool_version;
	list<mapping<string opt_name, string opt_value>> tool_opts; 
    }RNASeqExpression;

/*
      Id for expression sample
      @id ws KBaseRNASeq.RNASeqExpression

   */  
        typedef string ws_expression_sample_id;

/*
  Set object for RNASeqExpression objects
*/

  typedef structure {
	ws_alignmentSet_id alignmentSet_id;
        ws_Sampleset_id sampleset_id;
        ws_genome_annotation_id genome_id;
        list<string> sample_ids;
        list<string> condition;
        list<ws_expression_sample_id> sample_expression_ids;
        list<mapping<string read_sample_name , string expression_name>> mapped_expression_objects;
        list<mapping<string read_sample_id , ws_expression_sample_id expression_id>> mapped_expression_ids;
        }RNASeqExpressionSet;
/*
      Id for expression sample set
      @id ws KBaseRNASeq.RNASeqExpressionSet

   */
        typedef string ws_expressionSet_id;
/*
   Object RNASeqDifferentialExpression file structure
   @optional tool_opts tool_version sample_ids comments
*/
   typedef structure {
	string tool_used;
	string tool_version;
        list<mapping<string opt_name, string opt_value>> tool_opts;	
      	Handle file;
        list<string> sample_ids;
	list<string> condition;
	ws_expressionSet_id expressionSet_id;
	ws_alignmentSet_id alignmentSet_id;
        ws_Sampleset_id sampleset_id;
	string comments;
      	}RNASeqDifferentialExpression;

/*
    Object for the cummerbund plot
    @optional png_json_handle plot_title plot_description
*/
    typedef structure {
       Handle png_handle;
       Handle png_json_handle;
       string plot_title;
       string plot_description;
       }cummerbundplot;
/*
  List of cummerbundplot
*/


    typedef list<cummerbundplot> cummerbundplotSet;

/*
   Object type for the cummerbund_output
*/
    typedef structure {
       cummerbundplotSet cummerbundplotSet;
       string rnaseq_experiment_id;
       string cuffdiff_input_id;
       }cummerbund_output;


/*
     Object for Report type
*/
    typedef structure {
		string report_name;
		string report_ref;
    }ResultsToReport;

/* FUNCTIONS used in the service */
 	
   typedef structure{
	string ws_id;
   	string sampleset_id;
   	string sampleset_desc;
   	string domain;
        string platform;
	list<string> sample_ids;
   	list<string> condition;
   	string source;
   	string Library_type;
	string publication_id;
   	string external_source_date;
   	}CreateRNASeqSampleSetParams;
   	
  async funcdef CreateRNASeqSampleSet(CreateRNASeqSampleSetParams params)
	returns(RNASeqSampleSet) authentication required;
	
   typedef structure{
	string ws_id;
	string reference;
	string output_obj_name;
	}Bowtie2IndexParams;

  async funcdef BuildBowtie2Index(Bowtie2IndexParams params)
     returns(ResultsToReport) authentication required;

   typedef structure{
        string ws_id;
        ws_genome_annotation_id reference;
        string output_obj_name;
        }GetFeaturesToGTFParams;

  async funcdef GetFeaturesToGTF(GetFeaturesToGTFParams params)
     returns(ResultsToReport) authentication required;   
   
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
        string sampleset_id;
	string genome_id;
	string bowtie_index;
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
     returns(ResultsToReport) authentication required;

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
     string read_sample; 
     string genome_id;
     string bowtie2_index;
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
     returns (ResultsToReport) authentication required;

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

typedef structure{
	string ws_id;
	string sample_alignment;
	int num_threads;
        /*string library-type; */
        /*string library-norm-method; */
	int min-intron-length;
	int max-intron-length;
	int overhang-tolerance;
        }CufflinksParams;

  async funcdef CufflinksCall(CufflinksParams params)
    returns (ws_expression_sample_id) authentication required;

	typedef structure{
	string ws_id;
        RNASeqSampleSet rnaseq_exp_details;
        string output_obj_name;
        string time-series;
	string library-type;
        string library-norm-method;
	string multi-read-correct;
        int  min-alignment-count;
	string dispersion-method;
	string no-js-tests;
	int frag-len-mean;
	int frag-len-std-dev;
	int max-mle-iterations;
	string compatible-hits-norm;
	string no-length-correction;
        }CuffdiffParams;

  async funcdef CuffdiffCall(CuffdiffParams params)
   returns (RNASeqDifferentialExpression) authentication required;
};
