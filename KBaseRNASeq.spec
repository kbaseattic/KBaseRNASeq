#include <KBaseAssembly.spec>

 module KBaseRNASeq{

   /* Importing datatype objects from other modules */
   /* indicates true or false values, false <= 0, true >=1 */
	typedef int bool; 
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
       		/*ws_genome_annotation_id genome_id;*/
       		string genome_id;
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
                /*ws_genome_annotation_id genome_id;*/
                string genome_id;
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
    @optional aligner_opts aligner_version aligned_using replicate_id platform size mapped_sample_id sampleset_id alignment_stats bowtie2_index 
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
        /*ws_genome_annotation_id genome_id;*/
        string genome_id;
        ws_bowtieIndex_id bowtie2_index;
	mapping<string opt_name, string opt_value> aligner_opts;
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
  @optional condition sample_alignments bowtie2_index aligned_using aligner_version aligner_opts
  @metadata ws aligned_using
  @metadata ws aligner_version
  @metadata ws sampleset_id
*/

  typedef structure {
	string aligned_using;
	string aligner_version;
	mapping<string opt_name, string opt_value> aligner_opts;
	ws_Sampleset_id sampleset_id;
	/*ws_genome_annotation_id genome_id;*/
	string genome_id;
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
  @optional description platform tool_opts data_quality_level original_median external_source_date source file processing_comments mapped_sample_id tpm_expression_levels
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
        mapping<string feature_id,float feature_value> expression_levels; 
        mapping<string feature_id,float feature_value> tpm_expression_levels;
        /*ws_genome_annotation_id genome_id;*/
	string genome_id;
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
	mapping<string opt_name, string opt_value> tool_opts; 
    }RNASeqExpression;

/*
      Id for expression sample
      @id ws KBaseRNASeq.RNASeqExpression

   */  
        typedef string ws_expression_sample_id;

/*
  Set object for RNASeqExpression objects
  @optional sample_ids condition tool_used tool_version tool_opts
  @metadata ws tool_used
  @metadata ws tool_version
  @metadata ws alignmentSet_id 
*/

  typedef structure {
        string tool_used;
        string tool_version;
	mapping<string opt_name, string opt_value> tool_opts;
	ws_alignmentSet_id alignmentSet_id;
        ws_Sampleset_id sampleset_id;
        /*ws_genome_annotation_id genome_id;*/
	string genome_id;
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
        /*ws_genome_annotation_id genome_id;*/
	string genome_id;
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


    typedef structure {
               string function;
               string gene;
               string locus;
               float log2fc;
               float log2fc_f;
               float log2fc_fa;
               float p_value;
               float p_value_f;
               string significant;
               float value_1;
               float value_2;
                } gene_expression_stat;
 
    typedef list <gene_expression_stat> voldata;
    typedef structure {
               string condition_1;
               string condition_2;
               voldata voldata;
                } condition_pair_unit;
 
 
    typedef structure {
               list <condition_pair_unit> condition_pairs;
               list <string> unique_conditions;
                } DifferentialExpressionStat;
 
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
        string ws_id;
        string sampleset_id;
        string genome_id;
	int num_threads;
        string quality_score;
        int skip;
        int trim3;
        int trim5;
	int np;
	int minins;
	int maxins;
	string orientation;
        int min_intron_length;
	int max_intron_length;
	bool no_spliced_alignment;
	bool transcriptome_mapping_only;
        string tailor_alignments;
        }Hisat2Params;

  async funcdef Hisat2Call(Hisat2Params params)
     returns(ResultsToReport) authentication required;

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
        string ws_id;
        string sample_alignment;
        int num-threads;
	string label;
        float min_isoform_abundance;
	int a_juncs;
	int  min_length;
	float j_min_reads;
	float c_min_read_coverage;
	int gap_sep_value;
	bool disable_trimming;
	bool ballgown_mode;
	bool skip_reads_with_no_ref;
	string merge;
        }StringTieParams;

async funcdef StringTieCall(StringTieParams params)
    returns (ResultsToReport) authentication required;

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
    returns (ResultsToReport) authentication required;

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

        typedef structure{
        string ws_id;
        RNASeqExpressionSet expressionset_id;
        string output_obj_name;
        int num_threads;
        }DifferentialExpParams;

  async funcdef DiffExpCallforBallgown(DifferentialExpParams params)
   returns (RNASeqDifferentialExpression) authentication required;

};
