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
      Id for KBaseFile.SingleEndLibrary
      @id ws KBaseFile.SingleEndLibrary
   */
 
      	typedef string ws_singleEndLibrary_id;
   /*
      Id for KBaseGenomes.ContigSet
      @id ws KBaseGenomes.ContigSet
   */

        typedef string ws_reference_assembly_id;
 
   /* 
      Id for KBaseFile.PairedEndLibrary
      @id ws KBaseFile.PairedEndLibrary
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
      @optional genome_id
      @metadata ws handle.file_name
      @metadata ws handle.type
      @metadata ws genome_scientific_name
      @metadata ws genome_id
   */

   	typedef structure {
       		Handle handle;
       		ws_genome_id genome_id;
		string genome_scientific_name;
   	} ReferenceAnnotation;

  /*
      Id for ReferenceAnnotation
      @id ws ReferenceAnnotation
   */

        typedef string ws_referenceAnnotation_id;

  /*
      @optional genome_id
      @metadata ws handle.file_name
      @metadata ws handle.type
      @metadata ws genome_id
      @metadata ws genome_scientific_name
   */

        typedef structure {
                Handle handle;
                ws_genome_id genome_id;
		string genome_scientific_name;
        }BowtieIndexes;

  /*
      Id for BowtieIndexes
      @id ws BowtieIndexes
   */

        typedef string ws_bowtieIndex_id;

  /* Data structure for top level information for sample annotation and ontology */
 
    /* Kbase SampleAnnotation ID */ 
    typedef string sample_annotation_id; 

    /* Kbase OntologyID  */ 
    typedef string ontology_id; 
 
    /* Kbase OntologyName  */ 
    typedef string ontology_name; 

    /* Kbase Ontology Description  */ 
    typedef string ontology_definition; 
    
    typedef structure {
	sample_annotation_id sample_annotation_id;
	ontology_id ontology_id;
	ontology_name ontology_name;
	ontology_definition ontology_definition;
    } SampleAnnotation;
 

  /* List of Sample Annotations */

   typedef list<SampleAnnotation> sample_annotations;

  /* Specification for the RNASeqFastq Metadata
   
    Object for the RNASeq Metadata
    @optional platform source tissue condition source_id ext_source_date sample_desc title sample_annotations  */
   	
    typedef structure {
		string sample_id;
       		string library_type;
		string replicate_id;
       		string platform;
       		string sample_desc;
       		string title;
       		string source;
       		string source_id;
       		string ext_source_date;
       		string domain;
       		ws_genome_id genome_id;
		sample_annotations sample_annotations;
       		list<string> tissue;
       		list<string> condition;
    	}RNASeqSampleMetaData;

/*
      Complete List of RNASeq MetaData
    
*/
   	 typedef list<RNASeqSampleMetaData> RNASeqSamplesMetaData;
 

/*
     RNASeq fastq  object
     @optional  singleend_sample pairedend_sample metadata
     @metadata ws analysis_id
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
	 ws_singleEndLibrary_id singleend_sample;
	 ws_pairedEndLibrary_id pairedend_sample;
	 ws_rnaseq_analysis_id  analysis_id;
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
	Handle handle1;
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
  @optional experiment_desc num_replicates title platform genome_id source tissue condition annotation_id publication_id source_id external_source_date sample_ids alignments expression_values sample_annotations_map
  @metadata ws experiment_id
  @metadata ws title
  @metadata ws experiment_desc
  @metadata ws platform
  @metadata ws genome_id
  @metadata ws num_samples
  @metadata ws num_replicates 
  @metadata ws annotation_id 
*/

 typedef structure {
        string experiment_id;
        string title;
	string experiment_desc;
        string platform;
        ws_genome_id genome_id;
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
        string source_id;
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
	Handle handle1;
	RNASeqAnalysis analysis;
	}RNASeqCuffmergetranscriptome;

/*
   Object RNASeqDifferentialExpression file structure
*/
   typedef structure {
      	Handle handle1;
	RNASeqAnalysis analysis;
      	}RNASeqCuffdiffdifferentialExpression;
	
/* FUNCTIONS used in the service */

/* Function parameters to call tophat */

   typedef structure{
       ws_rnaseqSample_id sample_id;
       list<string> results;
       }fastqcParams;

funcdef CallFastqc(fastqcParams params)
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
     ws_rnaseqSample_id sample_id;
     string output_obj_name;
     ws_reference_assembly_id reference;
     ws_bowtieIndex_id bowtie_index;
     t_opts_str opts_dict;
     ws_referenceAnnotation_id annotation_gtf;
     }TophatParams;

funcdef TophatCall(TophatParams params)
     returns (string job_id) authentication required;

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
        ws_samplealignment_id alignment_sample_id;
        string output_obj_name;
        ws_referenceAnnotation_id annotation_gtf;
        cuff_opts opts_dict;
        }CufflinksParams;

funcdef CufflinksCall(CufflinksParams params)
    returns (string job_id) authentication required;

typedef structure{
        RNASeqAnalysis analysis;
        string output_obj_name;
        mapping <string Cuffmerge_opts, int num_threads> opts_dict;
        }CuffmergeParams;

 
funcdef CuffmergeCall(CuffmergeParams params)
    returns (string job_id) authentication required;


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
        RNASeqAnalysis rnaseq_exp_details;
        string output_obj_name;
        ws_referenceAnnotation_id annotation_gtf;
        cuffdiff_opts  opts_dict;
        }CuffdiffParams;

funcdef CuffdiffCall(CuffdiffParams params)
   returns (string job_id) authentication required;


typedef structure{
        ws_samplealignment_id alignment_sample_id;
        }AlignmentStatsParams;

funcdef getAlignmentStats(AlignmentStatsParams params)
   returns (string job_id) authentication required;


typedef structure{
        ws_expression_sample_id expression_sample;
        int number_of_bins;
        }ExpressionHistogramParams;
        
	
funcdef createExpressionHistogram(ExpressionHistogramParams params)
   returns (string job_id) authentication required;


typedef structure{
        RNASeqAnalysis analysis;
	string out_obj_name;
        }ExpressionSeriesParams;

funcdef createExpressionSeries(ExpressionSeriesParams params)
   returns (string job_id) authentication required;

typedef structure{
        /* RNASeqAnalysis rnaseq_exp_details;*/
	string in_id; /* the input pointer for rnaseq_exp_detail object? */
	string ws_id; /* the working workspace name */
	string out_id; /* final expression matrix name */
        }ExpressionMatrixParams;

funcdef createExpressionMatrix(ExpressionMatrixParams params)
   returns (string job_id) authentication required;


};
