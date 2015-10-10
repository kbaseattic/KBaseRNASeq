#BEGIN_HEADER
#END_HEADER


class KBaseRNASeq:
    '''
    Module Name:
    KBaseRNASeq

    Module Description:
    
    '''

    ######## WARNING FOR GEVENT USERS #######
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    #########################################
    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        #END_CONSTRUCTOR
        pass

    def CallFastqc(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CallFastqc
        #END CallFastqc

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CallFastqc return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def BuildBowtie2Index(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN BuildBowtie2Index
        #END BuildBowtie2Index

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method BuildBowtie2Index return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CallBowtie2(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CallBowtie2
        #END CallBowtie2

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CallBowtie2 return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def TophatCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN TophatCall
        #END TophatCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method TophatCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CufflinksCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CufflinksCall
        #END CufflinksCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CufflinksCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CuffmergeCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CuffmergeCall
        #END CuffmergeCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CuffmergeCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CuffdiffCall(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CuffdiffCall
        #END CuffdiffCall

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CuffdiffCall return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def getAlignmentStats(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN getAlignmentStats
        #END getAlignmentStats

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method getAlignmentStats return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionHistogram(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionHistogram
        #END createExpressionHistogram

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionHistogram return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def CallCummeRbund(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN CallCummeRbund
        #END CallCummeRbund

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method CallCummeRbund return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionSeries(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionSeries
        #END createExpressionSeries

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionSeries return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]

    def createExpressionMatrix(self, ctx, params):
        # ctx is the context object
        # return variables are: job_id
        #BEGIN createExpressionMatrix
        #END createExpressionMatrix

        # At some point might do deeper type checking...
        if not isinstance(job_id, basestring):
            raise ValueError('Method createExpressionMatrix return value ' +
                             'job_id is not type basestring as required.')
        # return the results
        return [job_id]
