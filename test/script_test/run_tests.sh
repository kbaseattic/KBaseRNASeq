#!/bin/bash
export KB_RUNTIME=/kb/runtime
export PYTHONPATH="/mnt/kb/dev_container/modules/KBaseRNASeq/lib"
python /mnt/kb/dev_container/modules/KBaseRNASeq/test/script_test/basic_test.py $1 $2 $3
