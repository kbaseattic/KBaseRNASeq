script_dir=$(dirname "$(readlink -f "$0")")
export PATH=/kb/runtime/bin:/kb/deployment/bin:$PATH
export PYTHONPATH=/kb/deployment/lib:$PYTHONPATH
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
WORK_DIR=/kb/deployment/services/KBaseRNASeq
cd $WORK_DIR
python /kb/deployment/lib/.py --host 0.0.0.0 --port 6606
# >$WORK_DIR/py_server.out 2>$WORK_DIR/py_server.err
# echo $pid > $WORK_DIR/pid.txt
