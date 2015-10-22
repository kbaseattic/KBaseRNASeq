script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
WORK_DIR=/kb/module/work
cat $WORK_DIR/token | xargs sh /kb/deployment/bin/run_KBaseRNASeq.sh $WORK_DIR/input.json $WORK_DIR/output.json
