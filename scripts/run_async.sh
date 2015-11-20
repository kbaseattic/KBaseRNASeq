script_dir=$(dirname "$(readlink -f "$0")")
export KB_DEPLOYMENT_CONFIG=$script_dir/../deploy.cfg
WORK_DIR=/kb/module/work
echo $script_dir
cat $WORK_DIR/token | xargs sh $script_dir/../bin/run_KBaseRNASeq.sh $WORK_DIR/input.json $WORK_DIR/output.json
