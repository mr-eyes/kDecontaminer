cd /groups/lorolab/mr-eyes/oveview_exp
SAMPLES="Ast25B Ast26B Ast27B Ast28B Ast29B Ast30A Ast34D Ast35D Ast36C Ast42B Ast44B Ast45B AW2C AW3D AW8D"
SAMPLES_DIR="/groups/lorolab/mr-eyes/final_experiment/samples_cDBGs"
OUTPUT_DIR="/groups/lorolab/mr-eyes/oveview_exp/samples_kfs"
kmerCount="/groups/lorolab/mr-eyes/oveview_exp/kDecontaminer/kmerCount.py"
for SAMPLE in $SAMPLES; 
do 
    CMD="python $kmerCount $SAMPLES_DIR/$SAMPLE/cDBG_k75_$SAMPLE.unitigs.fa $OUTPUT_DIR";
    clusterize -d -n 2 "${CMD}"; 
done
