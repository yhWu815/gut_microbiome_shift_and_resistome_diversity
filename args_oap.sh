#!/bin/bash
#!/bin/bash
echo "---Activate args_oap conda environment---"
echo "${USER}"

source /home/${USER}/miniconda3/etc/profile.d/conda.sh
conda activate args_oap

echo '---ARGS OAP---'
declare -A example=(["input dir"]="meta_clean" ["output dir"]="args_output" ["file format"]="fastq" ["number of threads"]="24")

args_oap stage_one -i ${example["input dir"]} -o ${example["output dir"]} -f ${example["file format"]} -t ${example["number of threads"]}
args_oap stage_two -i ${example["output dir"]} -t ${example["number of threads"]}