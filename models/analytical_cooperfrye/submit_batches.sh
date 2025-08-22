#!/bin/bash

mkdir -p $1

sbatch <<EOT
#!/bin/bash
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 10:00:00
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4550
#SBATCH --array=0-$3
#SBATCH --output="$1/job_%a.out"
#SBATCH --error="$1/job_%a.err"

echo This is thread = "\${SLURM_ARRAY_TASK_ID}"

cen=$1
mkdir -p \${cen}

INPUT_PARAMETERS_PATH="./input/dfinput_\${SLURM_ARRAY_TASK_ID}.dat"
INPUT_PARAMETERS_FILE="\$(basename -- \${INPUT_PARAMETERS_PATH})"
echo "Check: INPUT_PARAMETERS_PATH="\${INPUT_PARAMETERS_PATH}
echo "Check: INPUT_PARAMETERS_FILE="\${INPUT_PARAMETERS_FILE}

SLURM_SCRIPT_HS_PATH="$2/ev\${SLURM_ARRAY_TASK_ID}"
echo "Check: SLURM_SCRIPT_HS_PATH="\${SLURM_SCRIPT_HS_PATH}

cat input/dfinput.dat \
  | sed "s:SLURM_SCRIPT_HS_PATH:\${SLURM_SCRIPT_HS_PATH}:g" \
  > \${INPUT_PARAMETERS_PATH}

RESULTS_DIRECTORY=\${cen}/ev\${SLURM_ARRAY_TASK_ID}

mkdir -p \${RESULTS_DIRECTORY}

# save settings file for reference
cp \${INPUT_PARAMETERS_PATH} \${RESULTS_DIRECTORY}

./fo \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}
echo "Next step: run ./fo" \${INPUT_PARAMETERS_FILE} \${RESULTS_DIRECTORY}

exit 0
EOT
