#!/bin/bash
#-------------------------------------------------------------------------------
# set up results directory
RESULTS_DIRECTORY=out/hyper_0047_0000/Norm_69/0/shearBulk/trentoAS
mkdir -p $RESULTS_DIRECTORY

sbatch <<EOT
#!/bin/bash
#SBATCH --output="$RESULTS_DIRECTORY/job.out"
#SBATCH -A qgp
#SBATCH -p qgp
#SBATCH -t 06:00:00
#-------------------------------------------------------------------------------
# set up input parameters file
INPUT_PARAMETERS_FILE=../input/dfinput.dat
cp \${INPUT_PARAMETERS_FILE} $RESULTS_DIRECTORY

#-------------------------------------------------------------------------------
#OUTPUT_FILE=out/..
#cp \${OUTPUT_FILE} $RESULTS_DIRECTORY
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# run the code
echo "Running command: ./fo" \${INPUT_PARAMETERS_FILE} $RESULTS_DIRECTORY
./fo \${INPUT_PARAMETERS_FILE} $RESULTS_DIRECTORY

exit 0
EOT
