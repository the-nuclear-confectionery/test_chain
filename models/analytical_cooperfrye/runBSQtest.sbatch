#!/bin/bash                                                                                                                                                          
#SBATCH --partition=qgp             # Partition (job queue)                                                                                                          
#SBATCH --job-name=IDNUM_BINN            # Assign an short name to your job                                                                                                                                                                                                                                         
#SBATCH --requeue                    # Return job to the queue if preempted                                                                                          
#SBATCH --nodes=1                    # Number of nodes you require                                                                                                   
#SBATCH --ntasks=1                   # Total # of tasks across all nodes                                                                                             
#SBATCH --cpus-per-task=1            # Cores per task (>1 if multithread tasks)                                                                                      
#SBATCH --mem=600                   # Real memory (RAM) required (MB)                                                                                                
#SBATCH --time=30:00:00              # Total run time limit (HH:MM:SS)                                                                                               
#SBATCH --output=BSQtest_decays/TYP_BINN.out                                                                                                                                   
#SBATCH --error=BSQtest_decays/TYP_BINN.err                                                                                                                                    
#SBATCH --export=ALL                 # Export you current env to the job env                                                                                         


# set data and working directories                                                                                                                                   
cd /projects/jnorhos/v-USPhydro2/df



srun decays/reso BSQtest_decays/ev0bsqsbvc_dNdphidpp.dat BSQtest_decays/ev0dbsqsbvc_dNdphidpp.dat reso16p.inp > logfile.txt

sleep 3
sacct --format NTasks,MaxRSS,Elapsed,AveRSS,AveCPU 
sleep 2
                                                                                                                                                       
                                                                                                                                         


