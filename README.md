# Installation and Usage Guide

## 1. Install Anaconda or Miniconda

If you haven't installed Anaconda or Miniconda, do so first. On a cluster, you can check for an available module with:

module spider anaconda

Then, load the specific version with:

module load <module_name>

If this is your first time using Anaconda, you will need to initialize it by running:

conda init

Then, restart your terminal for the changes to take effect.

---

## 2. Install the Conda Environment

From the main project directory, install the Conda environment using:

conda env create -f environment.yaml

Activate the environment with:

conda activate ctop

---

## 3. Install the Project

Run the installation script:

sh install

Note: This process may take a while. If you are on a cluster where the main node has limited resources, consider using a SLURM job.

After installation, deactivate and reactivate the environment to apply environment variable settings:

conda deactivate
conda activate ctop

The installation sets up environment variables, for example, $WORKDIR will always point to the project root directory.

---

## 4. Run a Single Event with the Wrapper

Ensure the Conda environment is active:

conda activate ctop

Then run:

python $WORKDIR/wrapper/main.py <event_number> <path/to/database> <path/to/yml/config>

- The SQL database will store event metadata. If the specified database does not exist, it will be created automatically.
- Example YAML configuration files can be found in:

$WORKDIR/examples

---

## 5. Run Multiple Events with SLURM

To launch a batch job for multiple events, use the run_wrapper.sh script located in:

$WORKDIR/scripts

Modify the script as needed and submit it using SLURM.

