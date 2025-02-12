import warnings
import os
import subprocess
import re
import random
import shutil
from utils.db import insert_initial_condition
from stages.initial_condition import InitialCondition

class AmptInitialCondition(InitialCondition):
    def validate(self):
        """Perform validation specific to AMPT if it's enabled."""
        if self.config['input']['initial_conditions']['type'] != 'ampt':
            return  # Skip validation if AMPT is not the active IC
        # read energy EFRM from global config
        EFRM = self.config['global']['energy']['value']
        #check unit 
        if self.config['global']['energy']['unit'] == 'GeV':
            EFRM = EFRM
        elif self.config['global']['energy']['unit'] == 'TeV':
            EFRM = EFRM * 1000
        elif self.config['global']['energy']['unit'] == 'MeV':
            EFRM = EFRM / 1000
        else:
            raise ValueError(f"Unknown energy unit: {self.config['global']['energy']['units']}")
        #convert to int and save
        self.config['input']['initial_conditions']['parameters']['EFRM'] = int(EFRM)
        
        # read ion A and ion B
        ion_mapping = {
            "p":  (1, 1),
            "d":  (2, 1),
            "Cu": (64, 29),
            "Au": (197, 79),
            "Pb": (208, 82),
            "Xe": (131, 54),
            "Ru": (96, 40),
            "U":  (238, 92),
        }
        
        A = self.config['input']['initial_conditions']['parameters']['ion_A']['species']
        B = self.config['input']['initial_conditions']['parameters']['ion_B']['species']
        try:
            IAP, IZP = ion_mapping[A]
            IAT, IZT = ion_mapping[B]
        except KeyError:
            raise ValueError(f"Unknown ion type: {A} or {B}")
        
        params = self.config['input']['initial_conditions']['parameters']
        params['IAP'] = IAP
        params['IZP'] = IZP
        params['IAT'] = IAT
        params['IZT'] = IZT

    
        
        print("Base validation for AMPT completed.")

    def generate_seeds(self):
        """
        Check if seeds are provided; if set to 'random' generate new random seeds.
        For AMPT, we treat the HIJING and parton cascade seeds separately.
        """
        params = self.config['input']['initial_conditions']['parameters']
        
        # HIJING seed
        if params.get('hijing_seed') == 'random':
            hijing_seed = random.randint(0, 2**31 - 1)
            print(f"Generated random HIJING seed for AMPT: {hijing_seed}")
        else:
            hijing_seed = params.get('hijing_seed')

        # Parton cascade seed
        if params.get('parton_cascade_seed') == 'random':
            parton_cascade_seed = random.randint(0, 2**31 - 1)
            print(f"Generated random parton cascade seed for AMPT: {parton_cascade_seed}")
        else:
            parton_cascade_seed = params.get('parton_cascade_seed')

        return hijing_seed, parton_cascade_seed

    def create_temp_config(self, event_id):
        """
        Create a temporary AMPT configuration file using the YAML-provided parameters.
        The file is written with one parameter per line in the order expected by AMPT.
        """
        if self.config['input']['initial_conditions']['type'] != 'ampt':
            raise ValueError("create_temp_config should only be called when AMPT is the active IC.")

        params = self.config['input']['initial_conditions']['parameters']
        # Define output directories for AMPT
        output_dir = os.path.join(self.config['global']['output'], f"event_{event_id}", 'ampt')
        os.makedirs(output_dir, exist_ok=True)
        config_dir = os.path.join(self.config['global']['output'], f"event_{event_id}", 'configs')
        os.makedirs(config_dir, exist_ok=True)
        config_file_path = os.path.join(config_dir, 'input.ampt')

        # Update seeds (if needed)
        hijing_seed, parton_cascade_seed = self.generate_seeds()
        params['hijing_seed'] = hijing_seed
        params['parton_cascade_seed'] = parton_cascade_seed

        with open(config_file_path, 'w') as f:
            # Write parameters in the required order:
            f.write(f"{params.get('EFRM')}\n")         # EFRM: sqrt(S_NN) in GeV if FRAME is CMS
            f.write(f"{params.get('FRAME')}\n")        # FRAME
            f.write(f"{params.get('PROJ')}\n")           # PROJ
            f.write(f"{params.get('TARG')}\n")           # TARG
            f.write(f"{params.get('IAP')}\n")            # IAP (projectile A number)
            f.write(f"{params.get('IZP')}\n")             # IZP (projectile Z number)
            f.write(f"{params.get('IAT')}\n")            # IAT (target A number)
            f.write(f"{params.get('IZT')}\n")             # IZT (target Z number)
            f.write(f"{params.get('NEVNT')}\n")            # NEVNT (total number of events)
            f.write(f"{params.get('BMIN')}\n")           # BMIN (minimum impact parameter in fm)
            f.write(f"{params.get('BMAX')}\n")           # BMAX (maximum impact parameter in fm)
            f.write(f"{params.get('ISOFT')}\n")            # ISOFT (select Default AMPT or String Melting)
            f.write(f"{params.get('NTMAX')}\n")          # NTMAX (number of timesteps)
            f.write(f"{params.get('DT')}\n")             # DT (timestep in fm)
            f.write(f"{params.get('PARJ41')}\n")        # PARJ(41): parameter a in Lund symmetric splitting function
            f.write(f"{params.get('PARJ42')}\n")        # PARJ(42): parameter b in Lund symmetric splitting function
            f.write(f"{params.get('popcorn')}\n")          # popcorn flag (1: yes; 0: no) for netbaryon stopping
            f.write(f"{params.get('PARJ5')}\n")          # PARJ(5) to control BMBbar vs BBbar in popcorn
            f.write(f"{params.get('shadowing_flag')}\n")   # shadowing flag (1 for yes, 0 for no)
            f.write(f"{params.get('quenching_flag')}\n")   # quenching flag (0: no; 1: yes)
            f.write(f"{params.get('quenching_parameter')}\n")  # quenching parameter -dE/dx (GeV/fm)
            f.write(f"{params.get('p0_cutoff')}\n")      # p0 cutoff in HIJING for minijet productions
            f.write(f"{params.get('parton_screening_mass')}\n")  # parton screening mass in fm^(-1)
            f.write(f"{params.get('IZPC')}\n")             # IZPC: flag for forward-angle parton scatterings
            f.write(f"{params.get('alpha')}\n")         # alpha in parton cascade
            f.write(f"{params.get('dpcoal')}\n")         # dpcoal in GeV
            f.write(f"{params.get('drcoal')}\n")         # drcoal in fm
            f.write(f"{params.get('ihjsed')}\n")          # ihjsed: HIJING seed flag
            f.write(f"{params.get('hijing_seed')}\n")  # hijing_seed
            f.write(f"{params.get('parton_cascade_seed')}\n")  # parton_cascade_seed
            f.write(f"{params.get('k0s_decay_flag')}\n")    # flag for K0s weak decays
            f.write(f"{params.get('phi_decay_flag')}\n")    # flag for phi decays at end of hadron cascade
            f.write(f"{params.get('pi0_decay_flag')}\n")    # flag for pi0 decays at end of hadron cascade
            f.write(f"{params.get('oscar_output')}\n")      # optional OSCAR output flag
            f.write(f"{params.get('perturbative_deuteron_flag')}\n")  # perturbative deuteron flag
            f.write(f"{params.get('perturbative_deuteron_factor')}\n")  # integer factor for perturbative deuterons
            f.write(f"{params.get('deuteron_cross_section_choice')}\n") # cross section choice for deuteron reactions
            f.write(f"{params.get('pttrig')}\n")        # Pt threshold in GeV for minijet requirement
            f.write(f"{params.get('maxmiss')}\n")       # maxmiss: maximum # of tries to repeat HIJING event
            f.write(f"{params.get('initial_final_radiation_flag')}\n")  # flag on initial and final state radiation
            f.write(f"{params.get('kt_kick_flag')}\n")     # flag on Kt kick
            f.write(f"{params.get('quark_pair_embedding_flag')}\n")  # flag to turn on quark pair embedding
            # Embedded quark initial momentum (two values, separated by a comma)
            px = params.get('embedded_quark_initial_px')
            py = params.get('embedded_quark_initial_py')
            f.write(f"{px}, {py}\n")
            # Embedded quark initial position (two values, separated by a comma)
            x = params.get('embedded_quark_initial_x')
            y = params.get('embedded_quark_initial_y')
            f.write(f"{x}, {y}\n")
            # nsembd, psembd, tmaxembd on one line (comma separated)
            nsembd = params.get('nsembd')
            psembd = params.get('psembd')
            tmaxembd = params.get('tmaxembd')
            f.write(f"{nsembd}, {psembd}, {tmaxembd}\n")
            f.write(f"{params.get('ishadow')}\n")        # ishadow: flag for modifying shadowing
            f.write(f"{params.get('dshadow')}\n")       # dshadow: factor to modify nuclear shadowing
            f.write(f"{params.get('iphirp')}\n")          # iphirp: flag for random orientation of reaction plane


        print(f"AMPT config file created at {config_file_path}")
        return config_file_path
    

    def run(self, event_id):
        """Run the AMPT IC generator for a specific event in an isolated temporary environment."""
        if self.config['input']['initial_conditions']['type'] != 'ampt':
            print("Skipping AMPT execution because it is not the active initial condition.")
            return

        # Step 1: Create a temporary working directory for this event
        tmp_ampt_dir = os.path.join(self.config['global']['tmp'], f"event_{event_id}", 'ampt')
        os.makedirs(tmp_ampt_dir, exist_ok=True)
        print(f"Created temporary working directory: {tmp_ampt_dir}")

        # Step 2: Copy the entire models/AMPT directory to the temporary working directory
        ampt_base_dir = self.config['global']['basedir'] + '/models/AMPT'
        print(f"Copying AMPT base directory from {ampt_base_dir} to {tmp_ampt_dir}")
        shutil.copytree(ampt_base_dir, tmp_ampt_dir, dirs_exist_ok=True)
        print(f"Copied the entire AMPT directory from {ampt_base_dir} to {tmp_ampt_dir}")
        # remove the ana directory
        #shutil.rmtree(os.path.join(tmp_ampt_dir, 'ana'))

        # create empty ana directory
        os.makedirs(os.path.join(tmp_ampt_dir, 'ana'), exist_ok=True)

        # Step 3: Copy the generated config file to the temporary directory
        config_file_path = self.create_temp_config(event_id)
        shutil.copy(config_file_path, os.path.join(tmp_ampt_dir, 'input.ampt'))

        # Step 4: Run AMPT using ./exec in the temporary directory
        os.chdir(tmp_ampt_dir)
        print(f"Running AMPT for event {event_id} in {tmp_ampt_dir} ...")
        try:
            subprocess.run(['./exec'], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error during AMPT execution for event {event_id}: {e}")
            return

        # Step 5: Copy results to the final output directory
        output_dir = os.path.join(self.config['global']['output'], f"event_{event_id}", 'ampt')
        os.makedirs(output_dir, exist_ok=True)
        #copy parton-collisionsHistory.dat,parton-initial-afterPropagation.dat
        # ampt.dat dnde_ch.dat

        shutil.copy(os.path.join(tmp_ampt_dir, 'ana/parton-collisionsHistory.dat'), output_dir)
        shutil.copy(os.path.join(tmp_ampt_dir, 'ana/parton-initial-afterPropagation.dat'), output_dir)
        shutil.copy(os.path.join(tmp_ampt_dir, 'ana/ampt.dat'), output_dir)
        shutil.copy(os.path.join(tmp_ampt_dir, 'ana/dnde_ch.dat'), output_dir)


        print(f"AMPT output copied to {output_dir}")

        # Step 6: Clean up the temporary working directory
        #shutil.rmtree(tmp_ampt_dir)

        #insert initial condition into db
        hijing_seed = self.config['input']['initial_conditions']['parameters']['hijing_seed']
        parton_cascade_seed = self.config['input']['initial_conditions']['parameters']['parton_cascade_seed']
        print(f"Writing initial condition to db with hijing_seed={hijing_seed}, parton_cascade_seed={parton_cascade_seed}")
        insert_initial_condition(self.db_connection, event_id, hijing_seed, parton_cascade_seed, 0, 0, 0, "AMPT - e2 = cascade_seed")

        print(f"AMPT execution completed for event {event_id}.")


