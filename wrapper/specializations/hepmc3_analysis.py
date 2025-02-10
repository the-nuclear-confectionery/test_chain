from stages.analysis import Analysis
import os
from utils.db import insert_analysis


class HepMC3Analysis(Analysis):

    def __init__(self, config, db_connection):
        self.config = config
        self.db_connection = db_connection

    def validate(self, event_id):
        # Check for afterburner module
        if self.config['input']['afterburner']['type'] == None:
            if self.config['input']['analysis']['parameters']['input_file'] == "default":
                raise ValueError("Afterburner is disabled, input file must be specified.")

        elif self.config['input']['afterburner']['type'] == 'smash':
            # Check if HepMC_treeroot output is enabled
            formats = self.config['input']['afterburner']['parameters']['Output']['Particles']['Format']

            if 'HepMC_treeroot' not in formats:
                raise ValueError("HepMC_treeroot output is not enabled in SMASH afterburner configuration.")
            else:
                self.config['input']['analysis']['parameters']['input_file'] =  self.config['global']['output'] + '/event_' + str(event_id) + '/smash/SMASH_HepMC_particles.root'
        
        if self.config['input']['analysis']['parameters']['output_file'] == 'default':
            self.config['input']['analysis']['parameters']['output_file'] = self.config['global']['output'] + '/event_' + str(event_id) + '/histo_event_' + str(event_id) + '.root'

        print(f"Validation of HepMC3 configuration completed successfully.")


    

    def run(self, event_id):
        # Simple placeholder for the actual implementation
        print(f"Running HepMC3 analysis for event ID: {event_id}")

        output_file = self.config['input']['analysis']['parameters']['output_file']
        input_file = self.config['input']['analysis']['parameters']['input_file']
        rapidity_cut = self.config['input']['analysis']['parameters']['rapidity_cut']

        analysis_path= self.config['global']['basedir'] + '/analysis/build/spectra_hepmc3.exe'
        # Run the analysis executable with the provided parameters
        command = f"{analysis_path} {input_file} {output_file} {rapidity_cut}"
        print(f"Executing command: {command}")
        os.system(command)
        # Insert analysis results into the database
        insert_analysis(
            self.db_connection,
            event_id=event_id,
            analysis_type='HepMC3',
            rapidity_cut=rapidity_cut,
        )
