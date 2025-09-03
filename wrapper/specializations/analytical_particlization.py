import os
import shutil
import random
import subprocess
from stages.particlization import Particlization
from utils.db import insert_particlization


class AnalyticalParticlization(Particlization):

    def validate(self, event_id: int):
        part = self.config['input']['particlization']
        p = part['parameters']

        # Resolve input_file: if absent, try to use hydro output (CCAKE) or error out.
        if p.get('input_file', 'default') == 'default':
            hydro = self.config['input'].get('hydrodynamics', {})
            if hydro.get('type', None) == 'CCAKE':
                surface_file = os.path.join(
                    self.config['global']['output'],
                    f"event_{event_id}",
                    "ccake",
                    "freeze_out.dat",
                )
                p['input_file'] = surface_file
            else:
                raise ValueError(
                    "Analytical particlization requires 'input_file' unless "
                    "hydrodynamics: CCAKE produced freeze_out.dat"
                )
            
        #check defaults
        if p.get('folder', 'default') == 'default':
            p['folder'] =  f"event_{event_id}"
            
        if p.get('df_correction_file', 'default') == 'default':
            p['df_correction_file'] = "input/MHlistall.dat"
        if p.get('basis_files', 'default') == 'default':
            p['basis_files'] = ["input/gl15.dat", "input/gq20.dat"]
        if p.get('ptphi_grid_file', 'default') == 'default':
            p['ptphi_grid_file'] = "input/ptphipoints.dat"
        if p.get('hadron_list_file', 'default') == 'default':
            p['hadron_list_file'] = "input/resoweakPDG2016Plus.dat"
        if p.get('hadron_check_file', 'default') == 'default':
            p['hadron_check_file'] = "input/numbers16p.dat"

        print("Analytical particlization validation passed.")



    def run(self, event_id: int):
        seed = -1
        self.create_temp_config(seed, event_id)

        # Workspace
        tmpdir = os.path.join(self.config['global']['tmp'], f"event_{event_id}", "analytical")
        os.makedirs(tmpdir, exist_ok=True)
        original_input = os.path.join(self.config['global']['basedir'], "models", "analytical_cooperfrye", "input")
        output_dir = os.path.join(tmpdir, "results")
        #os.makedirs(input_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        # copy input directory to executable
        shutil.copytree(original_input, os.path.join(tmpdir, "input"))
        input_dir = os.path.join(tmpdir, "input")
        #mkdir folder inside input
        os.makedirs(os.path.join(input_dir, self.config['input']['particlization']['parameters']['folder']), exist_ok=True)

        # Copy main input file
        p = self.config['input']['particlization']['parameters']
        input_file = p['input_file']
        shutil.copy(input_file, os.path.join(input_dir, os.path.basename(input_file)))

        #copy input file to folder, and rename_it to freeze_out0.dat
        shutil.copy(p['input_file'], os.path.join(input_dir, self.config['input']['particlization']['parameters']['folder'], "freeze_out0.dat"))

        # Copy executable
        exe_src = os.path.join(self.config['global']['basedir'], "models", "analytical_cooperfrye", "fo")
        print(f"Analytical executable: {exe_src}")
        exe_tmp = os.path.join(tmpdir, os.path.basename(exe_src))
        shutil.copy(exe_src, exe_tmp)
        os.chmod(exe_tmp, 0o755)


        # Copy config into tmpdir
        cfg_src = os.path.join(self.config['global']['output'], f"event_{event_id}", "configs", "analytical_parameters.dat")
        shutil.copy(cfg_src, tmpdir)
        #copy to input and rename to dfinput.dat
        shutil.copy(cfg_src, os.path.join(input_dir, "dfinput.dat"))
        #create tmp output tmp/out/event_id
        os.makedirs(os.path.join(tmpdir, "out", f"event_{event_id}"), exist_ok=True)


        # Run
        command = ["./fo"]
        print(f"Running command: {' '.join(command)} in {tmpdir}")
        subprocess.run(command, cwd=tmpdir)

        # Insert DB record (no explicit sample count here → use -1)
        insert_particlization(
            self.db_connection,
            event_id=event_id,
            seed=seed,
            particlization_type=self.config['input']['particlization']['type'],
            nsamples=-1,
        )

        # Copy from tmp output to final output
        shutil.copytree(
            os.path.join(tmpdir, "out", f"event_{event_id}"),
            os.path.join(self.config['global']['output'], f"event_{event_id}", "analytical"),
            dirs_exist_ok=True
        )
        #rename the corrected and uncorrected
        os.rename(
            os.path.join(self.config['global']['output'], f"event_{event_id}", "analytical", "evv2c_dNdphidpp.dat"),
            os.path.join(self.config['global']['output'], f"event_{event_id}", "analytical", "corrected_spectra.dat")
        )
        os.rename(
            os.path.join(self.config['global']['output'], f"event_{event_id}", "analytical", "evv2_dNdphidpp.dat"),
            os.path.join(self.config['global']['output'], f"event_{event_id}", "analytical", "uncorrected_spectra.dat")
        )

    def create_temp_config(self, seed: int, event_id: int):

        out_cfg_dir = os.path.join(self.config['global']['output'], f"event_{event_id}", "configs")
        os.makedirs(out_cfg_dir, exist_ok=True)
        out_cfg_file = os.path.join(out_cfg_dir, "analytical_parameters.dat")

        p = self.config['input']['particlization']['parameters']
        input_file = p.get('input_file')
        folder_path  = p.get('folder')
        typeofequations = p.get('equations_version', 'v2')
        basis = p.get('basis_files', [])
        basis_pt  = basis[0] if len(basis) > 0 else ''
        basis_phi = basis[1] if len(basis) > 1 else ''
        decays = p.get('decays', 3)
        er = p.get('event_range', [0, 0])
        ev0, ev1 = (er + [0, 0])[:2] if isinstance(er, list) else (0, 0)
        freezeouttemp = p.get('temperature', 0.150)
        df_cor_file = p.get('df_correction_file', '')
        ptphi_grid = p.get('ptphi_grid_file', '')
        hadronl_list = p.get('hadron_list_file', '')      
        had_check = p.get('hadron_check_file', '')

        with open(out_cfg_file, 'w') as f:
            f.write(f"typeofequations: {typeofequations}\n")
            f.write(f"folder: {folder_path}\n")
            f.write(f"range(pt,phi): {basis_pt} {basis_phi}\n")
            f.write(f"decays: {decays}\n")
            f.write(f"rangeofevents: {ev0} {ev1}\n")
            f.write(f"freezeouttemp: {freezeouttemp}\n")
            f.write(f"df_cor_file: {df_cor_file}\n")
            f.write(f"range(ptmax,ptstepsize,phisteps): {ptphi_grid}\n")
            f.write(f"hadronl_list: {hadronl_list}\n")
            f.write(f"had_check: {had_check}\n")