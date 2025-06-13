import os
import yaml
import h5py
import numpy as np
import external_packages.freestream as freestream
from stages.preequilibrium import Preequilibrium
from utils.db import insert_preequilibrium

class Freestreaming(Preequilibrium):
    def __init__(self, config, db_connection):
        super().__init__(config, db_connection)

    def validate(self, event_id):
        if self.config['input']['preequilibrium']['type'] != 'freestreaming':
            raise ValueError("Invalid preequilibrium type specified.")
        print("Freestreaming validation passed.")

    def create_temp_config(self, event_id):
        """Create a temporary configuration file for freestreaming."""
        output_dir = os.path.join(self.config['global']['output'], f"event_{event_id}", 'freestream')
        config_dir = os.path.join(self.config['global']['output'], f"event_{event_id}", 'configs')
        os.makedirs(output_dir, exist_ok=True)

        config_data = {
            'grid-max': self.config['global']['grid']['x_max'],
            'tau_fs': self.config['global']['tau_hydro'],
            'label': f"event_{event_id}",
            'output': output_dir
        }

        config_file_path = os.path.join(config_dir, 'freestream_config.yaml')
        with open(config_file_path, 'w') as f:
            yaml.dump(config_data, f)

        print(f"Freestreaming config file created at {config_file_path}")
        return config_file_path

    def run(self, event_id):
        """Run the freestreaming module."""
        config_file_path = self.create_temp_config(event_id)
        with open(config_file_path, 'r') as f:
            parsed_file = yaml.safe_load(f)

        grid_max = parsed_file['grid-max']
        tau_fs = parsed_file['tau_fs']
        label = parsed_file['label']
        results_parent_path = parsed_file['output']

        initial_file = os.path.join(self.config['global']['output'], f"event_{event_id}", 'trento', 'ic.hdf')
        output_path = results_parent_path
        os.makedirs(output_path, exist_ok=True)

        with h5py.File(initial_file, 'r') as event_database:
            for event_name, event_dataset in event_database.items():
                event_path = os.path.join(output_path)
                event_file = os.path.join(event_path, f"fs.dat")

                ic = np.array(event_dataset)
                fs = freestream.FreeStreamer(ic, grid_max, tau_fs)

                x, y, eta = np.meshgrid(np.linspace(-grid_max, grid_max, ic.shape[0]),
                                        np.linspace(-grid_max, grid_max, ic.shape[1]),
                                        [0], indexing='ij')

                e = fs.energy_density()
                rhoB = np.zeros_like(e)
                rhoS = np.zeros_like(e)
                rhoQ = np.zeros_like(e)
                ux = fs.flow_velocity(1)
                uy = fs.flow_velocity(2)
                ueta = np.zeros_like(e)
                bulk = np.zeros_like(e)
                pixx = fs.shear_tensor(0,0)
                pixy = fs.shear_tensor(0,1)
                pixeta = np.zeros_like(e)
                piyy = fs.shear_tensor(1,1)
                piyeta = np.zeros_like(e)
                pietaeta = np.zeros_like(e)

                stepx = 2 * grid_max / e.shape[0]
                stepy = 2 * grid_max / e.shape[1]
                stepEta = 1
                xmin, ymin, etamin = -grid_max, -grid_max, 0

                self.write_preeq(event_file, x.flatten(), y.flatten(), eta.flatten(), e.flatten(), rhoB.flatten(), rhoS.flatten(), rhoQ.flatten(),
                                 ux.flatten(), uy.flatten(), ueta.flatten(), bulk.flatten(), pixx.flatten(), pixy.flatten(), pixeta.flatten(),
                                 piyy.flatten(), piyeta.flatten(), pietaeta.flatten(), stepx, stepy, stepEta, xmin, ymin, etamin)
        
        print(f"Freestreaming completed for event {event_id}.")
        insert_preequilibrium(self.db_connection, event_id, None, self.config['global']['tau_hydro'], "freestreaming")