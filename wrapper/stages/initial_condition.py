import os
from utils.db import insert_initial_condition, insert_event, insert_collision_system

class InitialCondition:
    centrality_estimator = 0.0

    def __init__(self, config, db_connection):
        self.config = config
        self.db_connection = db_connection

    
    def get_centrality_estimator(self):
        return self.centrality_estimator
    
    def write_ic(self, output_file, x, y, eta, e, rhoB, rhoS, rhoQ, ux, uy, ueta, 
                  bulk, pixx, pixy, pixeta, piyy, piyeta, pietaeta, 
                  stepx, stepy, stepEta, xmin, ymin, etamin):
        """
        Writes the ic in the CCAKE-compatible format.

        Parameters:
            output_file (str): Path to the output file.
            x, y, eta (list of floats): Coordinates of the grid points.
            e, rhoB, rhoS, rhoQ (list of floats): Thermodynamic quantities.
            ux, uy, ueta (list of floats): Contravariant flow velocity components.
            bulk, pixx, pixy, pixeta, piyy, piyeta, pietaeta (list of floats): 
                Viscous and hydrodynamic properties.
            stepx, stepy, stepEta (float): Step sizes in x, y, and eta directions.
            xmin, ymin, etamin (float): Minimum values for x, y, and eta coordinates.
        """
        # Open the output file for writing
        with open(output_file, 'w') as f:
            # Write the header line
            f.write(f"#0 {stepx} {stepy} {stepEta} 0 {xmin} {ymin} {etamin}\n")

            # Iterate over the lists of quantities and write each row
            for i in range(len(x)):
                f.write(f"{x[i]} {y[i]} {eta[i]} {e[i]} {rhoB[i]} {rhoS[i]} {rhoQ[i]} "
                        f"{ux[i]} {uy[i]} {ueta[i]} {bulk[i]} {pixx[i]} {pixy[i]} {pixeta[i]} "
                        f"{piyy[i]} {piyeta[i]} {pietaeta[i]}\n")
        
