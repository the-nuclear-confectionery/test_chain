class Preequilibrium:
    centrality_estimator = 0.0

    def __init__(self, config, db_connection):
        self.config = config
        self.db_connection = db_connection

    def validate(self, event_id):
        if 'preequilibrium' not in self.config:
            raise ValueError("Preequilibrium configuration is missing.")
        preequilibrium_type = self.config['preequilibrium'].get('type', None)
        if not preequilibrium_type:
            raise ValueError("Preequilibrium type must be specified.")

        print(f"Preequilibrium validation passed for type: {preequilibrium_type}")

    def run(self, event_id):
        raise NotImplementedError("Subclasses should implement this method.")
    
    def get_centrality_estimator(self):
        return self.centrality_estimator
        
    def write_preeq(self, output_file, x, y, eta, e, rhoB, rhoS, rhoQ, ux, uy, ueta, bulk, pixx, pixy, pixeta, piyy, piyeta, pietaeta, stepx, stepy, stepEta, xmin, ymin, etamin):
        """
        Writes the pre-equilibrium data in the CCAKE-compatible format.
        """
        with open(output_file, 'w') as f:
            f.write(f"#0 {stepx} {stepy} {stepEta} 0 {xmin} {ymin} {etamin}\n")
            for i in range(len(x)):
                f.write(f"{x[i]} {y[i]} {eta[i]} {e[i]} {rhoB[i]} {rhoS[i]} {rhoQ[i]} "
                        f"{ux[i]} {uy[i]} {ueta[i]} {bulk[i]} {pixx[i]} {pixy[i]} {pixeta[i]} "
                        f"{piyy[i]} {piyeta[i]} {pietaeta[i]}\n")