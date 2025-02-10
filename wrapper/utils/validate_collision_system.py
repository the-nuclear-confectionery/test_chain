import sys
from utils.db import insert_collision_system

class CollisionSystem:
    def __init__(self, config, db_connection):
        self.config = config
        self.db_connection = db_connection
        self.ic_type = config['input'].get('initial_conditions', {}).get('type', 'none').lower()

    def get_species_from_globals(self):
        """Retrieve species information from the globals section, if available."""
        return (
            self.config['global'].get('ion_A'),
            self.config['global'].get('ion_B')
        )

    def get_species_from_ic(self):
        """Retrieve species information from the initial conditions section, if available."""
        if self.ic_type != 'none':
            parameters = self.config['input'].get('initial_conditions', {}).get('parameters', {})
            return (
                parameters.get('ion_A', {}).get('species'),
                parameters.get('ion_B', {}).get('species')
            )
        return (None, None)

    def resolve_species(self):
        """Resolve collision system species by comparing globals and initial conditions."""
        globals_ion_A, globals_ion_B = self.get_species_from_globals()
        ic_ion_A, ic_ion_B = self.get_species_from_ic()

        # Case 1: If both globals and iC provide species info
        if ic_ion_A and ic_ion_B:
            if globals_ion_A and globals_ion_B:
                if (globals_ion_A == ic_ion_A and globals_ion_B == ic_ion_B):
                    # iC and globals agree, continue
                    return ic_ion_A, ic_ion_B
                else:
                    # iC and globals disagree, warn and overwrite globals
                    print(f"Warning: Initial conditions species ({ic_ion_A}, {ic_ion_B}) differ from globals "
                          f"({globals_ion_A}, {globals_ion_B}). Overwriting globals with initial conditions.")
                    return ic_ion_A, ic_ion_B
            else:
                # No globals, use iC directly
                return ic_ion_A, ic_ion_B

        # Case 2: If no species info in iC but globals has them
        if globals_ion_A and globals_ion_B:
            return globals_ion_A, globals_ion_B

        # Case 3: Neither iC nor globals have species info, error
        print("Error: No collision system species specified in either initial conditions or globals.")
        sys.exit(1)

    def insert_collision_system(self):
        """Insert collision system into the database after resolving species."""
        ion_A, ion_B = self.resolve_species()
        return insert_collision_system(
            self.db_connection,
            ion_A,
            ion_B,
            self.config['global']['energy']['value'],
            self.config['global']['energy']['unit']
        )
