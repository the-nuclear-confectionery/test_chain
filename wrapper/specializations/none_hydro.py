from stages.hydrodynamics import Hydrodynamics
from utils.db import insert_hydro

class NoneHydro(Hydrodynamics):
    def validate(self, event_id):
        print("No hydrodynamics specified. Skipping hydrodynamics stage.")

    def run(self, event_id):
        print(f"Skipping hydrodynamics stage for event {event_id}.")

        # Insert placeholder overlay into the database
        insert_hydro(
            self.db_connection,
            event_id=event_id,
            initial_time=None,
            freeze_out_temperature=None,
            dimensions='none',
            hydro_type='none',

        )