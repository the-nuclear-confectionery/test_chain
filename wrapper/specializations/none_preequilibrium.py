from stages.preequilibrium import Preequilibrium
from utils.db import insert_preequilibrium

class NonePreequilibrium(Preequilibrium):
    def validate(self, event_id):
        print("No Preequilibrium specified. Skipping preequilibrium stage.")

    def run(self, event_id):
        print(f"Skipping Preequilibrium stage for event {event_id}.")
        # Insert placeholder overlay into the database
        insert_preequilibrium(
            self.db_connection,
            event_id=event_id,
            seed=None,  # No seed for "none"
            time=None,
            preequilibrium_type='none'
        )

