from stages.afterburner import Afterburner
from utils.db import insert_afterburner

class NoneAfterburner(Afterburner):
    def validate(self, event_id):
        print("No afterburner specified. Skipping afterburner stage.")
        # No validation needed for "none" type
    
    def run(self, event_id):
        print(f"Skipping afterburner stage for event {event_id}.")
    
        # Insert placeholder afterburner into the database
        insert_afterburner(
            self.db_connection,
            event_id=event_id,
            seed=None,  # No seed for "none"
            afterburner_type=None,
        )