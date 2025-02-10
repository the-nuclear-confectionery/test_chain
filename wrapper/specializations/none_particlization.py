from stages.particlization import Particlization
from utils.db import insert_particlization

class NoneParticlization(Particlization):
    def validate(self, event_id):
        print("No particlization specified. Skipping particlization stage.")
        # No validation needed for "none" type
    
    def run(self, event_id):
        print(f"Skipping particlization stage for event {event_id}.")
    
        # Insert placeholder particlization into the database
        insert_particlization(
            self.db_connection,
            event_id=event_id,
            seed=None,  # No seed for "none"
            particlization_type=None,
            nsamples=None,
        )