from stages.overlay import Overlay
from utils.db import insert_overlay

class NoneOverlay(Overlay):
    def validate(self, event_id):
        print("No overlay specified. Skipping overlay stage.")

    def run(self, event_id):
        print(f"Skipping overlay stage for event {event_id}.")
        # Insert placeholder overlay into the database
        insert_overlay(
            self.db_connection,
            event_id=event_id,
            seed=None,  # No seed for "none"
            eps2=None,
            eps3=None,
            eps4=None,
            eps5=None,
            overlay_type='none'
        )


