from stages.analysis import Analysis
from utils.db import insert_analysis

class NoneAnalysis(Analysis):
    def validate(self, event_id):
        print("No analysis specified. Skipping analysis stage.")
        # No validation needed for "none" type
    
    def run(self, event_id):
        print(f"Skipping analysis stage for event {event_id}.")
    
        # Insert placeholder analysis into the database
        insert_analysis(
            self.db_connection,
            event_id=event_id,
            analysis_type=None,
            rapidity_cut=None,
        )