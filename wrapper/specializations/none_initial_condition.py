from stages.initial_condition import InitialCondition
from utils.db import insert_initial_condition

class NoneInitialCondition(InitialCondition):
    def validate(self):
        print("No initial condition specified. Skipping initial condition stage.")

    def run(self, event_id):
        print(f"Skipping initial condition stage for event {event_id}.")

        # Insert placeholder initial condition into the database
        insert_initial_condition(
            self.db_connection,
            event_id=event_id,
            seed=None,
            ic_type='none'
        )
