class Afterburner:

    def __init__(self, config, db_connection):
        self.config = config
        self.db_connection = db_connection

    def validate(self, event_id):
        if 'Afterburner' not in self.config:
            raise ValueError("Afterburner configuration is missing.")
        afterbuner_type = self.config['afterbuner'].get('type', None)
        if not afterbuner_type:
            raise ValueError("Afterburner type must be specified.")

        print(f"Afterburner validation passed for type: {afterbuner_type}")

    def run(self, event_id):
        raise NotImplementedError("Subclasses should implement this method.")
