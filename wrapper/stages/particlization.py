class Particlization:

    def __init__(self, config, db_connection):
        self.config = config
        self.db_connection = db_connection

    def validate(self, event_id):
        if 'Particlization' not in self.config:
            raise ValueError("Particlization configuration is missing.")
        hydro_type = self.config['particlization'].get('type', None)
        if not hydro_type:
            raise ValueError("Particlization type must be specified.")

        print(f"Particlization validation passed for type: {hydro_type}")

    def run(self, event_id):
        raise NotImplementedError("Subclasses should implement this method.")
