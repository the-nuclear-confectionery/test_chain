class Analysis:

    def __init__(self, config, db_connection):
        self.config = config
        self.db_connection = db_connection

    def validate(self, event_id):
        if 'analysis' not in self.config:
            raise ValueError("Analysis configuration is missing.")
        hydro_type = self.config['analysis'].get('type', None)
        if not hydro_type:
            raise ValueError("Analysis type must be specified.")

        print(f"Analysis validation passed for type: {hydro_type}")

    def run(self, event_id):
        raise NotImplementedError("Subclasses should implement this method.")
