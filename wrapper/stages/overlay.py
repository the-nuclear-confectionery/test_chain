class Overlay:
    centrality_estimator = 0.0

    def __init__(self, config, db_connection):
        self.config = config
        self.db_connection = db_connection

    def validate(self, event_id):
        if 'overlay' not in self.config:
            raise ValueError("Overlay configuration is missing.")
        overlay_type = self.config['overlay'].get('type', None)
        if not overlay_type:
            raise ValueError("Overlay type must be specified.")

        print(f"Overlay validation passed for type: {overlay_type}")

    def run(self, event_id):
        raise NotImplementedError("Subclasses should implement this method.")
    
    def get_centrality_estimator(self):
        return self.centrality_estimator
        
