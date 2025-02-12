import sqlite3

def initialize_database(db_path):
    connection = sqlite3.connect(db_path)
    cursor = connection.cursor()

    # Create main events table
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS events (
        event_id INTEGER PRIMARY KEY,
        collision_id INTEGER,
        centrality_estimator REAL,
        output_path TEXT,
        ic_type TEXT,
        overlay_type TEXT,
        hydro_type TEXT,
        particlization_type TEXT,
        afterburner_type TEXT,
        analysis_type TEXT,
        yaml_config_path TEXT,
        FOREIGN KEY(collision_id) REFERENCES collision_system(collision_id)
    )
    ''')

    # Create initial conditions table
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS initial_conditions (
        event_id INTEGER PRIMARY KEY,
        seed INTEGER,
        eps2 REAL,
        eps3 REAL,
        eps4 REAL,
        eps5 REAL,
        ic_type TEXT,
        FOREIGN KEY(event_id) REFERENCES events(event_id)
    )
    ''')

    # Create overlay table with same structure as initial conditions
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS overlays (
        event_id INTEGER PRIMARY KEY,
        seed INTEGER,
        eps2 REAL,
        eps3 REAL,
        eps4 REAL,
        eps5 REAL,
        overlay_type TEXT,
        FOREIGN KEY(event_id) REFERENCES events(event_id)
    )
    ''')

    # Create hydro table
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS hydro (
        event_id INTEGER PRIMARY KEY,
        hydro_type TEXT,
        dimensions INTEGER,
        FOREIGN KEY(event_id) REFERENCES events(event_id)
    )
    ''')

    # Create particlization table
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS particlization (
        event_id INTEGER PRIMARY KEY,
        particlization_type TEXT,
        seed INTEGER,
        nsamples INTEGER,
        FOREIGN KEY(event_id) REFERENCES events(event_id)
    )
    ''')

    #create afterburner table
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS afterburner (
        event_id INTEGER PRIMARY KEY,
        afterburner_type TEXT,
        seed INTEGER,
        FOREIGN KEY(event_id) REFERENCES events(event_id)
    )
    ''')

    #create analysis table
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS analysis (
        event_id INTEGER PRIMARY KEY,
        analysis_type TEXT,
        rapidity_cut REAL,
        FOREIGN KEY(event_id) REFERENCES events(event_id)
    )
    ''')

    # Create collision system table
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS collision_system (
        collision_id INTEGER PRIMARY KEY,
        projectile TEXT,
        target TEXT,
        energy REAL,
        energy_units TEXT CHECK(energy_units IN ('GeV', 'TeV'))
    )
    ''')

    connection.commit()
    return connection

def insert_collision_system(connection, projectile, target, energy, energy_units):
    cursor = connection.cursor()
    cursor.execute('''
    SELECT collision_id FROM collision_system 
    WHERE projectile = ? AND target = ? AND energy = ? AND energy_units = ?
    ''', (projectile, target, energy, energy_units))
    result = cursor.fetchone()

    if result:
        return result[0]  # Return existing collision_id

    cursor.execute('''
    INSERT INTO collision_system (projectile, target, energy, energy_units) 
    VALUES (?, ?, ?, ?)
    ''', (projectile, target, energy, energy_units))
    connection.commit()
    return cursor.lastrowid  # Return new collision_id

def insert_event(connection, event_id, output_path, collision_id, yaml_config_path):
    cursor = connection.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO events (event_id, collision_id,output_path , yaml_config_path)
    VALUES (?, ?, ?, ?)
    ''', (event_id, collision_id, output_path, yaml_config_path))
    connection.commit()

def insert_initial_condition(connection, event_id, seed, eps2, eps3, eps4, eps5, ic_type):
    cursor = connection.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO initial_conditions (event_id, seed, eps2, eps3, eps4, eps5, ic_type)
    VALUES (?, ?, ?, ?, ?, ?, ?)
    ''', (event_id, seed, eps2, eps3, eps4, eps5, ic_type))
    connection.commit()

def insert_overlay(connection, event_id, seed, eps2, eps3, eps4, eps5, overlay_type):
    cursor = connection.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO overlays (event_id, seed, eps2, eps3, eps4, eps5, overlay_type)
    VALUES (?, ?, ?, ?, ?, ?, ?)
    ''', (event_id, seed, eps2, eps3, eps4, eps5, overlay_type))
    connection.commit()

def insert_hydro(connection, event_id, hydro_type, dimensions):
    cursor = connection.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO hydro (event_id, hydro_type, dimensions)
    VALUES (?, ?, ?)
    ''', (event_id, hydro_type, dimensions))
    connection.commit()

def insert_particlization(connection, event_id, seed,particlization_type, nsamples):
    cursor = connection.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO particlization (event_id, seed, particlization_type, nsamples)
    VALUES (?, ?,?, ?)
    ''', (event_id, seed, particlization_type, nsamples))
    connection.commit()

def insert_afterburner(connection, event_id, seed, afterburner_type):
    cursor = connection.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO afterburner (event_id, seed, afterburner_type)
    VALUES (?, ?, ?)
    ''', (event_id, seed, afterburner_type))
    connection.commit()

def insert_analysis(connection, event_id, analysis_type, rapidity_cut):
    cursor = connection.cursor()
    cursor.execute('''
    INSERT OR REPLACE INTO analysis (event_id, analysis_type, rapidity_cut)
    VALUES (?, ?, ?)
    ''', (event_id, analysis_type, rapidity_cut))
    connection.commit()

def update_event_centrality_estimator(connection, event_id, centrality_estimator):
    cursor = connection.cursor()
    cursor.execute('''
    UPDATE events SET centrality_estimator = ? WHERE event_id = ?
    ''', (centrality_estimator, event_id))
    connection.commit()

def update_ic_type(connection, event_id, ic_type):
    cursor = connection.cursor()
    cursor.execute('''
    UPDATE events SET ic_type = ? WHERE event_id = ?
    ''', (ic_type, event_id))
    connection.commit()

def update_overlay_type(connection, event_id, overlay_type):
    cursor = connection.cursor()
    cursor.execute('''
    UPDATE events SET overlay_type = ? WHERE event_id = ?
    ''', (overlay_type, event_id))
    connection.commit()

def update_hydro_type(connection, event_id, hydro_type):
    cursor = connection.cursor()
    cursor.execute('''
    UPDATE events SET hydro_type = ? WHERE event_id = ?
    ''', (hydro_type, event_id))
    connection.commit()

def update_particlization_type(connection, event_id, particlization_type):
    cursor = connection.cursor()
    cursor.execute('''
    UPDATE events SET particlization_type = ? WHERE event_id = ?
    ''', (particlization_type, event_id))
    connection.commit()

def update_afterburner_type(connection, event_id, afterburner_type):
    cursor = connection.cursor()
    cursor.execute('''
    UPDATE events SET afterburner_type = ? WHERE event_id = ?
    ''', (afterburner_type, event_id))
    connection.commit()

def update_analysis_type(connection, event_id, analysis_type):
    cursor = connection.cursor()
    cursor.execute('''
    UPDATE events SET analysis_type = ? WHERE event_id = ?
    ''', (analysis_type, event_id))
    connection.commit()

