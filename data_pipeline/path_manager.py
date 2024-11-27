import os

class DatabasePathManager:
    _database_path = None  # Class-level variable to store the path

    @staticmethod
    def load_path_from_config():
        """Load the path from config.py."""
        config_path = "config.py"
        if not os.path.exists(config_path):
            # Create config.py with an empty path
            with open(config_path, "w") as f:
                f.write("database_path = None\n")
            DatabasePathManager._database_path = None
        else:
            # Load the path from config.py
            config_data = {}
            with open(config_path, "r") as f:
                exec(f.read(), config_data)
            DatabasePathManager._database_path = config_data.get("database_path")

        return DatabasePathManager._database_path

    @staticmethod
    def save_path_to_config(path):
        """Save the database path to config.py."""
        with open("config.py", "w") as f:
            f.write(f"database_path = '{path}'\n")
        DatabasePathManager._database_path = path

    @staticmethod
    def get_database_path():
        """Retrieve the current database path."""
        return DatabasePathManager._database_path

    @staticmethod
    def set_database_path(path):
        """Set the database path and save to config."""
        if os.path.exists(path) and os.path.isdir(path):
            DatabasePathManager.save_path_to_config(path)
        else:
            raise ValueError("Invalid path. Please ensure the path exists and is a directory.")