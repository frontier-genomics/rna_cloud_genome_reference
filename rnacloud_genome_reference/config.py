from configparser import ConfigParser

class Config:
    _config_instance = None

    @classmethod
    def _get_instance(cls):
        if cls._config_instance is None:
            cls._config_instance = ConfigParser()
            cls._config_instance.read('config/config.ini')
        return cls._config_instance

    @classmethod
    def get_str(cls, section: str, key: str):
        config = cls._get_instance()
        value = config.get(section, key)
        if not isinstance(value, str):
            raise TypeError(f"Section: {section}, Key: {key} is not a string")
        return value