"""Configuration and environment helpers extracted from shared.utils."""

import configparser as ConfigParser
import os


def configMap(cfgFile):
    """Reads a configuration file and returns a dictionary of all sections, options and values."""
    result = {}
    config = ConfigParser.ConfigParser()
    # Activates case-sensitivity:
    config.optionxform = str

    try:
        config.read(cfgFile)
    except ConfigParser.ParsingError:
        raise ValueError("Invalid configuration file %s" % cfgFile)

    for sect in config.sections():
        result[sect] = {}
        options = config.options(sect)
        for opt in options:
            try:
                result[sect][opt] = config.get(sect, opt)
            except Exception:
                result[sect][opt] = None
    return result


def getEnvironmentValue(name, default=None):
    """Returns the value for the given environment variable name, if found;
    otherwise returns default value."""
    try:
        return os.environ[name]
    except KeyError:
        return default
