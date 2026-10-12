# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, json, re

# Configuration files live in the base psst directory; my_config.json (not tracked
# by git) is used if it exists, otherwise default_config.json is used
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_CONFIG_FILE = os.path.join(BASE_DIR, "default_config.json")
USER_CONFIG_FILE = os.path.join(BASE_DIR, "my_config.json")

# Sections whose keys are chosen by the user rather than fixed by default_config.json
FREEFORM_SECTIONS = ["rpid_shortcuts"]

RPID_REGEX = re.compile(r"^[A-Z]{6}\d{4}$")

class ConfigError(Exception):
    pass

def read_config(configFile):
    '''
    Parameters:
        configFile -- a string indicating the location of a JSON configuration file.
    Returns:
        config -- a dict of the parsed configuration.
    '''
    try:
        with open(configFile, "r") as fileIn:
            config = json.load(fileIn)
    except json.JSONDecodeError as e:
        raise ConfigError(f"Configuration file '{configFile}' is not valid JSON; {e}")
    if not isinstance(config, dict):
        raise ConfigError(f"Configuration file '{configFile}' must contain a JSON object")
    return config

def validate_config(config, configFile, schema):
    '''
    Validates a configuration against the structure of default_config.json, which acts
    as the schema; the same sections and keys must be present, with no extras, except
    for FREEFORM_SECTIONS which may contain any keys.

    Parameters:
        config -- a dict of the parsed configuration to validate.
        configFile -- a string indicating the location of the configuration file.
        schema -- a dict of the parsed default_config.json.
    Returns:
        config -- the validated dict, with RPID shortcut keys and values uppercased.
    '''
    # Validate sections
    missing = [ s for s in schema if s not in config ]
    unknown = [ s for s in config if s not in schema ]
    if missing != [] or unknown != []:
        raise ConfigError(f"Configuration file '{configFile}' sections do not match " +
                          f"'{DEFAULT_CONFIG_FILE}'; missing = {missing}, unknown = {unknown}")

    # Validate keys and values within each section
    for section, values in config.items():
        if not isinstance(values, dict):
            raise ConfigError(f"Section '{section}' in '{configFile}' must be a JSON object")
        if section not in FREEFORM_SECTIONS:
            missing = [ k for k in schema[section] if k not in values ]
            unknown = [ k for k in values if k not in schema[section] ]
            if missing != [] or unknown != []:
                raise ConfigError(f"Section '{section}' in '{configFile}' keys do not match " +
                                  f"'{DEFAULT_CONFIG_FILE}'; missing = {missing}, unknown = {unknown}")
        for key, value in values.items():
            if not isinstance(value, str) or value == "":
                raise ConfigError(f"'{section}' > '{key}' in '{configFile}' must be a non-empty string")

    # Validate RPID shortcuts
    shortcuts = {}
    for key, value in config["rpid_shortcuts"].items():
        if RPID_REGEX.match(value.upper()) == None:
            raise ConfigError(f"RPID shortcut '{key}' in '{configFile}' must map to 6 letters " +
                              f"followed by 4 numbers, not '{value}'")
        shortcuts[key.upper()] = value.upper()
    config["rpid_shortcuts"] = shortcuts

    return config

def load_config():
    '''
    Returns:
        config -- a dict of the validated configuration from my_config.json if it
                  exists, otherwise from default_config.json.
        configFile -- a string of the configuration file location that was used.
    '''
    if not os.path.isfile(DEFAULT_CONFIG_FILE):
        raise ConfigError(f"Default configuration file '{DEFAULT_CONFIG_FILE}' does not exist!")
    schema = read_config(DEFAULT_CONFIG_FILE)

    configFile = USER_CONFIG_FILE if os.path.isfile(USER_CONFIG_FILE) else DEFAULT_CONFIG_FILE
    config = validate_config(read_config(configFile), configFile, schema)
    return config, configFile

CONFIG, CONFIG_FILE = load_config()
