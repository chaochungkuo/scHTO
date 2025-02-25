import os
import yaml
import logging

logger = logging.getLogger(__name__)


class ConfigValidationError(Exception):
    """Custom exception for configuration validation errors."""
    pass


def load_config(config_path: str) -> dict:
    """Load and validate the YAML configuration file."""
    with open(config_path, 'r') as file:
        config = yaml.safe_load(file)
    validate_config(config)
    return config


def validate_config(config: dict):
    """Validate that the configuration dictionary contains all required sections and keys."""
    required_sections = [
        "libraries_with_HTOs",
        "positions",
        "libraries_to_be_demultiplexed",
        "HTO_sequences",
        "expected_cell_number"
    ]
    for section in required_sections:
        if section not in config:
            logger.error("Missing required section: %s", section)
            raise ConfigValidationError(f"Missing required section: {section}")
        else:
            logger.debug("Found required section: %s", section)

    # Validate libraries sections
    for lib_section in ["libraries_with_HTOs", "libraries_to_be_demultiplexed"]:
        for entry in config.get(lib_section, []):
            for key in ["htolib_name", "R12", "path"]:
                if key not in entry:
                    logger.error("Missing key '%s' in %s entry: %s", key, lib_section, entry)
                    raise ConfigValidationError(f"Missing key '{key}' in {lib_section} entry: {entry}")
            if not os.path.exists(entry["path"]):
                logger.error("File does not exist: %s", entry["path"])
                raise ConfigValidationError(f"File does not exist: {entry['path']}")

    # Validate positions: check that each position entry has the required keys and values.
    position_keys = [
        "cell_barcode_R12", "cell_barcode_start", "cell_barcode_end",
        "umi_R12", "umi_start", "umi_end",
        "hto_R12", "hto_start", "hto_end"
    ]
    for key in position_keys:
        if key not in config["positions"]:
            logger.error("Missing key '%s' in positions section", key)
            raise ConfigValidationError(f"Missing key '{key}' in positions section")

    logger.info("Configuration file validated successfully.")