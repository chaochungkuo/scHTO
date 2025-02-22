import os
import tempfile
import yaml
import pytest
from config_validator import load_config, ConfigValidationError


def test_missing_section():
    config_data = {"positions": []}  # Missing required sections.
    with tempfile.NamedTemporaryFile("w", delete=False) as tf:
        yaml.dump(config_data, tf)
        temp_name = tf.name
    with pytest.raises(ConfigValidationError):
        load_config(temp_name)
    os.remove(temp_name)


def test_valid_config(tmp_path):
    config_data = {
        "libraries_with_HTOs": [
            {"htolib_name": "pool1", "R12": "R1", "path": __file__},
            {"htolib_name": "pool1", "R12": "R2", "path": __file__}
        ],
        "positions": [
            {"name": "cell_barcode_start", "position": 1, "R12": "R1"},
            {"name": "cell_barcode_end", "position": 16, "R12": "R1"},
            {"name": "umi_start", "position": 17, "R12": "R1"},
            {"name": "umi_end", "position": 26, "R12": "R1"},
            {"name": "hto_start", "position": 1, "R12": "R2"},
            {"name": "hto_end", "position": 11, "R12": "R2"}
        ],
        "libraries_to_be_demultiplexed": [
            {"htolib_name": "pool1", "R12": "R1", "path": __file__},
            {"htolib_name": "pool1", "R12": "R2", "path": __file__}
        ],
        "HTO_sequences": [
            {"htolib_name": "pool1", "sample_name": "sample1", "hto_sequence": "AAA"}
        ],
        "cutoffs": {"min_umi": 5},
        "expected_cell_number": [
            {"htolib_name": "pool1", "sample_name": "sample1", "estimate_number": 30000}
        ]
    }
    config_file = tmp_path / "config.yaml"
    config_file.write_text(yaml.dump(config_data))
    config = load_config(str(config_file))
    assert "libraries_with_HTOs" in config