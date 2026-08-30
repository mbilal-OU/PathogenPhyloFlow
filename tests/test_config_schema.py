from pathlib import Path

import jsonschema
import yaml


def test_default_config_matches_schema():
    root = Path(__file__).resolve().parents[1]
    config = yaml.safe_load((root / "config" / "config.yaml").read_text(encoding="utf-8"))
    schema = yaml.safe_load((root / "config" / "config.schema.yaml").read_text(encoding="utf-8"))
    jsonschema.validate(config, schema)
