import yaml
from pathlib import Path

def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    for key in ["alphafold_dir", "tmp_dir", "features_dir"]:
        p = Path(cfg["paths"][key]).expanduser().resolve()
        cfg["paths"][key] = str(p)

    return cfg
