import argparse
from pathlib import Path
import yaml

def main(args):
    filepaths = {
        'dirs': {},
        'subdirs': {
            'results': [
                "processed_expansions",
                "raw_expansions",
            ],
            'artifacts': [
                "cofactors",
                "rules",
                "starters_targets",
                "imgs",
                "operator_mapping",
            ],
            'data': [
                "sprhea",
            ],
        },
    }

    # Add user selected top-level direcotries
    filepaths['dirs']['artifacts'] = args.artifacts
    filepaths['dirs']['data'] = args.data
    filepaths['dirs']['results'] = args.results

    # Make sub directories for directories external to project directory
    for sub in filepaths['subdirs']['data']:
        (Path(args.data) / sub).mkdir(parents=True, exist_ok=False)

    for sub in filepaths['subdirs']['results']:
        (Path(args.results) / sub).mkdir(parents=True, exist_ok=False)

    # Write the filepaths to config.yaml
    config_path = Path(__file__).parent.parent / "config.yaml"
    with open(config_path, 'w') as f:
        yaml.dump(filepaths, f, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Configs initializer.")
    parser.add_argument('artifacts', type=str, help='Filepath to artifacts directory')
    parser.add_argument('data', type=str, help='Filepath to data directory')
    parser.add_argument('results', type=str, help='Filepath to results directory')

    args = parser.parse_args()
    
    main(args)