import argparse
import pathway_converter
import constants
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--download', action='store_true',
                        help='''Download the latest database dumps.
                        If omitted, database dumps must already be present in 'raw' folder.''')
    return parser.parse_args()


def main():
    """
    The importer pipeline consists in two phases:
    I. Convert Pathway Data from different databases from XML to JSON (Optionally download latest database dumps)
        Right now, Wikipathways and KEGG are used
    II. Insert the Pathways and Pathway-Gene Mappings into ProteomicsDB
    """

    parsed_arguments = parse_arguments()

    for organism in constants.ORGANISMS:
        print('Importing data for organism: {}'.format(organism['identifier']))
        raw_dir = Path('..', 'raw', str(organism['taxcode']))
        kegg_raw_dir = raw_dir.joinpath('kegg')
        wikipathways_raw_dir = raw_dir.joinpath('wikipathways')
        # Optionally download the latest database dumps
        if parsed_arguments.download:
            pathway_converter.download_kegg_pathways(kegg_raw_dir, organism['identifier'])
            pathway_converter.download_wikipathways_pathways(wikipathways_raw_dir, organism['identifier'])

        # Convert the various XML formats to a unified JSON representation
        json_dir = Path('..', 'json', str(organism['taxcode']))
        kegg_json_dir = json_dir.joinpath('kegg')
        wikipathways_json_dir = json_dir.joinpath('wikipathways')
        pathway_converter.process_kegg_pathways(kegg_raw_dir, kegg_json_dir, organism['identifier'])
        pathway_converter.process_wikipathways_pathways(wikipathways_raw_dir, wikipathways_json_dir,
                                                        organism)


if __name__ == '__main__':
    main()
