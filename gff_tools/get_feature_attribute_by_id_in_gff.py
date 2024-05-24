import argparse
import fileinput
import itertools
import os
import sys

from pathlib import Path
from typing import List

import gffutils
import gffutils.exceptions


missing_ids: int = 0


def main(gff_file: Path, id_list: List[Path], attribute_keys: List[str], ignore_id_error: bool, recursive: bool):
    # Create GFF db
    gff_db = gffutils.create_db(str(gff_file), ':memory:', merge_strategy="create_unique")

    with fileinput.input(id_list) as ids:
        # For each ID in id_list search for the attributes and instantly write them to STDOUT
        for id_line in ids:
            id_line = id_line.strip()
            attributes: List[str] = []
            # gffutils databases seem to have no "in" funtion for keys/IDs, sp try/escept strucutre is necessary
            try:
                id_feature: gffutils.Feature = gff_db[id_line]
                for attribute_key in attribute_keys:
                    # Search also in parents and children
                    if recursive:
                        attribute_results: List[str] = []
                        # Search in ID feature, all parents and all children
                        feature: gffutils.Feature
                        for feature in itertools.chain([id_feature], gff_db.parents(id_line), gff_db.children(id_line)):
                            if attribute_key in feature.attributes:
                                attribute_results.append(f'{feature.featuretype}:{feature.attributes[attribute_key][0]}')
                        if len(attribute_results) >= 1:
                            attributes.append(";".join(attribute_results))
                        else:
                            attributes.append("-")
                    
                    # Search only in ID feature
                    else:
                        attributes.append(id_feature.attributes.get(attribute_key, ["-"])[0])
            # ID was not fount in db
            except gffutils.exceptions.FeatureNotFoundError:
                if ignore_id_error:
                    attributes = ["-"] * len(attribute_keys)
                else:
                    print('\033[91m' + f"ERROR:\nUnknown ID: {id_line}" + '\033[0m', file=sys.stderr)
                    sys.exit(1)
            
            print("{}\t{}".format(id_line, "\t".join(attributes)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog=os.path.basename(__file__),
                    description='Get GFF feature attribute entries for ID list')
    parser.add_argument('-g', '--gff_file', type=Path, help="Input GFF file", required=True)
    parser.add_argument('-a', '--attributes', help="Requested attributes.", nargs='+', required=True)
    parser.add_argument('-i', '--id_list', type=Path, help="Input IDs file. STDIN if not given.", nargs='+')
    parser.add_argument('--ignore-id-error', dest='ignore_id_error', help="Ignore unknown IDs.", action='store_true')
    parser.add_argument('--recursive', help="Search for attributes also in parents and children", action='store_true')
    args = parser.parse_args()
    
    if args.id_list is None:
        args.id_list = []
    
    main(args.gff_file, args.id_list, args.attributes, args.ignore_id_error, args.recursive)