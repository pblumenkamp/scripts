import argparse
import sys

from pathlib import Path
from typing import Dict, Set

import gffutils
from memory_profiler import profile


missing_ids: int = 0

@profile
def main(gff_file: Path):
    used_ids: Set[str] = set()
    duplicate_ids_counter: Dict[str, int] = {}
    with open(gff_file) as gff:
        #Keep all comments at the beginning of the GFF
        for line in gff:
            if line.startswith("#"):
                print(line, end="")
            else:
                break
    
    # Create GFF db
    gff_db = gffutils.create_db(str(gff_file), ':memory:', merge_strategy="create_unique", transform=add_id)
    
    # Find all ID duplicates
    for entry in gff_db.all_features():
        entry_id: str = entry.attributes["ID"][0]
        if entry_id in used_ids:
            duplicate_ids_counter[entry_id] = 1
        else:
            used_ids.add(entry_id)
    del used_ids
    
    # Write GFF and add numbers to ID duplicates
    for entry in gff_db.all_features():
        entry_id = entry.attributes["ID"][0]
        if entry_id in duplicate_ids_counter:
            entry.attributes["ID"] = f'{entry_id}_{duplicate_ids_counter[entry_id]}'
            duplicate_ids_counter[entry_id] += 1
        #print(entry)
    
    print(f'{missing_ids} IDs added.', file=sys.stderr)


def add_id(x: gffutils.Feature) -> gffutils.Feature:
    if "ID" not in x.attributes:
        global missing_ids 
        missing_ids += 1
        if "Parent" in x.attributes:
            x.attributes["ID"] = f'{x.attributes["Parent"][0]}_{x.featuretype}'
        else:
            x.attributes["ID"] = f'{x.featuretype}'
    return x


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='add_missing_ids_in_gff.py',
                    description='Add missing IDs to GFF features (IDs are added at the end of attributes). Output: STDOUT')
    parser.add_argument('gff_file', type=Path, help="Input GFF file")
    args = parser.parse_args()
    main(args.gff_file)