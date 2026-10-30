"""Secondary structure (DSSP codes) of rank_001 structure files in a directory.
Usage: ost sec_struct.py <directory> <output.json>

Finds .pdb, .pdb.gz, .cif and .cif.gz files containing 'rank_001' in their name
within <directory> and writes {filename: secondary structure string} as JSON.
"""
import json
import os
import sys

from ost import io, mol

PDB = (".pdb", ".pdb.gz")
CIF = (".cif", ".cif.gz")


def load(path):
    if path.endswith(CIF):
        return io.LoadMMCIF(path)
    return io.LoadPDB(path)


def main(directory, out_path):
    result = {}
    for root, _, files in os.walk(directory):
        for name in sorted(files):
            # 1. Dateiendung prüfen
            if not name.endswith(PDB + CIF):
                continue

            # 2. Nur Rank 1 berücksichtigen
            if "rank_001" not in name:
                continue

            path = os.path.join(root, name)
            key = os.path.relpath(path, directory)
            try:
                ent = load(path)
                mol.alg.AssignSecStruct(ent)
                result[key] = "".join(str(r.GetSecStructure()) for r in ent.residues)
            except Exception as e:
                print(f"skipped {key}: {e}", file=sys.stderr)

    with open(out_path, "w") as f:
        json.dump(result, f, indent=1)

    print(f"{len(result)} structures written to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])