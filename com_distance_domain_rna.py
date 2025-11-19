#!/usr/bin/env python3
import argparse, os, sys, glob, math, csv
from collections import defaultdict
from typing import Dict, List, Tuple, Iterable
from Bio.PDB import MMCIFParser, PDBParser

# ---------- Element masses (atomic mass units) ----------
MASS = {
    "H": 1.0079, "D": 2.0141, "C": 12.0107, "N": 14.0067, "O": 15.999,
    "P": 30.9738, "S": 32.065, "SE": 78.96, "MG": 24.305, "K": 39.0983,
    "NA": 22.9897, "CL": 35.453, "F": 18.998, "BR": 79.904, "I": 126.904,
    "ZN": 65.38, "FE": 55.845, "CA": 40.078, "MN": 54.938, "B": 10.81
}

NUC_SET = {"A","C","G","U","DA","DC","DG","DT","I","DI","URA","ADE","GUA","CYT","THY"}

# ---------- Utility ----------
def expand_paths(structure_specs: Iterable[str]) -> List[str]:
    paths = []
    for spec in structure_specs:
        matches = glob.glob(spec)
        if matches:
            paths.extend(sorted(matches))
        else:
            # if it's a direct path that doesn't expand, still include and let existence check handle it
            paths.append(spec)
    # uniq, keep order
    seen = set()
    out = []
    for p in paths:
        if p not in seen:
            out.append(p)
            seen.add(p)
    return out

def parse_domain_spec(s: str) -> Tuple[str, Dict[str, List[Tuple[int,int]]]]:
    """
    'Label=A:126-155;262-332;B:450-511' ->
      ("Label", {"A":[(126,155),(262,332)], "B":[(450,511)]})

    Rules:
      - Separate label from ranges with '='
      - Within ranges, separate spans by ';'
      - Each span is CHAIN:LO-HI or CHAIN:IDX (single position)
    """
    if "=" not in s:
        raise ValueError(
            f"Domain spec must look like 'Label=A:126-155;B:450-511'. Got: {s}"
        )
    label, rest = s.split("=", 1)
    label = label.strip()
    if not label:
        raise ValueError(f"Empty domain label in spec: {s}")

    ranges_by_chain: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    for block in rest.split(";"):
        block = block.strip()
        if not block:
            continue
        if ":" not in block:
            raise ValueError(f"Missing ':' in span '{block}' (spec: {s})")
        ch, span = block.split(":", 1)
        ch = ch.strip()
        span = span.strip()
        if "-" in span:
            lo, hi = span.split("-", 1)
            lo, hi = int(lo), int(hi)
        else:
            lo = hi = int(span)
        ranges_by_chain[ch].append((lo, hi))
    return label, dict(ranges_by_chain)

def parse_simple_range(r: str) -> Tuple[int,int]:
    r = r.strip()
    if "-" in r:
        a,b = r.split("-", 1)
        return int(a), int(b)
    i = int(r)
    return i, i

def load_structure(path: str):
    _, ext = os.path.splitext(path.lower())
    if ext in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure("complex", path)

def atom_mass(atom) -> float:
    el = atom.element.upper().strip() if atom.element else ""
    if not el:
        el = atom.get_name()[0].upper()
    return MASS.get(el, 12.0107)  # default to carbon

def is_nucleotide_res(res) -> bool:
    return res.get_resname().upper().strip() in NUC_SET

def select_domain_atoms(structure, ranges_by_chain: Dict[str, List[Tuple[int,int]]]):
    atoms = []
    for model in structure:
        for chain in model:
            ch = chain.id
            if ch not in ranges_by_chain:
                continue
            for res in chain:
                hetflag, resseq, icode = res.id
                if hetflag != " ":
                    continue
                if any(lo <= resseq <= hi for (lo,hi) in ranges_by_chain[ch]):
                    for at in res:
                        if (at.element or "").upper().strip() == "H":
                            continue
                        atoms.append(at)
    return atoms

def select_rna_segment_atoms(structure, chain_id: str, res_range: Tuple[int,int]):
    lo, hi = res_range
    atoms = []
    for model in structure:
        chain = model[chain_id] if chain_id in [c.id for c in model] else None
        if chain is None:
            continue
        for res in chain:
            hetflag, resseq, icode = res.id
            if hetflag != " ":
                continue
            if not (lo <= resseq <= hi):
                continue
            # Accept modified bases too; just keep heavy atoms
            for at in res:
                if (at.element or "").upper().strip() == "H":
                    continue
                atoms.append(at)
    return atoms

def center_of_mass(atoms):
    if not atoms:
        return None
    msum = 0.0
    cx = cy = cz = 0.0
    for a in atoms:
        m = atom_mass(a)
        x,y,z = a.coord
        msum += m
        cx += m*x; cy += m*y; cz += m*z
    if msum == 0.0:
        # fallback to geometric center
        n = float(len(atoms))
        return (
            sum(a.coord[0] for a in atoms)/n,
            sum(a.coord[1] for a in atoms)/n,
            sum(a.coord[2] for a in atoms)/n
        )
    return (cx/msum, cy/msum, cz/msum)

def euclid(a, b) -> float:
    dx = a[0]-b[0]; dy = a[1]-b[1]; dz = a[2]-b[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser(
        description="Batch COM distance between multiple protein domains and an RNA segment across multiple structures."
    )
    ap.add_argument("--structures", nargs="+", required=True,
                    help="One or more file paths or glob patterns (e.g. 'AF_cifs/*.cif').")
    ap.add_argument("--domain", action="append", required=True,
                    help="Repeatable. Format: Label=CHAIN:LO-HI;CHAIN:LO-HI;...\n"
                         "Example: --domain NTD=A:1-90 --domain CTD=A:200-320;A:330-360")
    ap.add_argument("--rna_chain", default="B", help="RNA chain ID (default: B)")
    ap.add_argument("--rna_res", default="11-13", help='RNA residue range (default: "11-13")')
    ap.add_argument("--out_csv", default="com_distances.csv", help="Output CSV path (default: com_distances.csv)")
    args = ap.parse_args()

    files = expand_paths(args.structures)
    if not files:
        sys.exit("No input structures matched.")
    rna_range = parse_simple_range(args.rna_res)

    # Parse all domain specs
    domain_specs = []
    for dspec in args.domain:
        label, ranges_by_chain = parse_domain_spec(dspec)
        domain_specs.append((label, ranges_by_chain))

    # Prepare CSV
    cols = [
        "file", "domain_label",
        "rna_chain", "rna_res_start", "rna_res_end",
        "dom_com_x", "dom_com_y", "dom_com_z",
        "rna_com_x", "rna_com_y", "rna_com_z",
        "com_distance_A"
    ]
    out_rows = []

    for path in files:
        if not os.path.exists(path):
            print(f"[WARN] Not found: {path}", file=sys.stderr)
            continue
        try:
            structure = load_structure(path)
        except Exception as e:
            print(f"[ERROR] Failed to parse {path}: {e}", file=sys.stderr)
            continue

        # RNA atoms once per file
        rna_atoms = select_rna_segment_atoms(structure, args.rna_chain, rna_range)
        if not rna_atoms:
            print(f"[WARN] No RNA atoms found for {path} chain {args.rna_chain} {rna_range[0]}-{rna_range[1]}", file=sys.stderr)
            # still produce rows with NaNs to signal missing RNA segment
            rna_com = None
        else:
            rna_com = center_of_mass(rna_atoms)

        for label, ranges_by_chain in domain_specs:
            dom_atoms = select_domain_atoms(structure, ranges_by_chain)
            if not dom_atoms:
                print(f"[WARN] No domain atoms for '{label}' in {path}", file=sys.stderr)
                dom_com = None
                d = float("nan")
            else:
                dom_com = center_of_mass(dom_atoms)
                if (dom_com is None) or (rna_com is None):
                    d = float("nan")
                else:
                    d = euclid(dom_com, rna_com)

            out_rows.append({
                "file": os.path.basename(path),
                "domain_label": label,
                "rna_chain": args.rna_chain,
                "rna_res_start": rna_range[0],
                "rna_res_end": rna_range[1],
                "dom_com_x": (None if dom_com is None else round(dom_com[0], 3)),
                "dom_com_y": (None if dom_com is None else round(dom_com[1], 3)),
                "dom_com_z": (None if dom_com is None else round(dom_com[2], 3)),
                "rna_com_x": (None if rna_com is None else round(rna_com[0], 3)),
                "rna_com_y": (None if rna_com is None else round(rna_com[1], 3)),
                "rna_com_z": (None if rna_com is None else round(rna_com[2], 3)),
                "com_distance_A": (None if d != d else round(d, 3)),  # keep NaN as None
            })

    # Write CSV
    os.makedirs(os.path.dirname(args.out_csv) or ".", exist_ok=True)
    with open(args.out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    # Console summary
    print(f"\nWrote {len(out_rows)} rows to {args.out_csv}")
    # quick per-domain min distance summary
    by_dom = defaultdict(list)
    for r in out_rows:
        by_dom[r["domain_label"]].append((r["file"], r["com_distance_A"]))
    for label, items in by_dom.items():
        vals = [v for (_, v) in items if isinstance(v, (int,float)) and v==v]
        if vals:
            print(f"  {label}: N={len(vals)}  mean={sum(vals)/len(vals):.3f} Å  min={min(vals):.3f} Å  max={max(vals):.3f} Å")
        else:
            print(f"  {label}: no valid distances")

if __name__ == "__main__":
    main()
