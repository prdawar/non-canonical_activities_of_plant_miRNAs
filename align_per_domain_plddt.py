#!/usr/bin/env python3
# align_per_domain_plddt.py
import os, sys, glob, csv, argparse, math, copy
from typing import Dict, List, Tuple, Optional
from Bio.PDB import MMCIFParser, PDBParser, Superimposer, is_aa, Select, PDBIO

# ----------------------------- Parsing helpers -----------------------------

def load_structure(path: str):
    ext = os.path.splitext(path)[1].lower()
    if ext in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True).get_structure(os.path.basename(path), path)
    if ext in (".pdb", ".ent"):
        return PDBParser(QUIET=True).get_structure(os.path.basename(path), path)
    # auto-try CIF then PDB
    try:
        return MMCIFParser(QUIET=True).get_structure(os.path.basename(path), path)
    except Exception:
        return PDBParser(QUIET=True).get_structure(os.path.basename(path), path)

def parse_domain_spec(domains: str) -> Dict[str, Dict[Optional[str], List[Tuple[int,int]]]]:
    """
    Parse a domain specification string into:
      { domain_name: { chain_id_or_None: [(start,end), ...], ... }, ... }

    Syntax (flexible, commas separate domains, semicolons separate ranges):
      "RBD=A:50-120;130-160, CTD=B:10-80"
      "Core=50-200"  (applies to all chains)
      "D1=A:10-40;60-90,B:15-25"
    """
    out: Dict[str, Dict[Optional[str], List[Tuple[int,int]]]] = {}
    if not domains.strip():
        return out
    for chunk in domains.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise ValueError(f"Domain chunk missing '=': {chunk}")
        name, spec = chunk.split("=", 1)
        name = name.strip()
        spec = spec.strip()
        chain_map: Dict[Optional[str], List[Tuple[int,int]]] = {}
        if ":" not in spec and ";" not in spec and "-" in spec:
            # single range for all chains
            a, b = spec.split("-")
            chain_map[None] = [(int(a), int(b))]
        else:
            # e.g., "A:50-120;130-160" or "A:10-40;60-90 B:.."
            # Allow multiple chain groups separated by spaces or semicolons? We'll split on spaces and handle commas before.
            for grp in spec.split(","):
                grp = grp.strip()
                if not grp:
                    continue
                # Allow multiple ranges per chain separated by ';'
                if ":" not in grp:
                    # treat as all-chain ranges
                    for seg in grp.split(";"):
                        seg = seg.strip()
                        if seg:
                            a, b = seg.split("-")
                            chain_map.setdefault(None, []).append((int(a), int(b)))
                else:
                    chain, segs = grp.split(":", 1)
                    chain = chain.strip()
                    for seg in segs.split(";"):
                        seg = seg.strip()
                        if seg:
                            a, b = seg.split("-")
                            chain_map.setdefault(chain, []).append((int(a), int(b)))
        out[name] = chain_map
    return out

def res_in_ranges(resseq: int, ranges: List[Tuple[int,int]]) -> bool:
    if not ranges:  # empty means unrestricted
        return True
    for a, b in ranges:
        if a <= resseq <= b:
            return True
    return False

# ----------------------------- Atom collection -----------------------------

def collect_atoms_for_domain(structure,
                             ranges_by_chain: Dict[Optional[str], List[Tuple[int,int]]],
                             atom_mode: str,
                             bmin: Optional[float] = None,
                             bfilter_on: bool = False):
    """
    Build atom list keyed by (chain, resseq, icode, atom_name or 'CA').
    For AF set bfilter_on=True to enforce pLDDT (B-factor) >= bmin.
    """
    atom_mode = atom_mode.upper()
    atoms = {}
    for model in structure:
        for chain in model:
            chain_id = chain.id if chain.id is not None else ""
            # domain ranges for this chain or 'None'
            if any(k is not None for k in ranges_by_chain.keys()):
                ranges = ranges_by_chain.get(chain_id, [])
            else:
                ranges = ranges_by_chain.get(None, [])
            for res in chain:
                het, resseq, icode = res.id
                if het.strip():
                    continue
                if not is_aa(res, standard=True):
                    continue
                if not res_in_ranges(int(resseq), ranges):
                    continue
                if atom_mode == "CA":
                    if "CA" in res:
                        a = res["CA"]
                        if bfilter_on and bmin is not None and a.bfactor < bmin:
                            continue
                        key = (chain_id, int(resseq), (icode or "").strip(), "CA")
                        atoms[key] = a
                else:
                    for a in res.get_atoms():
                        if a.element == "H":
                            continue
                        if bfilter_on and bmin is not None and a.bfactor < bmin:
                            continue
                        key = (chain_id, int(resseq), (icode or "").strip(), a.get_name().strip())
                        atoms[key] = a
    return atoms

# ----------------------------- Alignment & RMSD -----------------------------

def align_rmsd(ref_map, af_map) -> Tuple[float, int, Superimposer, List[Tuple]]:
    keys = sorted(set(ref_map.keys()) & set(af_map.keys()))
    if not keys:
        return float("nan"), 0, None, []
    ref_atoms = [ref_map[k] for k in keys]
    af_atoms  = [af_map[k] for k in keys]
    sup = Superimposer()
    sup.set_atoms(ref_atoms, af_atoms)  # computes transform AF->Ref
    # sup.rms is reported; we also return count
    return float(sup.rms), len(keys), sup, keys

# ----------------------------- Selectors for saving -----------------------------

class DomainSelector(Select):
    def __init__(self, ranges_by_chain, atom_mode, bmin, bfilter_on):
        self.ranges_by_chain = ranges_by_chain
        self.atom_mode = atom_mode.upper()
        self.bmin = bmin
        self.bfilter_on = bfilter_on

    def accept_residue(self, residue):
        het, resseq, icode = residue.id
        if het.strip():
            return 0
        chain = residue.get_parent()
        chain_id = chain.id if chain is not None else ""
        if any(k is not None for k in self.ranges_by_chain.keys()):
            ranges = self.ranges_by_chain.get(chain_id, [])
        else:
            ranges = self.ranges_by_chain.get(None, [])
        return 1 if res_in_ranges(int(resseq), ranges) else 0

    def accept_atom(self, atom):
        if self.atom_mode == "CA" and atom.get_name().strip() != "CA":
            return 0
        if atom.element == "H":
            return 0 if self.atom_mode == "HEAVY" else 1
        if self.bfilter_on and self.bmin is not None and atom.bfactor < self.bmin:
            return 0
        return 1

# ----------------------------- Per-file processing -----------------------------

def process_file(ref_path: str,
                 af_path: str,
                 domains: Dict[str, Dict[Optional[str], List[Tuple[int,int]]]],
                 bmin: float,
                 atom_mode: str,
                 outdir: str,
                 save_pdb: bool):

    ref = load_structure(ref_path)
    af  = load_structure(af_path)

    # Collect per-domain RMSDs
    per_domain = []
    all_ref_keys = []
    all_af_keys  = []
    ref_maps = {}
    af_maps  = {}

    # Build maps per domain
    for dname, chain_ranges in domains.items():
        ref_map = collect_atoms_for_domain(ref, chain_ranges, atom_mode, bmin=None, bfilter_on=False)
        af_map  = collect_atoms_for_domain(af,  chain_ranges, atom_mode, bmin=bmin, bfilter_on=True)
        ref_maps[dname] = ref_map
        af_maps[dname]  = af_map

        rms, n, sup, keys = align_rmsd(ref_map, af_map)
        per_domain.append((dname, rms, n))

        # Save aligned PDB for this domain, if requested
        if save_pdb and n > 0:
            # apply transform to a *copy* so we don't disturb other domains
            af_copy = copy.deepcopy(af)
            # Collect same atoms in the copy by building a map again (coordinates correspond)
            af_map_copy  = collect_atoms_for_domain(af_copy, chain_ranges, atom_mode, bmin=bmin, bfilter_on=True)
            # Apply transform to those atoms
            sup.apply([af_map_copy[k] for k in keys])

            io = PDBIO()
            io.set_structure(af_copy)
            sel = DomainSelector(chain_ranges, atom_mode, bmin, True)
            base = os.path.splitext(os.path.basename(af_path))[0]
            out_path = os.path.join(outdir, f"{base}__{dname}_aligned.pdb")
            io.save(out_path, sel)

        # For the overall RMSD, accumulate keys (union over all domains)
        all_ref_keys.extend(list(ref_map.keys()))
        all_af_keys.extend(list(af_map.keys()))

    # Compute overall RMSD on the union of domain selections (intersection with AF present)
    overall_ref = {}
    overall_af  = {}
    for dname in domains.keys():
        overall_ref.update(ref_maps[dname])
        overall_af.update(af_maps[dname])

    ov_rms, ov_n, _, _ = align_rmsd(overall_ref, overall_af)

    return per_domain, (ov_rms, ov_n)

# ----------------------------- CLI -----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Per-domain alignment (AF→Ref) with pLDDT filter (B ≥ bmin) on AF and explicit domain residue ranges."
    )
    ap.add_argument("--ref", required=True, help="Reference structure (mmCIF or PDB)")
    ap.add_argument("--af_glob", required=True, help="Glob for AlphaFold mmCIF/PDB files, e.g. 'AF/*.cif'")
    ap.add_argument("--domains", required=True,
                    help=("Named domains with ranges, e.g. "
                          "'RBD=A:50-120;130-160, CTD=B:10-80' or 'Core=50-200'"))
    ap.add_argument("--bmin", type=float, default=60.0, help="pLDDT threshold on AF B-factor (keep B ≥ bmin)")
    ap.add_argument("--atoms", choices=["CA","heavy"], default="CA", help="Fit on C-alpha or all heavy atoms")
    ap.add_argument("--outdir", default="per_domain_out", help="Output directory")
    ap.add_argument("--save_pdb", action="store_true", help="Save per-domain aligned AF PDBs")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    dom_spec = parse_domain_spec(args.domains)

    af_files = sorted(glob.glob(args.af_glob))
    if not af_files:
        print(f"[WARN] No AF files matched: {args.af_glob}")
        sys.exit(0)

    csv_path = os.path.join(args.outdir, "rmsd_per_domain.csv")
    with open(csv_path, "w", newline="") as f:
        fieldnames = ["af_file","atoms","bmin","domain","range_spec","n_atoms","rmsd"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for af in af_files:
            try:
                per_dom, overall = process_file(
                    ref_path=args.ref,
                    af_path=af,
                    domains=dom_spec,
                    bmin=args.bmin,
                    atom_mode=args.atoms,
                    outdir=args.outdir,
                    save_pdb=args.save_pdb
                )
                # per-domain rows
                for dname, rms, n in per_dom:
                    writer.writerow({
                        "af_file": af,
                        "atoms": args.atoms,
                        "bmin": args.bmin,
                        "domain": dname,
                        "range_spec": args.domains,
                        "n_atoms": n,
                        "rmsd": f"{rms:.6f}" if n > 0 else "nan",
                    })
                    print(f"{os.path.basename(af)}  [{dname}]  RMSD={rms:.3f} Å  (n={n})")

                # overall row
                ov_rms, ov_n = overall
                writer.writerow({
                    "af_file": af,
                    "atoms": args.atoms,
                    "bmin": args.bmin,
                    "domain": "ALL",
                    "range_spec": args.domains,
                    "n_atoms": ov_n,
                    "rmsd": f"{ov_rms:.6f}" if ov_n > 0 else "nan",
                })
                print(f"{os.path.basename(af)}  [ALL]     RMSD={ov_rms:.3f} Å  (n={ov_n})")

            except Exception as e:
                print(f"[ERROR] {af}: {e}")
                # write an error line for visibility
                writer.writerow({
                    "af_file": af, "atoms": args.atoms, "bmin": args.bmin,
                    "domain": "ERROR", "range_spec": args.domains,
                    "n_atoms": 0, "rmsd": "nan"
                })

    print(f"\n[DONE] Wrote: {csv_path}")

if __name__ == "__main__":
    main()
