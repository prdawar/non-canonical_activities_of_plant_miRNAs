#!/usr/bin/env python3
import os, sys, glob, csv
from typing import Iterable, Tuple, List, Dict
from Bio.PDB import MMCIFParser, Superimposer

# --------- CONFIG (edit if needed) ----------
CIF_DIR = "/Users/jaye143/Desktop/Protein_RMSD/DOMAINS/Non-canonical/AF_cifs"
CHAIN_IDS = {"B", "C"}             # RNA chains to use
OUT_CSV   = "rna_chainBC_rmsd.csv" # output summary
ALLOW_MISSING = False              # set True to skip files missing any chain

# --------- Helpers ----------
RNA_ATOM_KEYS = {"P", "OP1", "OP2", "OP3", "O5'", "O3'", "C1'", "C2'", "C3'", "C4'", "C5'", "N1", "N9"}
RNA_NAME_HINT = {"A","U","G","C","DA","DT","DG","DC","URA","PSU","I","DI","OMC","5MC","H2U","M2G","M1G","M7G"}

def is_heavy(atom) -> bool:
    # Biopython sets atom.element when available; fall back to name check
    el = getattr(atom, "element", None)
    if el:
        return el.upper() != "H"
    return not atom.get_id().startswith("H")

def is_rna_residue(res) -> bool:
    # Heuristic: resname indicates nucleotide *or* has typical RNA atoms
    rn = res.get_resname().strip().upper()
    if rn in RNA_NAME_HINT: 
        return True
    atom_names = {a.get_id().upper() for a in res.get_atoms()}
    if "P" in atom_names or ("C1'" in atom_names) or ("N1" in atom_names) or ("N9" in atom_names):
        return True
    return False

def chain_residues(structure, chain_ids: Iterable[str]):
    """Return residues from specified chains, preserving order."""
    model = structure[0]  # first model
    residues = []
    for cid in chain_ids:
        if cid not in model:
            if ALLOW_MISSING:
                continue
            else:
                raise KeyError(f"Chain '{cid}' not found.")
        ch = model[cid]
        for res in ch.get_residues():
            # skip waters and hetero that aren’t standard residues unless recognized as RNA
            hetflag, resseq, icode = res.id
            if hetflag != " " and not is_rna_residue(res):
                continue
            residues.append((cid, res))
    return residues

def build_rna_atom_list(structure, chain_ids: Iterable[str]) -> List[Tuple[Tuple[str, int, str, str], object]]:
    """
    Collect ordered heavy atoms from RNA residues in selected chains.
    Returns list of ((chain, res_idx_seq, resname, atomname), Atom)
    Ordering is by (chain order in CHAIN_IDS) then residue order then atom name.
    """
    # stable chain order: preserve the order given by CHAIN_IDS
    ordered_chains = [cid for cid in CHAIN_IDS if cid in structure[0]]
    items = []
    for cid in ordered_chains:
        ch = structure[0][cid]
        for res in ch.get_residues():
            if not is_rna_residue(res): 
                continue
            # residue key by running order (resseq may differ across files)
            # To keep things robust, we rely on list position; we’ll match by relative order later
            # But we need something stable for atom pairing; we’ll only use atom names intersection.
            resname = res.get_resname().strip()
            # gather heavy atoms only
            heavy_atoms = [a for a in res.get_atoms() if is_heavy(a)]
            # sort atoms by name for deterministic order
            for a in sorted(heavy_atoms, key=lambda x: x.get_id()):
                items.append(((cid, -1, resname, a.get_id().upper()), a))
    return items

def build_residue_blocks(structure, chain_ids: Iterable[str]) -> List[Tuple[str, object]]:
    """Return list of RNA residues (chain-tagged) in order, for block-wise matching."""
    blocks = []
    for cid in [c for c in CHAIN_IDS if c in structure[0]]:
        ch = structure[0][cid]
        for res in ch.get_residues():
            if is_rna_residue(res):
                blocks.append((cid, res))
    return blocks

def pair_atoms_by_order_and_name(ref_blocks, mob_blocks) -> Tuple[List[object], List[object]]:
    """
    Pair atoms between reference and mobile by:
      1) walking residues in order (min common length),
      2) intersecting atom-name sets within each residue,
      3) appending atoms in sorted(atom name) order.
    """
    n = min(len(ref_blocks), len(mob_blocks))
    ref_atoms, mob_atoms = [], []
    for i in range(n):
        _, rres = ref_blocks[i]
        _, mres = mob_blocks[i]
        rnames = {a.get_id().upper(): a for a in rres.get_atoms() if is_heavy(a)}
        mnames = {a.get_id().upper(): a for a in mres.get_atoms() if is_heavy(a)}
        common = sorted(set(rnames.keys()) & set(mnames.keys()))
        for an in common:
            ref_atoms.append(rnames[an])
            mob_atoms.append(mnames[an])
    return ref_atoms, mob_atoms

def load_structure(path: str):
    parser = MMCIFParser(QUIET=True)
    struct_id = os.path.basename(path)
    return parser.get_structure(struct_id, path)

def compute_rna_rmsd(ref_path: str, target_path: str) -> Tuple[float, int]:
    ref = load_structure(ref_path)
    mob = load_structure(target_path)
    # Build block-wise residue order (RNA only)
    ref_blocks = build_residue_blocks(ref, CHAIN_IDS)
    mob_blocks = build_residue_blocks(mob, CHAIN_IDS)

    if len(ref_blocks) == 0 or len(mob_blocks) == 0:
        raise ValueError("No RNA residues found in one or both structures for chains " + ",".join(CHAIN_IDS))

    ref_atoms, mob_atoms = pair_atoms_by_order_and_name(ref_blocks, mob_blocks)
    if len(ref_atoms) < 3:
        raise ValueError(f"Not enough common heavy atoms to align (found {len(ref_atoms)}).")

    # Align mobile to reference using RNA atoms only
    si = Superimposer()
    si.set_atoms(ref_atoms, mob_atoms)  # pass Atom objects directly
    si.apply(mob_atoms)                 # apply transform to mobile atoms
    return float(si.rms), len(ref_atoms)

def main():
    cif_files = sorted(glob.glob(os.path.join(CIF_DIR, "*.cif")))
    if len(cif_files) < 2:
        print(f"Need at least 2 CIFs in {CIF_DIR}")
        sys.exit(1)

    ref_path = cif_files[0]
    print(f"[INFO] Reference: {os.path.basename(ref_path)} (RNA chains: {','.join(sorted(CHAIN_IDS))})")
    rows = []
    for path in cif_files:
        try:
            rmsd, natoms = compute_rna_rmsd(ref_path, path)
            rows.append({"file": os.path.basename(path), "rmsd_A": rmsd, "n_atoms": natoms})
            print(f"{os.path.basename(path):50s}  RMSD = {rmsd:7.3f} Å  (atoms: {natoms})")
        except Exception as e:
            rows.append({"file": os.path.basename(path), "rmsd_A": None, "n_atoms": 0, "error": str(e)})
            print(f"{os.path.basename(path):50s}  ERROR: {e}")

    out_csv = os.path.join(CIF_DIR, OUT_CSV)
    fieldnames = ["file", "rmsd_A", "n_atoms", "error"]
    # Ensure error column exists for all rows
    for r in rows:
        if "error" not in r:
            r["error"] = ""
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n[OK] Wrote summary: {out_csv}")

if __name__ == "__main__":
    main()
