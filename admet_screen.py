#!/usr/bin/env python3
"""
ADMET Property Screening for AutoDock Vina .pdbqt ligand files
=================================================================
Screens molecules across batch1-4 folders for:
  - Lipinski's Rule of Five (drug-likeness)
  - ADMET descriptors via RDKit
  - Veber rules (oral bioavailability)
  - PAINS/Pan-assay interference filters
  - BBB permeability, P-gp substrate prediction
  - Toxicity flags (hERG, AMES)

ADMET-PASSED .pdbqt files are automatically copied to:
  <passed_dir>/all/           <- every passed file in one flat folder
  <passed_dir>/batch1/        <- per-batch sub-folders (mirrors input)
  <passed_dir>/batch2/
  ...

Usage:
    pip install rdkit pandas openpyxl
    python admet_screen.py --input_dir /path/to/3D --output results_admet.xlsx
    python admet_screen.py --input_dir /path/to/3D --output results_admet.xlsx \
                           --passed_dir /path/to/passed_molecules
"""

import os
import sys
import argparse
import glob
import shutil
import warnings
warnings.filterwarnings("ignore")

try:
    import pandas as pd
except ImportError:
    sys.exit("Install pandas:  pip install pandas openpyxl")

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, FilterCatalog
    from rdkit.Chem import QED
    from rdkit.Chem.FilterCatalog import FilterCatalogParams
except ImportError:
    sys.exit("Install RDKit:  pip install rdkit")


# ─────────────────────────────────────────────────────────────────────────────
# 1.  Parse .pdbqt -> RDKit Mol
# ─────────────────────────────────────────────────────────────────────────────

def pdbqt_to_mol(filepath):
    """Convert a .pdbqt file to an RDKit Mol object."""
    pdb_lines = []
    with open(filepath, "r") as fh:
        for line in fh:
            tag = line[:6].strip()
            if tag in ("ATOM", "HETATM", "END", "ENDMDL"):
                pdb_lines.append(line[:60].rstrip() + "\n")
    if not pdb_lines:
        return None
    pdb_block = "".join(pdb_lines) + "END\n"
    mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=True)
    if mol is None:
        return None
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    return mol


# ─────────────────────────────────────────────────────────────────────────────
# 2.  ADMET descriptor calculation
# ─────────────────────────────────────────────────────────────────────────────

_pains_params = FilterCatalogParams()
_pains_params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
_pains_catalog = FilterCatalog.FilterCatalog(_pains_params)

_brenk_params = FilterCatalogParams()
_brenk_params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
_brenk_catalog = FilterCatalog.FilterCatalog(_brenk_params)


def calc_admet(mol, name, batch):
    """Return a dict of ADMET descriptors for a single molecule."""
    d = {"Molecule": name, "Batch": batch}

    mw       = Descriptors.MolWt(mol)
    logp     = Descriptors.MolLogP(mol)
    hbd      = Lipinski.NumHDonors(mol)
    hba      = Lipinski.NumHAcceptors(mol)
    tpsa     = Descriptors.TPSA(mol)
    rb       = Lipinski.NumRotatableBonds(mol)
    rings    = rdMolDescriptors.CalcNumRings(mol)
    arom     = rdMolDescriptors.CalcNumAromaticRings(mol)
    hac      = mol.GetNumHeavyAtoms()
    fsp3     = rdMolDescriptors.CalcFractionCSP3(mol)
    mw_exact = Descriptors.ExactMolWt(mol)

    d.update({
        "MW (Da)"        : round(mw, 2),
        "Exact_MW (Da)"  : round(mw_exact, 4),
        "LogP"           : round(logp, 2),
        "HBD"            : hbd,
        "HBA"            : hba,
        "TPSA (A2)"      : round(tpsa, 2),
        "RotBonds"       : rb,
        "Rings"          : rings,
        "AromaticRings"  : arom,
        "HeavyAtomCount" : hac,
        "Fsp3"           : round(fsp3, 3),
    })

    try:
        qed_val = QED.qed(mol)
    except Exception:
        qed_val = None
    d["QED"] = round(qed_val, 3) if qed_val is not None else "N/A"

    # Lipinski Ro5
    ro5_v = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    d["Ro5_Violations"] = ro5_v
    d["Lipinski_Pass"]  = "YES" if ro5_v <= 1 else "NO"

    # Veber
    d["Veber_Pass"] = "YES" if (tpsa <= 140 and rb <= 10) else "NO"

    # Ghose
    ghose = (160 <= mw <= 480) and (-0.4 <= logp <= 5.6) and (20 <= hac <= 70)
    d["Ghose_Pass"] = "YES" if ghose else "NO"

    # Egan
    d["Egan_Pass"] = "YES" if (logp <= 5.88 and tpsa <= 131.6) else "NO"

    # Muegge
    muegge = ((200 <= mw <= 600) and (-2 <= logp <= 5) and (tpsa <= 150)
              and (rb <= 15) and (rings <= 7) and (hbd <= 5) and (hba <= 10))
    d["Muegge_Pass"] = "YES" if muegge else "NO"

    # BBB
    bbb = (mw < 450) and (1 <= logp <= 4) and (hbd <= 3) and (tpsa <= 90)
    d["BBB_Permeable"] = "YES" if bbb else "NO"

    # GI absorption
    gi_high = (tpsa <= 75) or (tpsa <= 140 and rb <= 10)
    d["GI_Absorption"] = "HIGH" if (gi_high and logp < 5) else "LOW"

    # P-gp substrate
    pgp = (mw > 400) and (tpsa > 60) and ((hbd + hba) > 8)
    d["Pgp_Substrate"] = "YES" if pgp else "NO"

    # hERG
    basic_n = any(
        atom.GetAtomicNum() == 7 and atom.GetTotalDegree() <= 3
        for atom in mol.GetAtoms()
    )
    d["hERG_Risk"] = "HIGH" if (logp > 3 and basic_n) else "LOW"

    # Ames
    ames_smarts = [
        "[$(c1ccc([N+](=O)[O-])cc1)]",
        "[$(c1ccccc1N)]",
        "[$(CC(=O)O)]",
        "[N;H2][c]",
    ]
    ames_flag = any(
        mol.HasSubstructMatch(Chem.MolFromSmarts(s))
        for s in ames_smarts
        if Chem.MolFromSmarts(s) is not None
    )
    d["Ames_Alert"] = "YES" if ames_flag else "NO"

    # PAINS
    pains_m = _pains_catalog.GetMatches(mol)
    d["PAINS_Flag"]   = "YES" if pains_m else "NO"
    d["PAINS_Alerts"] = "; ".join(m.GetDescription() for m in pains_m) if pains_m else "None"

    # Brenk
    brenk_m = _brenk_catalog.GetMatches(mol)
    d["Brenk_Flag"]   = "YES" if brenk_m else "NO"
    d["Brenk_Alerts"] = "; ".join(m.GetDescription() for m in brenk_m) if brenk_m else "None"

    # SA complexity
    sa_proxy = rings + (1 if fsp3 < 0.2 else 0) + (1 if mw > 500 else 0)
    d["SA_Complexity"] = "HIGH" if sa_proxy >= 4 else ("MED" if sa_proxy >= 2 else "LOW")

    # Lead-likeness
    d["LeadLike"] = "YES" if ((250 <= mw <= 350) and (logp <= 3.5) and (rb <= 7)) else "NO"

    # Overall
    passed = (
        d["Lipinski_Pass"] == "YES"
        and d["Veber_Pass"]   == "YES"
        and d["PAINS_Flag"]   == "NO"
        and d["Brenk_Flag"]   == "NO"
        and d["Ames_Alert"]   == "NO"
        and d["hERG_Risk"]    == "LOW"
    )
    d["Overall_ADMET_Pass"] = "PASS" if passed else "FAIL"
    return d


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Copy passed .pdbqt files to a dedicated folder
# ─────────────────────────────────────────────────────────────────────────────

def copy_passed_files(df, source_paths, passed_dir):
    """
    Copy every ADMET-PASS molecule's original .pdbqt file into:
      passed_dir/all/           <- flat pool of ALL passed files
      passed_dir/<batch>/       <- per-batch sub-folders
    """
    passed_df = df[df["Overall_ADMET_Pass"] == "PASS"]

    if passed_df.empty:
        print("\n  [!] No molecules passed ADMET filters - nothing to copy.")
        return 0

    all_dir = os.path.join(passed_dir, "all")
    os.makedirs(all_dir, exist_ok=True)

    copied  = 0
    missing = 0

    for _, row in passed_df.iterrows():
        mol_name = row["Molecule"]
        batch    = row["Batch"]
        src      = source_paths.get(mol_name)

        if not src or not os.path.isfile(src):
            print(f"  [!] Source not found for {mol_name}, skipping.")
            missing += 1
            continue

        filename = os.path.basename(src)

        # flat "all" copy
        shutil.copy2(src, os.path.join(all_dir, filename))

        # per-batch copy
        batch_out = os.path.join(passed_dir, batch)
        os.makedirs(batch_out, exist_ok=True)
        shutil.copy2(src, os.path.join(batch_out, filename))

        copied += 1

    print(f"\n  [OK] {copied} passed .pdbqt files copied to:")
    print(f"       {all_dir}  (flat / all batches combined)")
    print(f"       {passed_dir}/<batch>/  (per-batch folders)")
    if missing:
        print(f"  [!] {missing} source files could not be located.")
    return copied


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Main pipeline
# ─────────────────────────────────────────────────────────────────────────────

def run(input_dir, output_file, passed_dir=None):
    records      = []
    source_paths = {}   # mol_name -> original .pdbqt path

    batches = sorted(glob.glob(os.path.join(input_dir, "batch*")))
    if not batches:
        batches = [input_dir]   # flat folder fallback

    total = skipped = 0

    for batch_path in batches:
        batch_name  = os.path.basename(batch_path)
        pdbqt_files = sorted(glob.glob(os.path.join(batch_path, "*.pdbqt")))
        print(f"  {batch_name}: {len(pdbqt_files)} files found")

        for fp in pdbqt_files:
            name = os.path.splitext(os.path.basename(fp))[0]
            source_paths[name] = fp          # record BEFORE any filtering
            mol  = pdbqt_to_mol(fp)
            total += 1

            if mol is None:
                skipped += 1
                records.append({"Molecule": name, "Batch": batch_name,
                                 "Overall_ADMET_Pass": "PARSE_ERROR"})
                continue
            try:
                records.append(calc_admet(mol, name, batch_name))
            except Exception as e:
                skipped += 1
                records.append({"Molecule": name, "Batch": batch_name,
                                 "Overall_ADMET_Pass": f"ERROR: {e}"})

    df = pd.DataFrame(records)

    priority_cols = [
        "Molecule", "Batch", "Overall_ADMET_Pass",
        "MW (Da)", "LogP", "HBD", "HBA", "TPSA (A2)", "RotBonds",
        "QED", "Lipinski_Pass", "Veber_Pass", "Ghose_Pass", "Egan_Pass",
        "Muegge_Pass", "BBB_Permeable", "GI_Absorption", "Pgp_Substrate",
        "hERG_Risk", "Ames_Alert", "PAINS_Flag", "PAINS_Alerts",
        "Brenk_Flag", "Brenk_Alerts", "SA_Complexity", "LeadLike",
        "Rings", "AromaticRings", "Fsp3", "HeavyAtomCount", "Exact_MW (Da)",
    ]
    cols = [c for c in priority_cols if c in df.columns] + \
           [c for c in df.columns   if c not in priority_cols]
    df = df[cols]

    # ── Excel output ──────────────────────────────────────────────────────
    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        df.to_excel(writer, index=False, sheet_name="ADMET_Results")
        ws = writer.sheets["ADMET_Results"]

        from openpyxl.styles import PatternFill, Font, Alignment
        from openpyxl.utils import get_column_letter

        green       = PatternFill("solid", fgColor="C6EFCE")
        red         = PatternFill("solid", fgColor="FFC7CE")
        amber       = PatternFill("solid", fgColor="FFEB9C")
        header_fill = PatternFill("solid", fgColor="1F4E79")

        for cell in ws[1]:
            cell.fill      = header_fill
            cell.font      = Font(bold=True, color="FFFFFF", size=9)
            cell.alignment = Alignment(horizontal="center", wrap_text=True)

        for col_idx, col in enumerate(df.columns, 1):
            max_len = max(len(str(col)), df[col].astype(str).str.len().max())
            ws.column_dimensions[get_column_letter(col_idx)].width = min(max_len + 2, 30)

        for row in ws.iter_rows(min_row=2):
            for cell in row:
                hdr = ws.cell(row=1, column=cell.column).value
                val = str(cell.value) if cell.value else ""
                if hdr == "Overall_ADMET_Pass":
                    cell.fill = green if val == "PASS" else (amber if "ERROR" in val else red)
                elif hdr in ("Lipinski_Pass","Veber_Pass","Ghose_Pass","Egan_Pass","Muegge_Pass","LeadLike"):
                    cell.fill = green if val == "YES" else red
                elif hdr == "BBB_Permeable":
                    cell.fill = green if val == "YES" else amber
                elif hdr == "GI_Absorption":
                    cell.fill = green if val == "HIGH" else amber
                elif hdr in ("PAINS_Flag","Brenk_Flag","Ames_Alert"):
                    cell.fill = red if val == "YES" else green
                elif hdr == "hERG_Risk":
                    cell.fill = red if val == "HIGH" else green
                elif hdr == "Pgp_Substrate":
                    cell.fill = amber if val == "YES" else green

        # Summary sheet
        n_pass = int((df["Overall_ADMET_Pass"] == "PASS").sum())
        n_fail = int((df["Overall_ADMET_Pass"] == "FAIL").sum())
        summary = {
            "Total molecules processed": [total],
            "Parse errors / skipped":    [skipped],
            "ADMET PASS":                [n_pass],
            "ADMET FAIL":                [n_fail],
            "Lipinski Pass":  [int((df.get("Lipinski_Pass", pd.Series()) == "YES").sum())],
            "Veber Pass":     [int((df.get("Veber_Pass",    pd.Series()) == "YES").sum())],
            "PAINS flagged":  [int((df.get("PAINS_Flag",    pd.Series()) == "YES").sum())],
            "hERG HIGH risk": [int((df.get("hERG_Risk",     pd.Series()) == "HIGH").sum())],
            "BBB Permeable":  [int((df.get("BBB_Permeable", pd.Series()) == "YES").sum())],
        }
        pd.DataFrame(summary).T.rename(columns={0: "Count"}).to_excel(
            writer, sheet_name="Summary")
        ws2 = writer.sheets["Summary"]
        for cell in ws2[1]:
            cell.fill = header_fill
            cell.font = Font(bold=True, color="FFFFFF")

    # ── Copy passed files ──────────────────────────────────────────────────
    if passed_dir is None:
        passed_dir = os.path.join(os.path.dirname(os.path.abspath(output_file)),
                                  "admet_passed")

    copied = copy_passed_files(df, source_paths, passed_dir)

    print(f"\n{'='*62}")
    print(f"  Processed   : {total} molecules  ({skipped} skipped/errored)")
    print(f"  ADMET PASS  : {n_pass}")
    print(f"  ADMET FAIL  : {n_fail}")
    print(f"  Files copied: {copied}")
    print(f"  Excel output: {output_file}")
    print(f"  Passed dir  : {passed_dir}/")
    print(f"{'='*62}\n")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# 5.  CLI
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="ADMET screening for AutoDock Vina .pdbqt ligand files")
    parser.add_argument("--input_dir",  required=True,
        help="Root folder containing batch1/, batch2/, ... (or flat .pdbqt files)")
    parser.add_argument("--output",     default="admet_results.xlsx",
        help="Output Excel file  (default: admet_results.xlsx)")
    parser.add_argument("--passed_dir", default=None,
        help="Destination folder for ADMET-passed .pdbqt files  "
             "(default: <output_dir>/admet_passed/)")
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        sys.exit(f"ERROR: '{args.input_dir}' is not a valid directory.")

    print(f"\nADMET Screening Pipeline")
    print(f"  Input      : {args.input_dir}")
    print(f"  Excel out  : {args.output}")
    print(f"  Passed dir : {args.passed_dir or '<output_dir>/admet_passed/'}\n")

    run(args.input_dir, args.output, args.passed_dir)
