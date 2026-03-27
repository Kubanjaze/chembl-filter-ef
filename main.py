import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os, warnings
warnings.filterwarnings("ignore")
import pandas as pd
import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
RDLogger.DisableLog("rdApp.*")

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


# ── Phase 80 pattern: ChEMBL fetch ─────────────────────────────────────
def fetch_bioactivity(target_id: str, limit: int = 200) -> list[dict]:
    url = f"{CHEMBL_BASE}/activity.json"
    params = {"target_chembl_id": target_id, "limit": min(limit, 1000), "offset": 0}
    records = []
    while len(records) < limit:
        resp = requests.get(url, params=params, timeout=60)
        resp.raise_for_status()
        data = resp.json()
        for act in data.get("activities", []):
            smiles = act.get("canonical_smiles")
            pchembl = act.get("pchembl_value")
            stype = act.get("standard_type")
            if smiles and pchembl and stype in ("IC50", "Ki", "Kd"):
                records.append({
                    "molecule_chembl_id": act.get("molecule_chembl_id"),
                    "canonical_smiles": smiles,
                    "standard_type": stype,
                    "pchembl_value": float(pchembl),
                })
            if len(records) >= limit:
                break
        if not data.get("page_meta", {}).get("next"):
            break
        params["offset"] += params["limit"]
    return records[:limit]


# ── Phase 09 pattern: RDKit property filter ─────────────────────────────
def compute_and_filter(df: pd.DataFrame, mw_max=500, logp_max=5, tpsa_max=140) -> pd.DataFrame:
    props = []
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None:
            props.append({"valid": False})
            continue
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = rdMolDescriptors.CalcTPSA(mol)
        passes = mw <= mw_max and logp <= logp_max and tpsa <= tpsa_max
        props.append({"mw": round(mw, 2), "logp": round(logp, 2), "tpsa": round(tpsa, 2),
                       "valid": True, "passes_filter": passes})
    props_df = pd.DataFrame(props)
    return pd.concat([df.reset_index(drop=True), props_df], axis=1)


# ── Phase 35 pattern: EF@K ──────────────────────────────────────────────
def enrichment_factor(y_true, y_score, k_frac):
    n = len(y_true)
    k = max(1, int(n * k_frac))
    total_hits = sum(y_true)
    if total_hits == 0 or n == 0:
        return 0.0
    ranked = sorted(zip(y_score, y_true), key=lambda x: -x[0])
    hits_topk = sum(label for _, label in ranked[:k])
    ef = (hits_topk / k) / (total_hits / n)
    return round(ef, 3)


def main():
    parser = argparse.ArgumentParser(
        description="Phase 95 — ChEMBL → RDKit Filter → EF@K Integration",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--target", default="CHEMBL2189121", help="ChEMBL target ID")
    parser.add_argument("--limit", type=int, default=200)
    parser.add_argument("--activity-threshold", type=float, default=6.0, help="pChEMBL cutoff for 'active'")
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nPhase 95 — ChEMBL → RDKit Filter → EF@K")
    print(f"Target: {args.target} | Limit: {args.limit}\n")

    # Step 1: Fetch
    print("Step 1: Fetching from ChEMBL...")
    records = fetch_bioactivity(args.target, args.limit)
    df = pd.DataFrame(records)
    print(f"  Records fetched: {len(df)}")

    if df.empty:
        print("No records found.")
        return

    # Deduplicate by molecule
    df = df.sort_values("pchembl_value", ascending=False).drop_duplicates("molecule_chembl_id").reset_index(drop=True)
    print(f"  Unique molecules: {len(df)}")

    # Step 2: Filter
    print("\nStep 2: RDKit property filter...")
    df = compute_and_filter(df)
    valid = df[df["valid"] == True]
    passed = df[df.get("passes_filter", False) == True]
    print(f"  Valid SMILES: {len(valid)}/{len(df)}")
    print(f"  Pass drug-likeness: {len(passed)}/{len(valid)}")

    # Step 3: EF@K
    print(f"\nStep 3: Enrichment Factor (active = pChEMBL >= {args.activity_threshold})...")
    eval_df = passed.copy() if len(passed) > 0 else valid.copy()
    eval_df["active"] = (eval_df["pchembl_value"] >= args.activity_threshold).astype(int)
    n_active = eval_df["active"].sum()
    print(f"  Active compounds: {n_active}/{len(eval_df)} ({n_active/max(len(eval_df),1):.0%})")

    for k_frac in [0.10, 0.20, 0.50]:
        ef = enrichment_factor(eval_df["active"].values, eval_df["pchembl_value"].values, k_frac)
        print(f"  EF@{k_frac:.0%}: {ef}")

    # Save
    eval_df.to_csv(os.path.join(args.output_dir, "filtered_compounds.csv"), index=False)

    report = (
        f"Phase 95 — Integration Pipeline\n{'='*45}\n"
        f"Target: {args.target}\n"
        f"Fetched: {len(records)} records → {len(df)} unique molecules\n"
        f"Filter pass: {len(passed)}/{len(valid)}\n"
        f"Active (pChEMBL>={args.activity_threshold}): {n_active}/{len(eval_df)}\n"
    )
    for k_frac in [0.10, 0.20, 0.50]:
        ef = enrichment_factor(eval_df["active"].values, eval_df["pchembl_value"].values, k_frac)
        report += f"EF@{k_frac:.0%}: {ef}\n"
    report += f"Cost: $0.00\n"
    print(f"\n{report}")

    with open(os.path.join(args.output_dir, "enrichment_report.txt"), "w") as f:
        f.write(report)
    print("Done.")


if __name__ == "__main__":
    main()
