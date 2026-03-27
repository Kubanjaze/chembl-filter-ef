# Phase 95 — ChEMBL Fetch → RDKit Filter → EF@K Evaluation

**Version:** 1.0 | **Tier:** Standard | **Date:** 2026-03-27

## Goal
Integrate three earlier phases into a single pipeline: (1) fetch KRAS bioactivity from ChEMBL (Phase 80), (2) apply RDKit property filters (Phase 09 pattern), (3) evaluate enrichment factor (Phase 35 pattern). Demonstrates end-to-end data→filter→evaluate workflow.

CLI: `python main.py --target CHEMBL2189121 --limit 200`

Outputs: filtered_compounds.csv, enrichment_report.txt

## Logic
1. Fetch bioactivity data from ChEMBL for KRAS (reuse Phase 80 fetch logic)
2. Filter: keep IC50 records with pChEMBL values, compute RDKit properties (MW, LogP, TPSA)
3. Apply drug-likeness filter: MW<500, LogP<5, TPSA<140 (RO5-like)
4. Define "active" as pChEMBL >= 6.0 (IC50 ≤ 1 μM)
5. Compute EF@K at top-10%, top-20%, top-50% of filtered compounds (ranked by pChEMBL)
6. Report: filter pass rate, active count, EF values

## Key Concepts
- **Integration pattern**: chaining Phase 80 (fetch) → Phase 09 (filter) → Phase 35 (EF@K)
- **EF@K**: (hits in top K / K) / (total hits / total N) — measures enrichment vs random
- **Drug-likeness filter**: MW, LogP, TPSA thresholds applied via RDKit descriptors
- **pChEMBL as activity score**: higher = more potent, threshold ≥6.0 for "active"

## Verification Checklist
- [ ] ChEMBL data fetched successfully
- [ ] RDKit filters applied (count pass/fail)
- [ ] EF@K computed at multiple cutoffs
- [ ] EF > 1.0 indicates enrichment above random
- [ ] --help works

## Risks
- ChEMBL may return few IC50 records for KRAS — use --limit 200 to get more
- Most KRAS compounds may pass RO5 (small molecules) — filter may not be very selective
- EF may equal 1.0 if actives are uniformly distributed in score space
