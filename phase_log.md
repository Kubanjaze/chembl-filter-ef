# Phase 95 — Integration: ChEMBL → RDKit Filter → EF@K
## Phase Log

**Status:** ✅ Complete
**Started:** 2026-03-27
**Completed:** 2026-03-27
**Repo:** https://github.com/Kubanjaze/chembl-filter-ef

---

## Log

### 2026-03-27 12:15 — Plan written, initial push
- Implementation plan v1.0 written
- Combines Phases 80 + 09 + 35 patterns

### 2026-03-27 12:23 — Build complete
- 200 records → 118 unique → 8 pass drug-likeness → 4 active
- EF@10%=2.0 (double random enrichment)
- Key insight: RO5 filters too aggressive for KRAS covalent inhibitor space
- Cost: $0.00
