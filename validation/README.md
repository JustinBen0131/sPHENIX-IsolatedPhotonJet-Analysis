# Migration and content comparison

This directory is separate from the permanent event-production path. It supports the THE358 → THE352 closure, not a second production framework.

`field_mapping.csv` retains all **1,398** historical branch rows and their original `proposed_disposition`. Every actual disposition is **REVIEW**; none has been scientifically accepted. Names alone do not establish equivalence. Reviewers must bind source evidence, canonical location and validity; DERIVED needs exact source fields/formula, POST-TREE needs an owner, and ARCHIVE needs a reason.

`convert_legacy_to_canonical.py` currently audits that ledger. `--audit-only` reports missing review requirements and returns nonzero while unresolved. Actual conversion refuses before creating ROOT output. Layout-specific identity reconstruction, event/object reshaping, joins, validity translation and a complete schema writer are still required. The first pass's generic first-alias column copy was removed because it could silently lose information.

`compare_canonical.py` joins rows by all required identity columns. It compares every table/field/type, exact integers/enums, jagged shapes and values, explicit floating tolerances, metadata and completion accounting. Missing tables/fields/keys, duplicate keys, unsupported tables and differences cause a nonzero result. NaN, positive infinity and negative infinity remain distinct. File SHA256 equality is reported separately from normalized content equality. Build-specific provenance is retained in the report with a small explicit exemption list; additional allowed differences require review, not wildcard ignores.

These tools have had source/static checks only. No ROOT comparison or migration was executed.

## Closure order

1. Freeze the final repaired production reference and accepted model/science bindings.
2. Review every ledger disposition; no required information disappears by inference.
3. Implement one reviewed converter with exact layout/identity/validity adapters.
4. Prove required information survives conversion.
5. Compile/run the direct producer on representative pp DATA, pp photon/inclusive SIM, both pp DI lanes, AuAu DATA and both AuAu embedded SIM lanes.
6. Compare converted reference versus direct output: populations, identities, types, values, validity, weights, associations, provenance and downstream physics.
7. Only after acceptance change the active source of truth; archive historical repair routes separately.

Byte-identical ROOT serialization is not the closure criterion. A content-comparator pass also does not establish scientific acceptance without the reviewed mapping and downstream checks.
