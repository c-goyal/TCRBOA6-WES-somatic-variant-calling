import csv

INPUT = "/mnt/c/variant-calling-work/annotation/results/PASS_variants_annotated_summary.tsv"
OUTPUT = "/mnt/c/variant-calling-work/annotation/results/PASS_variants_prioritized.tsv"

def impact_rank(consequence):
    c = consequence.upper()

    if "FRAME_SHIFT" in c:
        return 4, "HIGH: frameshift"
    if "NONSENSE" in c:
        return 4, "HIGH: premature stop"
    if "SPLICE_SITE" in c:
        return 4, "HIGH: splice-site"
    if "MISSENSE" in c:
        return 2, "MODERATE: missense"
    if "SILENT" in c:
        return 1, "LOW: silent"
    if "UTR" in c or "INTRON" in c or "RNA" in c:
        return 1, "LOW/UNCERTAIN: non-coding"
    return 1, "LOW/UNCERTAIN"

with open(INPUT, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

for r in rows:
    impact, reason = impact_rank(r["VARIANT_CLASSIFICATION"])

    has_cgc = bool(r["CGC_NAME"].strip())
    has_clinvar = bool(r["CLINVAR_SIGNIFICANCE"].strip())
    has_cosmic = bool(r["COSMIC_OVERLAPPING"].strip())

    evidence_count = sum([has_cgc, has_clinvar, has_cosmic])

    if impact == 4 and evidence_count > 0:
        tier = 1
        tier_reason = reason + " + cancer/database evidence"
    elif impact == 4 or (impact == 2 and evidence_count > 0):
        tier = 2
        tier_reason = reason + (" + database evidence" if evidence_count > 0 else "")
    elif impact == 2:
        tier = 3
        tier_reason = reason
    else:
        tier = 4
        tier_reason = reason

    r["PRIORITY_TIER"] = str(tier)
    r["PRIORITY_REASON"] = tier_reason

rows.sort(
    key=lambda r: (
        int(r["PRIORITY_TIER"]),
        -float(r["TUMOR_AF"] or 0),
        -float(r["TLOD"] or 0)
    )
)

fieldnames = list(rows[0].keys())

with open(OUTPUT, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} prioritized variants to:")
print(OUTPUT)
