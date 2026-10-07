import csv
from collections import Counter
from urllib.parse import unquote

results_dir = "/mnt/c/variant-calling-work/annotation/results"

annotated_file = results_dir + "/PASS_variants_annotated_summary.tsv"
priority_file = results_dir + "/high_priority_candidates.tsv"

# --------------------------------------------------
# TABLE 1: variant_consequence_summary.tsv
# --------------------------------------------------

with open(annotated_file, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

consequences = Counter(row["VARIANT_CLASSIFICATION"] for row in rows)
total = len(rows)

summary_file = results_dir + "/variant_consequence_summary.tsv"

with open(summary_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(["VARIANT_CLASSIFICATION", "COUNT", "PERCENT"])

    for consequence, count in sorted(
        consequences.items(),
        key=lambda x: x[1],
        reverse=True
    ):
        percent = count / total * 100
        writer.writerow([consequence, count, f"{percent:.1f}"])

# --------------------------------------------------
# TABLE 2: final_variant_interpretation.tsv
# --------------------------------------------------

clinvar = {
    "MEN1": (
        "No exact match found",
        "No exact ClinVar interpretation identified"
    ),
    "VHL": (
        "No exact match found",
        "No exact ClinVar interpretation identified"
    ),
    "CEACAM8": (
        "No exact match found",
        "No exact ClinVar interpretation identified"
    ),
    "TMC2": (
        "No exact match found",
        "No exact ClinVar interpretation identified"
    ),
    "LOXL3": (
        "No exact match found",
        "No exact ClinVar interpretation identified"
    ),
    "CLTC": (
        "No exact match found",
        "No exact ClinVar interpretation identified"
    ),
    "TMCO1": (
        "Exact record found",
        "No germline classification; no somatic clinical impact; "
        "no oncogenicity submitted"
    ),
}

assessment = {
    "MEN1":
        "High-priority somatic candidate; frameshift with strong "
        "tumor enrichment; cancer/database support; exact clinical "
        "interpretation not established.",

    "VHL":
        "High-priority somatic candidate; frameshift with strong "
        "tumor enrichment; cancer/database support; exact clinical "
        "interpretation not established.",

    "CEACAM8":
        "Candidate of interest; frameshift with strong tumor "
        "enrichment; no exact ClinVar interpretation identified.",

    "TMC2":
        "Candidate of interest; nonsense variant with strong tumor "
        "enrichment; no exact ClinVar interpretation identified.",

    "LOXL3":
        "Candidate of interest; frameshift with tumor enrichment; "
        "no exact ClinVar interpretation identified.",

    "CLTC":
        "Candidate of interest; missense variant with lower tumor "
        "AF and CGC gene-level annotation; exact clinical "
        "interpretation not established.",

    "TMCO1":
        "Candidate of interest; canonical splice-site variant with "
        "tumor enrichment; exact ClinVar record exists but has no "
        "submitted clinical interpretation.",
}

priority_genes = [
    "MEN1",
    "VHL",
    "CEACAM8",
    "TMC2",
    "LOXL3",
    "CLTC",
    "TMCO1",
]

with open(priority_file, "r") as f:
    reader = csv.DictReader(f, delimiter="\t")
    priority_rows = list(reader)

priority_rows.sort(
    key=lambda row: priority_genes.index(row["GENE"])
)

final_file = results_dir + "/final_variant_interpretation.tsv"

with open(final_file, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")

    writer.writerow([
        "GENE",
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "VARIANT_CLASSIFICATION",
        "CODING_CHANGE",
        "PROTEIN_CHANGE",
        "TUMOR_AF",
        "NORMAL_AF",
        "TLOD",
        "NLOD",
        "CGC_EVIDENCE",
        "COSMIC_OVERLAP",
        "CLINVAR_EXACT_MATCH",
        "CLINVAR_INTERPRETATION",
        "FINAL_ASSESSMENT"
    ])

    for row in priority_rows:

        gene = row["GENE"]

        cgc = "Yes" if row["CGC_NAME"].strip() else "No"
        cosmic = "Yes" if row["COSMIC_OVERLAPPING"].strip() else "No"

        clinvar_match, clinvar_interpretation = clinvar[gene]

        writer.writerow([
            gene,
            row["CHROM"],
            row["POS"],
            row["REF"],
            row["ALT"],
            row["VARIANT_CLASSIFICATION"],
            row["CODING_CHANGE"],
            row["PROTEIN_CHANGE"],
            row["TUMOR_AF"],
            row["NORMAL_AF"],
            row["TLOD"],
            row["NLOD"],
            cgc,
            cosmic,
            clinvar_match,
            clinvar_interpretation,
            assessment[gene]
        ])

print("Created:", summary_file)
print("Created:", final_file)
print("PASS variants:", total)
print("Prioritized candidates:", len(priority_rows))
