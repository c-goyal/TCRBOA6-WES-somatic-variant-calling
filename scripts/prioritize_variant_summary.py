import pandas as pd

# Input: the table containing all 45 annotated variants
input_file = (
    "/mnt/c/variant-calling-work/annotation/results/"
    "variant_summary.tsv"
)

# Output: one final table with priority tiers and reasons
output_file = (
    "/mnt/c/variant-calling-work/annotation/results/"
    "variant_summary_final.tsv"
)

# Read the input table
df = pd.read_csv(input_file, sep="\t")


# Assign a consequence rank and explanation to each variant
def impact_rank(consequence):
    if consequence in [
        "FRAME_SHIFT_DEL",
        "FRAME_SHIFT_INS",
        "NONSENSE",
        "SPLICE_SITE",
    ]:
        return 4, "HIGH: predicted consequence"

    if consequence == "MISSENSE":
        return 2, "MODERATE: missense"

    if consequence == "SILENT":
        return 1, "LOW: silent"

    return 1, "LOW/UNCERTAIN: other consequence"


# Check whether a database annotation is actually present
def has_annotation(value):
    return (
        pd.notna(value)
        and str(value).strip() not in {"", ".", "NA"}
    )


# Assign a priority tier to every variant
tiers = []
reasons = []

for _, row in df.iterrows():

    impact, reason = impact_rank(
        row["Variant_Classification"]
    )

    evidence_count = sum([
        has_annotation(row["CGC_Name"]),
        has_annotation(row["ClinVar_Significance"]),
        has_annotation(row["COSMIC_Annotation"]),
    ])

    if impact == 4 and evidence_count > 0:
        tier = 1
        tier_reason = reason + " + database annotation"

    elif impact == 4 or (
        impact == 2 and evidence_count > 0
    ):
        tier = 2
        tier_reason = reason

        if evidence_count > 0:
            tier_reason += " + database annotation"

    elif impact == 2:
        tier = 3
        tier_reason = reason

    else:
        tier = 4
        tier_reason = reason

    tiers.append(tier)
    reasons.append(tier_reason)


# Add the priority information to the table
df["PRIORITY_TIER"] = tiers
df["PRIORITY_REASON"] = reasons

# Ensure numerical columns are treated as numbers
df["TUMOR_AF"] = pd.to_numeric(
    df["TUMOR_AF"], errors="coerce"
)

df["TLOD"] = pd.to_numeric(
    df["TLOD"], errors="coerce"
)

# Sort by priority, then tumour allele frequency,
# then TLOD, with higher values first within each tier
df = df.sort_values(
    by=["PRIORITY_TIER", "TUMOR_AF", "TLOD"],
    ascending=[True, False, False],
    na_position="last",
)

# Save the final table
df.to_csv(
    output_file,
    sep="\t",
    index=False,
)

# Print a brief report
print(f"Final table saved: {output_file}")
print(f"Total variants: {len(df)}")

print("\nVariants per priority tier:")
print(
    df["PRIORITY_TIER"]
    .value_counts()
    .sort_index()
    .to_string()
)

print("\nFirst 10 rows:")
print(
    df[[
        "Gene",
        "Variant_Classification",
        "TUMOR_AF",
        "TLOD",
        "PRIORITY_TIER",
        "PRIORITY_REASON",
    ]]
    .head(10)
    .to_string(index=False)
)
