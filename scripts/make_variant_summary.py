
import pandas as pd

# 1. Read the extracted VCF data
df = pd.read_csv(
    "/mnt/c/variant-calling-work/annotation/results/variant_summary_raw.tsv",
    sep="\t"
)

# 2. Remove square brackets from Funcotator annotations
annotation = (
    df["FUNCOTATION"]
    .str.removeprefix("[")
    .str.removesuffix("]")
)

# 3. Split Funcotator into individual fields
ann = annotation.str.split("|", regex=False, expand=True)

# 4. Select the useful annotation fields
selected = pd.DataFrame({
    "Gene": ann[0],
    "Variant_Classification": ann[5],
    "Transcript": ann[12],
    "cDNA_Change": ann[16],
    "Protein_Change": ann[18],
    "CGC_Name": ann[23],
    "ClinVar_Significance": ann[54],
    "COSMIC_Annotation": ann[66],
})

# 5. Combine annotations with the existing variant metrics
summary = pd.concat(
    [df.drop(columns=["FUNCOTATION"]), selected],
    axis=1,
)

# 6. Save the final table
summary.to_csv(
    "/mnt/c/variant-calling-work/annotation/results/variant_summary.tsv",
    sep="\t",
    index=False
)

print(f"Saved {len(summary)} variants.")
print(summary.head(3).to_string(index=False))
