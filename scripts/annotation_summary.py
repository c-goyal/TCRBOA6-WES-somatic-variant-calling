import gzip
import csv
import re

VCF = "/mnt/c/variant-calling-work/annotation/results/tumor.annotated.vcf.gz"
OUTPUT = "/mnt/c/variant-calling-work/annotation/results/PASS_variants_annotated_summary.tsv"

wanted_fields = [
    "Gencode_43_hugoSymbol",
    "Gencode_43_variantClassification",
    "Gencode_43_cDnaChange",
    "Gencode_43_proteinChange",
    "CGC_Name",
    "ClinVar_VCF_CLNSIG",
    "Cosmic_overlapping_mutations",
    "dbSNP_RS",
]

output_header = [
    "CHROM", "POS", "REF", "ALT",
    "NORMAL_DP", "NORMAL_AD", "NORMAL_AF",
    "TUMOR_DP", "TUMOR_AD", "TUMOR_AF",
    "TLOD", "NLOD",
    "GENE", "VARIANT_CLASSIFICATION",
    "CODING_CHANGE", "PROTEIN_CHANGE",
    "CGC_NAME", "CLINVAR_SIGNIFICANCE",
    "COSMIC_OVERLAPPING", "dbSNP_RS"
]

func_indexes = {}
normal_idx = None
tumor_idx = None
count = 0

with gzip.open(VCF, "rt") as vcf, open(OUTPUT, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(output_header)

    for line in vcf:
        line = line.rstrip("\n")

        # Skip blank lines
        if not line.strip():
            continue

        # Read the Funcotator field definitions from the header
        if line.startswith("##INFO=<ID=FUNCOTATION"):
            match = re.search(r"Funcotation fields are: (.*)\">", line)
            if match:
                func_fields = match.group(1).split("|")
                func_indexes = {
                    name: i for i, name in enumerate(func_fields)
                }
            continue

        # Find sample column positions
        if line.startswith("#CHROM"):
            columns = line.split("\t")
            normal_idx = columns.index("normal")
            tumor_idx = columns.index("tumor")
            continue

        # Skip all other header lines
        if line.startswith("#"):
            continue

        fields = line.split("\t")

        # Skip malformed lines
        if len(fields) < 11:
            continue

        # Keep PASS variants only
        if fields[6] != "PASS":
            continue

        # -------------------------
        # Parse INFO
        # -------------------------
        info = {}
        for item in fields[7].split(";"):
            if "=" in item:
                key, value = item.split("=", 1)
                info[key] = value

        # -------------------------
        # Parse FUNCOTATION
        # -------------------------
        raw_func = info.get("FUNCOTATION", "")
        raw_func = raw_func.strip("[]")

        # For this project we have one ALT per record.
        # Take the first Funcotation entry if multiple are present.
        first_annotation = raw_func.split("],[")[0].strip("[]")
        func_values = first_annotation.split("|")

        def get_func(name):
            idx = func_indexes.get(name)
            if idx is None:
                return ""
            if idx >= len(func_values):
                return ""
            value = func_values[idx]
            if value in ("", "NONE", "NA"):
                return ""
            return value

        # -------------------------
        # Parse normal/tumor sample fields
        # -------------------------
        format_keys = fields[8].split(":")
        normal_values = fields[normal_idx].split(":")
        tumor_values = fields[tumor_idx].split(":")

        normal = dict(zip(format_keys, normal_values))
        tumor = dict(zip(format_keys, tumor_values))

        row = [
            fields[0],
            fields[1],
            fields[3],
            fields[4],

            normal.get("DP", ""),
            normal.get("AD", ""),
            normal.get("AF", ""),

            tumor.get("DP", ""),
            tumor.get("AD", ""),
            tumor.get("AF", ""),

            info.get("TLOD", ""),
            info.get("NLOD", ""),

            get_func("Gencode_43_hugoSymbol"),
            get_func("Gencode_43_variantClassification"),
            get_func("Gencode_43_cDnaChange"),
            get_func("Gencode_43_proteinChange"),
            get_func("CGC_Name"),
            get_func("ClinVar_VCF_CLNSIG"),
            get_func("Cosmic_overlapping_mutations"),
            get_func("dbSNP_RS"),
        ]

        writer.writerow(row)
        count += 1

print(f"Wrote {count} PASS variants to:")
print(OUTPUT)
