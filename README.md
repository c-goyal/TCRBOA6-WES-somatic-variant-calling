# TCRBOA6 WES Somatic Variant Calling

## Overview

This project analyzes paired tumor-normal whole-exome sequencing (WES) data from the TCRBOA6 sample to identify, filter, annotate, and prioritize candidate somatic variants.

The workflow was implemented using commonly used bioinformatics tools, including BWA, samtools, GATK, Funcotator, Python, and Linux command-line utilities.

### Samples

- Tumor: `SRR2184243`
- Matched normal: `SRR2184242`
- Sequencing: Paired-end WES
- Capture platform: NimbleGen VCRome 2.1
- Reference genome: Broad b37 / `Homo_sapiens_assembly19.fasta`

---

## Workflow

```text
Raw FASTQ
   ↓
FastQC + MultiQC
   ↓
BWA-MEM alignment
   ↓
Coordinate sorting
   ↓
MarkDuplicates
   ↓
Base Quality Score Recalibration (BQSR)
   ↓
Mutect2 somatic variant calling
   ↓
Orientation bias + contamination estimation
   ↓
FilterMutectCalls
   ↓
Funcotator annotation
   ↓
Variant prioritization
   ↓
Candidate interpretation
```

---

## Quality Control

FastQC and MultiQC were used to evaluate the tumor and normal FASTQ files.

All four FASTQ files passed the major FastQC quality modules. Reads were 101 bp long, adapter contamination was negligible, and sequence duplication was approximately 20%.

A consistent warning was observed for per-sequence GC content across the four libraries. Because the warning was systematic and the other major QC metrics passed, the data were considered suitable for downstream analysis.

---

## Somatic Variant Calling

Tumor and matched-normal reads were aligned to the Broad b37 reference genome using BWA-MEM.

GATK Mutect2 was then used for tumor-normal somatic variant calling using:

- gnomAD population allele-frequency resource
- Mutect2 panel of normals
- VCRome 2.1 target intervals

Mutect2 initially identified:

**1,251 candidate variant records**

After filtering with `FilterMutectCalls`:

**45 variants passed the somatic filtering criteria.**

Estimated tumor contamination was:

**1.47%**

---

## Functional Annotation and Prioritization

The 45 PASS variants were annotated using GATK Funcotator with functional, cancer-related, and population resources including GENCODE, Cancer Gene Census (CGC), ClinVar, COSMIC, HGNC, dbSNP, and additional annotation datasets.

A consequence- and annotation-based heuristic was then used to prioritize variants with potentially important functional effects.

This prioritization is a research-oriented ranking strategy and is not a clinical pathogenicity classification.

### Consequence Summary

| Consequence | Count | Percentage |
|---|---:|---:|
| MISSENSE | 27 | 60.0% |
| INTRON | 4 | 8.9% |
| FRAME_SHIFT_DEL | 3 | 6.7% |
| FIVE_PRIME_UTR | 3 | 6.7% |
| THREE_PRIME_UTR | 2 | 4.4% |
| SILENT | 2 | 4.4% |
| SPLICE_SITE | 1 | 2.2% |
| FRAME_SHIFT_INS | 1 | 2.2% |
| NONSENSE | 1 | 2.2% |
| RNA | 1 | 2.2% |

---

## High-Priority Candidates

Seven variants were selected for further interpretation.

| Gene | Consequence | Coding change | Protein change | Tumor AF | Normal AF |
|---|---|---|---|---:|---:|
| MEN1 | Frameshift insertion | c.402_403insC | p.K135fs | 0.588 | 0.00479 |
| VHL | Frameshift deletion | c.375_379delCGATG | p.H125fs | 0.541 | 0.016 |
| CEACAM8 | Frameshift deletion | c.48_61delGGGGCTCCTGCTCA | p.Q16fs | 0.397 | 0.019 |
| TMC2 | Nonsense | c.493A>T | p.K165* | 0.322 | 0.00439 |
| LOXL3 | Frameshift deletion | c.937_1015del... | p.A313fs | 0.176 | 0.011 |
| CLTC | Missense | c.3160G>C | p.D1054H | 0.075 | 0.012 |
| TMCO1 | Splice-site | c.e9+1G>A | — | 0.055 | 0.011 |

---

## Candidate Interpretation

### MEN1

A frameshift insertion in MEN1 (`p.K135fs`) showed strong tumor enrichment and strong Mutect2 evidence. Funcotator provided Cancer Gene Census annotation and COSMIC-overlapping mutation information.

An exact ClinVar match was not identified for this allele.

**Interpretation:** High-priority somatic candidate with a potentially disruptive frameshift consequence. The available evidence does not establish clinical pathogenicity or actionability for the exact allele.

### VHL

A frameshift deletion in VHL (`p.H125fs`) showed strong tumor enrichment and strong Mutect2 evidence. Cancer Gene Census and COSMIC-overlapping annotations were also present.

An exact ClinVar match was not identified for this allele.

**Interpretation:** High-priority somatic candidate with a potentially disruptive frameshift consequence. Exact clinical significance remains unestablished.

### CEACAM8

A frameshift deletion (`p.Q16fs`) showed substantial tumor enrichment.

No exact ClinVar interpretation was identified.

**Interpretation:** Candidate of interest based primarily on the predicted loss-of-function consequence and tumor enrichment.

### TMC2

A nonsense variant (`p.K165*`) showed substantial tumor enrichment.

No exact ClinVar interpretation was identified.

**Interpretation:** Candidate of interest based on the predicted premature stop consequence and tumor enrichment.

### LOXL3

A frameshift deletion (`p.A313fs`) showed tumor enrichment.

No exact ClinVar interpretation was identified.

**Interpretation:** Candidate of interest requiring further investigation.

### CLTC

A missense variant (`p.D1054H`) showed lower tumor allele fraction than the other leading candidates. CGC gene-level annotation was present.

No exact ClinVar interpretation was identified.

**Interpretation:** Lower-priority candidate relative to the predicted loss-of-function variants.

### TMCO1

A canonical splice-site variant was identified with tumor enrichment. The exact genomic variant was represented in ClinVar, but the record contained:

- no germline classification
- no submitted somatic clinical impact
- no submitted oncogenicity interpretation

**Interpretation:** Candidate of interest because of its canonical splice-site consequence, but its clinical and oncogenic significance remains unestablished.

---

## Overall Results

The workflow produced:

```text
1,251 Mutect2 candidate records
            ↓
45 PASS variants
            ↓
7 prioritized candidates
```

MEN1 and VHL were the strongest computationally prioritized candidates based on their predicted frameshift consequences, tumor enrichment, Mutect2 evidence, and cancer-related database annotations.

However, these findings represent **computationally prioritized somatic candidates**, not confirmed pathogenic or clinically actionable variants.

---

## Limitations

- This analysis is based on a single tumor-normal WES pair.
- WES does not capture the complete genome.
- Variant prioritization was performed using a heuristic rather than a validated clinical classification framework.
- Absence from ClinVar does not imply that a variant is benign.
- Cancer Gene Census and COSMIC annotations do not by themselves establish pathogenicity of an exact allele.
- Orthogonal validation and additional biological or clinical evidence would be required for clinical interpretation.
- IGV read-level inspection was not included in the initial computational analysis and is planned as an additional validation step.

---

## Repository Structure

```text
TCRBOA6-WES-somatic-variant-calling/
├── README.md
├── scripts/
│   ├── annotation_summary.py
│   ├── prioritize_variants.py
│   └── make_final_tables.py
├── results/
│   ├── PASS_variants_annotated_summary.tsv
│   ├── PASS_variants_prioritized.tsv
│   ├── high_priority_candidates.tsv
│   ├── variant_consequence_summary.tsv
│   ├── final_variant_interpretation.tsv
│   └── contamination.table
└── docs/
```

Large sequencing, BAM, reference, database, and temporary files are intentionally excluded from the repository.

---

## Tools

- FastQC
- MultiQC
- BWA-MEM
- samtools
- GATK 4.7.0.0
- Funcotator
- Python
- Bash/Linux command-line tools

## Project Status

**Computational somatic variant calling, annotation, and prioritization completed.**

IGV-based read-level validation is planned as an additional analysis step.
