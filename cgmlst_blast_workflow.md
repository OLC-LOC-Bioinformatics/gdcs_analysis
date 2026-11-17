# cgMLST BLAST workflow: database search to result interpretation

This guide documents the workflow from the curated allele database through BLAST searches, group-level summaries, and the final per-genome low-identity/absent-locus report.

## 1. Inputs and paths

| Input/output | Path |
|---|---|
| Genome assemblies | `/mnt/nas2/processed_sequence_data/ncbi/salmonella_benchmark/BestAssemblies/*.fasta and included reference GCF_002761035.1` |
| Curated allele FASTA | `/mnt/nas/users/madhu/salmonella_ncbi_bds_blast/cgmlst-centroid.fasta` |
| BLAST database prefix | `/mnt/nas/users/madhu/salmonella_ncbi_bds_blast/cgmlst_loci_db` |
| Per-genome BLAST results | `/mnt/nas/users/madhu/salmonella_ncbi_bds_blast/cgmlst_results/` |
| Summary outputs | `/mnt/nas/users/madhu/salmonella_ncbi_bds_blast/cgmlst_summary/` |

The uploaded curated FASTA contained 14,766 unique sequence records across 330 cgMLST groups. Each FASTA identifier has this form:

```text
group_id|curated_allele_id
```

For example, `NZ_CM001471.1_3941|1891423980` belongs to group `NZ_CM001471.1_3941`; `1891423980` is that reference sequence's curated allele identifier. All records sharing the prefix before `|` are grouped together. The script uses the complete FASTA ID to match BLAST's `sseqid`.

## 2. Software

- BLAST+ (`makeblastdb`, `blastn`)
- Python 3 (the summary script uses only the standard library)

## 3. Build the nucleotide database

The database already exists in the project directory. Rebuild it only if needed, using the same FASTA that will be supplied to the summary script:

```bash
cd /mnt/nas/users/madhu/salmonella_ncbi_bds_blast

makeblastdb \
  -in cgmlst-centroid.fasta \
  -dbtype nucl \
  -out cgmlst_loci_db
```

The full FASTA sequence IDs must be preserved; they are used to map BLAST subjects to groups and curated allele IDs.

## 4. Run BLAST for each genome

The following is the provided `blast.sh`. Each genome contig is a BLAST query, and each curated allele sequence in the database is a subject. The command creates one headerless, tab-separated result file per genome.

```bash
#!/bin/bash

GENOME_DIR="/mnt/nas2/processed_sequence_data/ncbi/salmonella_benchmark/BestAssemblies"
DB="/mnt/nas/users/madhu/salmonella_ncbi_bds_blast/cgmlst_loci_db"
OUTDIR="/mnt/nas/users/madhu/salmonella_ncbi_bds_blast/cgmlst_results"
THREADS=7

mkdir -p "$OUTDIR"

shopt -s nullglob
genomes=("$GENOME_DIR"/*.fasta)

if [[ ${#genomes[@]} -eq 0 ]]; then
    echo "ERROR: No FASTA files found in $GENOME_DIR"
    exit 1
fi

for genome in "${genomes[@]}"; do
    filename=$(basename "$genome" .fasta)
    output="$OUTDIR/${filename}.tsv"

    echo "[$(date '+%F %T')] Processing $filename"

    # Skip existing results
    if [[ -s "$output" ]]; then
        echo "Already completed: $filename"
        continue
    fi

    blastn \
        -query "$genome" \
        -db "$DB" \
        -task blastn \
        -evalue 1e-10 \
        -num_threads "$THREADS" \
        -max_target_seqs 20000 \
        -max_hsps 1 \
        -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore qcovhsp" \
        -out "$output"

    if [[ $? -eq 0 ]]; then
        echo "Completed: $filename"
    else
        echo "ERROR: BLAST failed for $filename"
        rm -f "$output"
    fi
done
```

Run it if the per-genome BLAST files have not already been generated:

```bash
bash blast.sh
```

The script skips any nonempty existing result file. If a previous result may be incomplete, check it and the run log before relying on it; the script will not rerun it automatically.

### BLAST output columns

| Column | Meaning |
|---|---|
| `qseqid` | Genome contig identifier |
| `sseqid` | Full curated allele sequence ID from the database |
| `pident` | Percent identity over this local alignment |
| `length` | Alignment length |
| `qlen` | Genome contig length |
| `qstart`, `qend` | Alignment coordinates on the genome contig |
| `sstart`, `send` | Alignment coordinates on the curated allele subject |
| `evalue` | BLAST E-value |
| `bitscore` | BLAST bit score |
| `qcovhsp` | Percentage of the genome contig covered by this HSP; it is **not** allele coverage |

The summary script calculates subject-coordinate span coverage as:

```text
100 × (abs(send - sstart) + 1) / curated_allele_sequence_length
```

This is the span covered on one HSP, not coverage assembled from multiple alignments.

## 5. Summarize BLAST results

Copy the current `summarize_cgmlst_by_group.py` into the project directory and run:

```bash
cd /mnt/nas/users/madhu/salmonella_ncbi_bds_blast

python3 summarize_cgmlst_by_group.py \
  --alleles cgmlst-centroid.fasta \
  --results cgmlst_results \
  --out cgmlst_summary
```

The script reads every `*.tsv` in `cgmlst_results`. Each filename without `.tsv` becomes the genome name. Keep unrelated TSV files out of this directory. The number of result TSVs should match the genomes you intend to summarize; missing files are not included in the denominator.

The script groups reference sequences using the part of the FASTA ID before `|`. For each genome and group it processes the reported BLAST hits. The `--genome-min-identity` setting defaults to `80`; in `genome_summary.tsv`, a group is counted as found only if the highest reported identity for that group is **strictly greater than 80%**. Thus, exactly 80% is counted as missing in that file. Change the threshold, if needed, with (for example) `--genome-min-identity 90`.

The optional `--min-identity` and `--min-coverage` filters are separate. Both default to zero, so any BLAST-reported hit is eligible for the group-level calls and summaries. Use these filters only if you intend to require them:

```bash
python3 summarize_cgmlst_by_group.py \
  --alleles cgmlst-centroid.fasta \
  --results cgmlst_results \
  --out cgmlst_summary \
  --min-identity 90 \
  --min-coverage 90 \
  --genome-min-identity 80
```

If filters are added, they change the group-level calls and some output counts. `genome_summary.tsv` and `genome_missing_groups.tsv` use `--genome-min-identity` and the highest reported identity per group; their threshold is independent of the optional group-level filters.

## 6. Output files

The script writes 11 tab-separated files into `cgmlst_summary/`. Re-running it overwrites files with the same names.

| File | Contents |
|---|---|
| `allele_group_summary.tsv` | One row per cgMLST group: number of reference sequences, genomes with a qualifying reported hit, genomes with no reported hit, and identity range. With default zero filters, any reported hit counts here. |
| `sequence_variant_summary.tsv` | One row per curated FASTA sequence, with its group, allele ID, number of genomes with a reported hit to that sequence, and identity range. |
| `variants_no_hits_any_genome.tsv` | Curated reference sequences with no reported hit in any genome TSV. |
| `genome_allele_calls.tsv` | The selected reference hit for each genome and group, with identity and subject coverage. Selection prioritizes coverage, then identity, then bit score. |
| `incomplete_allele_groups.tsv` | Reported group hits excluded by `--min-identity` or `--min-coverage`. With both filters at zero, it should have no data rows. |
| `no_hit_allele_groups.tsv` | Genome–group pairs where no BLAST hit was reported to any reference sequence in the group. |
| `genome_summary.tsv` | One row per genome: number and proportion of 330 groups whose highest reported identity is **greater than** `--genome-min-identity` (default 80). Lowest proportions first. |
| `genome_missing_groups.tsv` | Genome–group pairs with no hit above the genome identity threshold. Includes best reported identity if a hit exists; blank means no hit. Exactly 80% is included as missing at the default threshold. |
| `genome_missing_group_variants.tsv` | The same threshold-missing pairs, with the highest-identity reference sequence ID and curated allele ID when a hit was reported, or `NO_HIT`. |
| `genome_low_identity_or_absent.tsv` | Compact genome–group list where best identity is **below** the genome identity threshold or there is no hit. Unlike `genome_missing_groups.tsv`, an identity exactly equal to 80% is not included. |
| `genome_flagged_loci.tsv` | Compact review table for groups with best identity below the threshold or no hit, including all curated reference allele IDs for each group. |

### Main per-genome summary

Example `genome_summary.tsv`:

```text
genome	groups_total	groups_found_identity_gt_80	groups_missing	proportion_found
GCA_000000001	330	328	2	0.9939
```

`proportion_found` is `groups_found_identity_gt_80 / groups_total`. Multiply it by 100 for a percentage. `groups_missing` is `groups_total - groups_found_identity_gt_80`.

### Missing-locus detail

`genome_flagged_loci.tsv` is the compact report for reviewing low-identity or absent loci:

```text
genome	locus	allele_id	identity	coverage	reason	reference_allele_ids
GCA_000000001	NZ_AOXE01000059.1_38	1948337259	79.200	95.000	identity_below_threshold	1948337259;1362204802
GCA_000000001	NZ_AOXE01000059.1_372	absent			no_reported_hit	1362204802;2492322876
```

- `genome`: genome filename without `.tsv`.
- `locus`: cgMLST group prefix from the FASTA ID.
- `allele_id`: curated ID after `|` for the highest-identity reported hit; `absent` if no hit was reported to the group.
- `identity`: percent identity for that hit.
- `coverage`: percent subject-coordinate span for that hit.
- `reason`: `identity_below_threshold` or `no_reported_hit`.
- `reference_allele_ids`: all curated reference variant IDs in the group, separated by semicolons.

The table includes a locus if no reported hit has identity above the threshold, or if BLAST reports no hit. It does not include loci with identity exactly equal to 80% at the default threshold, although `genome_missing_groups.tsv` does include them. Coverage is shown for context but does not determine inclusion. A short hit with identity above 80% is not listed in this flagged file.

For a hit below threshold, `allele_id` identifies the highest-identity **reference hit**, not necessarily the variant actually present in the genome. For `no_reported_hit`, the reference IDs are possible database variants in that group; none was reported as a BLAST hit. This does not establish which allele the genome should contain.

