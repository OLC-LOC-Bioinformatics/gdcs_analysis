#!/usr/bin/env python3
"""Summarize BLAST hits by cgMLST allele group (prefix before '|' in FASTA IDs)."""

import argparse
import csv
from pathlib import Path


def read_fasta_lengths(fasta):
    lengths = {}
    sequence_id = None
    sequence_length = 0

    def save_record(record_id, length):
        if record_id is None:
            return
        if length == 0:
            raise ValueError(f"FASTA record {record_id!r} has no sequence")
        if record_id in lengths:
            raise ValueError(f"Duplicate FASTA ID: {record_id}")
        lengths[record_id] = length

    with open(fasta) as handle:
        for line_number, line in enumerate(handle, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save_record(sequence_id, sequence_length)
                header = line[1:].split()
                if not header:
                    raise ValueError(f"{fasta}:{line_number}: empty FASTA header")
                sequence_id = header[0]
                sequence_length = 0
            else:
                if sequence_id is None:
                    raise ValueError(f"{fasta}:{line_number}: sequence before FASTA header")
                sequence_length += len(line)

    save_record(sequence_id, sequence_length)
    if not lengths:
        raise ValueError(f"No allele sequences found in {fasta}")
    return lengths


def main():
    parser = argparse.ArgumentParser(
        description="Summarize BLAST hits by the cgMLST group prefix before '|' in FASTA IDs."
    )
    parser.add_argument("--alleles", required=True, help="FASTA used to build the BLAST database")
    parser.add_argument("--results", required=True, help="Directory with one BLAST TSV per genome")
    parser.add_argument("--out", required=True, help="Output directory")
    parser.add_argument(
        "--min-identity", type=float, default=0.0,
        help="Minimum percent identity for a complete call (default: 0; no identity cutoff)",
    )
    parser.add_argument(
        "--genome-min-identity", type=float, default=80.0,
        help="Identity used for the per-genome summary: a group counts as found if any "
        "member has a hit above this percent identity (default: 80)",
    )
    parser.add_argument(
        "--min-coverage", type=float, default=0.0,
        help="Minimum allele coverage for a complete call (default: 0; any reported HSP)",
    )
    args = parser.parse_args()

    if not 0 <= args.min_identity <= 100:
        parser.error("--min-identity must be between 0 and 100")
    if not 0 <= args.min_coverage <= 100:
        parser.error("--min-coverage must be between 0 and 100")

    sequence_lengths = read_fasta_lengths(args.alleles)
    groups = {}
    sequence_to_group = {}
    for sequence_id in sequence_lengths:
        group_id, separator, curated_allele_id = sequence_id.partition("|")
        if not separator or not group_id or not curated_allele_id:
            raise ValueError(
                f"FASTA ID {sequence_id!r} must have format group_id|curated_allele_id"
            )
        sequence_to_group[sequence_id] = group_id
        groups.setdefault(group_id, set()).add(sequence_id)
    group_ids = sorted(groups)

    result_files = sorted(Path(args.results).glob("*.tsv"))
    if not result_files:
        parser.error(f"No .tsv files found in {args.results}")

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    genomes = [path.stem for path in result_files]
    best_by_genome = {}
    best_any_by_genome = {}
    max_identity_by_genome = {}
    best_identity_hit_by_genome = {}
    variant_best_hits = {sequence_id: {} for sequence_id in sequence_to_group}

    for path in result_files:
        best_by_genome[path.stem] = {}
        best_any_by_genome[path.stem] = {}
        max_identity_by_genome[path.stem] = {}
        best_identity_hit_by_genome[path.stem] = {}
        with path.open() as handle:
            for line_number, line in enumerate(handle, 1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    raise ValueError(f"{path}:{line_number}: expected 12 BLAST columns")

                allele_id = fields[1]
                if allele_id not in sequence_lengths:
                    raise ValueError(
                        f"{path}:{line_number}: subject ID {allele_id!r} "
                        "is not in the allele FASTA"
                    )

                group_id = sequence_to_group[allele_id]
                allele_length = sequence_lengths[allele_id]
                identity = float(fields[2])
                subject_start, subject_end = int(fields[7]), int(fields[8])
                coverage = 100.0 * (abs(subject_end - subject_start) + 1) / allele_length
                bitscore = float(fields[10])

                # Prefer the most complete alignment, then identity and bit score.
                candidate = (coverage, identity, bitscore, allele_id)
                previous_variant = variant_best_hits[allele_id].get(path.stem)
                if previous_variant is None or candidate > previous_variant:
                    variant_best_hits[allele_id][path.stem] = candidate

                if identity > max_identity_by_genome[path.stem].get(group_id, -1.0):
                    max_identity_by_genome[path.stem][group_id] = identity
                    best_identity_hit_by_genome[path.stem][group_id] = candidate

                previous_any = best_any_by_genome[path.stem].get(group_id)
                if previous_any is None or candidate > previous_any:
                    best_any_by_genome[path.stem][group_id] = candidate

                if identity < args.min_identity or coverage < args.min_coverage:
                    continue

                previous = best_by_genome[path.stem].get(group_id)
                if previous is None or candidate > previous:
                    best_by_genome[path.stem][group_id] = candidate

    variant_genomes = {sequence_id: [] for sequence_id in sequence_to_group}
    for genome in genomes:
        for group_id, candidate in best_by_genome[genome].items():
            variant_genomes[candidate[3]].append((genome, candidate))

    with (outdir / "allele_group_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "allele_group_id", "sequence_ids_in_group", "genomes_complete",
            "genomes_incomplete_or_low_identity", "genomes_no_hit",
            "genomes_total", "percent_complete",
            "min_identity", "max_identity",
        ])
        for group_id in group_ids:
            hits = [
                best_by_genome[genome][group_id]
                for genome in genomes if group_id in best_by_genome[genome]
            ]
            identities = [hit[1] for hit in hits]
            complete = len(hits)
            incomplete = sum(
                group_id in best_any_by_genome[genome]
                and group_id not in best_by_genome[genome]
                for genome in genomes
            )
            no_hit = len(genomes) - complete - incomplete
            writer.writerow([
                group_id,
                len(groups[group_id]),
                complete,
                incomplete,
                no_hit,
                len(genomes),
                f"{100 * complete / len(genomes):.2f}",
                f"{min(identities):.3f}" if identities else "",
                f"{max(identities):.3f}" if identities else "",
            ])

    with (outdir / "sequence_variant_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "allele_group_id", "curated_allele_id", "sequence_id", "genomes_with_hit",
            "min_identity", "max_identity",
        ])
        for sequence_id, group_id in sorted(sequence_to_group.items()):
            curated_allele_id = sequence_id.partition("|")[2]
            identities = [
                candidate[1]
                for candidate in variant_best_hits[sequence_id].values()
            ]
            writer.writerow([
                group_id,
                curated_allele_id,
                sequence_id,
                len(variant_best_hits[sequence_id]),
                f"{min(identities):.3f}" if identities else "",
                f"{max(identities):.3f}" if identities else "",
            ])

    with (outdir / "variants_no_hits_any_genome.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["allele_group_id", "curated_allele_id", "sequence_id"])
        for sequence_id, group_id in sorted(sequence_to_group.items()):
            if not variant_best_hits[sequence_id]:
                writer.writerow([
                    group_id,
                    sequence_id.partition("|")[2],
                    sequence_id,
                ])

    with (outdir / "genome_allele_calls.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "genome", "allele_group_id", "best_matching_sequence_id", "percent_identity",
            "allele_coverage",
        ])
        for genome in genomes:
            for group_id, candidate in sorted(best_by_genome[genome].items()):
                coverage, identity, _, sequence_id = candidate
                writer.writerow([
                    genome, group_id, sequence_id, f"{identity:.3f}", f"{coverage:.3f}",
                ])

    with (outdir / "incomplete_allele_groups.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "allele_group_id", "genome", "best_matching_sequence_id",
            "percent_identity", "allele_coverage", "reason",
        ])
        for genome in genomes:
            for group_id, candidate in sorted(best_any_by_genome[genome].items()):
                if group_id in best_by_genome[genome]:
                    continue
                coverage, identity, _, sequence_id = candidate
                reasons = []
                if identity < args.min_identity:
                    reasons.append("identity_below_threshold")
                if coverage < args.min_coverage:
                    reasons.append("coverage_below_threshold")
                writer.writerow([
                    group_id, genome, sequence_id, f"{identity:.3f}",
                    f"{coverage:.3f}", ";".join(reasons),
                ])

    with (outdir / "no_hit_allele_groups.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["allele_group_id", "genome"])
        for genome in genomes:
            for group_id in group_ids:
                if group_id not in best_any_by_genome[genome]:
                    writer.writerow([group_id, genome])

    threshold = args.genome_min_identity
    genome_rows = []
    for genome in genomes:
        found = [
            g for g in group_ids
            if max_identity_by_genome[genome].get(g, -1.0) > threshold
        ]
        genome_rows.append((genome, len(found)))

    with (outdir / "genome_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "genome", "groups_total", f"groups_found_identity_gt_{threshold:g}",
            "groups_missing", "proportion_found",
        ])
        # Lowest proportion first so the genomes to investigate appear at the top.
        for genome, count in sorted(genome_rows, key=lambda row: (row[1], row[0])):
            writer.writerow([
                genome, len(group_ids), count, len(group_ids) - count,
                f"{count / len(group_ids):.4f}",
            ])

    with (outdir / "genome_missing_groups.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome", "allele_group_id", "best_identity_found"])
        for genome, _ in genome_rows:
            for group_id in group_ids:
                best = max_identity_by_genome[genome].get(group_id)
                if best is None or best <= threshold:
                    writer.writerow([
                        genome, group_id, "" if best is None else f"{best:.3f}",
                    ])

    with (outdir / "genome_missing_group_variants.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "genome", "allele_group_id", "best_matching_sequence_id",
            "curated_allele_id", "percent_identity", "status",
        ])
        for genome, _ in genome_rows:
            for group_id in group_ids:
                best = max_identity_by_genome[genome].get(group_id)
                if best is not None and best > threshold:
                    continue
                candidate = best_identity_hit_by_genome[genome].get(group_id)
                if candidate is None:
                    writer.writerow([genome, group_id, "", "", "", "NO_HIT"])
                else:
                    sequence_id = candidate[3]
                    writer.writerow([
                        genome, group_id, sequence_id,
                        sequence_id.partition("|")[2], f"{candidate[1]:.3f}",
                        "AT_OR_BELOW_THRESHOLD",
                    ])

    with (outdir / "genome_low_identity_or_absent.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["genome", "locus", "allele_id", "identity", "coverage"])
        for genome, _ in genome_rows:
            for group_id in group_ids:
                candidate = best_identity_hit_by_genome[genome].get(group_id)
                if candidate is None:
                    writer.writerow([genome, group_id, "absent", "", ""])
                elif candidate[1] < threshold:
                    coverage, identity, _, sequence_id = candidate
                    writer.writerow([
                        genome, group_id, sequence_id.partition("|")[2],
                        f"{identity:.3f}", f"{coverage:.3f}",
                    ])

    with (outdir / "genome_flagged_loci.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow([
            "genome", "locus", "allele_id", "identity", "coverage", "reason",
            "reference_allele_ids",
        ])
        for genome, _ in genome_rows:
            for group_id in group_ids:
                candidate = best_identity_hit_by_genome[genome].get(group_id)
                if candidate is not None and candidate[1] >= threshold:
                    continue
                reference_allele_ids = ";".join(
                    sequence_id.partition("|")[2]
                    for sequence_id in sorted(groups[group_id])
                )
                if candidate is None:
                    writer.writerow([
                        genome, group_id, "absent", "", "", "no_reported_hit",
                        reference_allele_ids,
                    ])
                    continue
                coverage, identity, _, sequence_id = candidate
                writer.writerow([
                    genome, group_id, sequence_id.partition("|")[2],
                    f"{identity:.3f}", f"{coverage:.3f}", "identity_below_threshold",
                    reference_allele_ids,
                ])


if __name__ == "__main__":
    main()

