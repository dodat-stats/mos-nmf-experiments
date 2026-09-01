#!/usr/bin/env python3
"""Create a compact chromosome-1 dosage table through the VCF tabix index."""

import argparse
import csv
import gzip
import os
import tempfile

import pysam


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-maf", type=float, default=0.05)
    parser.add_argument("--thin-bp", type=int, default=250_000)
    parser.add_argument("--chromosome", default="1")
    return parser.parse_args()


def is_pass_biallelic_snp(record):
    passes_filter = not record.filter.keys() or "PASS" in record.filter.keys()
    return (
        passes_filter
        and len(record.ref) == 1
        and record.alts is not None
        and len(record.alts) == 1
        and len(record.alts[0]) == 1
    )


def alternate_frequency(record):
    value = record.info.get("AF")
    if value is None:
        return None
    if isinstance(value, tuple):
        if len(value) != 1:
            return None
        value = value[0]
    return float(value)


def dosage(sample):
    genotype = sample.get("GT")
    if genotype is None or len(genotype) != 2 or None in genotype:
        return "NA"
    return str(int(genotype[0]) + int(genotype[1]))


def main():
    args = parse_args()
    vcf = pysam.VariantFile(args.vcf)
    samples = list(vcf.header.samples)
    chromosome_length = vcf.header.contigs[args.chromosome].length
    if chromosome_length is None:
        chromosome_length = 249_250_621

    output_directory = os.path.dirname(os.path.abspath(args.output))
    os.makedirs(output_directory, exist_ok=True)
    descriptor, temporary_path = tempfile.mkstemp(
        prefix="chr1-dosage-", suffix=".tsv.gz", dir=output_directory
    )
    os.close(descriptor)

    retained = 0
    try:
        with gzip.open(temporary_path, "wt", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["variant", "chrom", "position", "ref", "alt", "alt_frequency"]
                + samples
            )

            for start in range(0, chromosome_length, args.thin_bp):
                end = min(start + args.thin_bp, chromosome_length)
                midpoint = 0.5 * (start + end)
                best_record = None
                best_af = None
                best_distance = float("inf")

                for record in vcf.fetch(args.chromosome, start, end):
                    if not is_pass_biallelic_snp(record):
                        continue
                    af = alternate_frequency(record)
                    if af is None:
                        continue
                    maf = min(af, 1.0 - af)
                    distance = abs(record.pos - midpoint)
                    if maf < args.min_maf or distance >= best_distance:
                        continue
                    best_record = record.copy()
                    best_af = af
                    best_distance = distance

                if best_record is None:
                    continue

                variant_id = best_record.id
                if variant_id is None:
                    variant_id = (
                        f"{best_record.chrom}:{best_record.pos}:"
                        f"{best_record.ref}:{best_record.alts[0]}"
                    )
                writer.writerow(
                    [
                        variant_id,
                        best_record.chrom,
                        best_record.pos,
                        best_record.ref,
                        best_record.alts[0],
                        f"{best_af:.8g}",
                    ]
                    + [dosage(best_record.samples[sample]) for sample in samples]
                )
                retained += 1

        os.replace(temporary_path, args.output)
    except BaseException:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)
        raise
    finally:
        vcf.close()

    print(f"Retained {retained} SNPs in {args.output}")


if __name__ == "__main__":
    main()
