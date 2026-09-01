#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "$0")/../.." && pwd)"
data_dir="$project_root/data/1000-genomes/phase3-chr1"
vcf_name="ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
vcf_path="$data_dir/$vcf_name"
output_path="$data_dir/chr1.phase3.maf05.thin250kb.dosage.tsv.gz"
filter_script="$project_root/analysis/helpers/filter-1000g-chr1.py"
python_path="$project_root/.venv-genetics/bin/python"
base_url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

mkdir -p "$data_dir"

download_if_missing() {
  local filename="$1"
  local minimum_bytes="${2:-1}"
  local observed_bytes=0
  if [[ -f "$data_dir/$filename" ]]; then
    observed_bytes="$(wc -c < "$data_dir/$filename")"
  fi
  if [[ "$observed_bytes" -lt "$minimum_bytes" ]]; then
    curl -L --fail --continue-at - \
      --output "$data_dir/$filename" "$base_url/$filename"
  fi
}

download_if_missing "$vcf_name" 1100000000
download_if_missing "$vcf_name.tbi"
download_if_missing "integrated_call_samples_v3.20130502.ALL.panel"
download_if_missing "20140625_related_individuals.txt"

gzip -t "$vcf_path"

if [[ ! -s "$output_path" || "$vcf_path" -nt "$output_path" ]]; then
  if [[ ! -x "$python_path" ]]; then
    python3 -m venv "$project_root/.venv-genetics"
  fi
  if ! "$python_path" -c 'import pysam' >/dev/null 2>&1; then
    "$python_path" -m pip install pysam
  fi
  "$python_path" "$filter_script" \
    --vcf "$vcf_path" \
    --output "$output_path" \
    --min-maf 0.05 \
    --thin-bp 250000
fi

echo "$output_path"
