#!/bin/bash

set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
wrapper="${repo_dir}/modules/nf-demultiplex/bin/souporcell.sh"
test_dir=$(mktemp -d)
trap 'rm -rf "${test_dir}"' EXIT

mkdir "${test_dir}/bin"
printf 'BC1\nBC2\n' | gzip > "${test_dir}/barcodes.tsv.gz"
: > "${test_dir}/input.bam"
: > "${test_dir}/variants.vcf"
: > "${test_dir}/reference.fa"

cat > "${test_dir}/bin/souporcell_pipeline.py" <<'EOF'
#!/bin/bash
printf '%s\n' "$@" > captured_args.txt
EOF
chmod +x "${test_dir}/bin/souporcell_pipeline.py"

run_case() {
  local sample_id=$1
  local skip_remap=$2
  local no_umi=$3

  (
    cd "${test_dir}"
    PATH="${test_dir}/bin:${PATH}" "${wrapper}" \
      "${sample_id}" barcodes.tsv.gz input.bam false variants.vcf reference.fa 2 \
      "${skip_remap}" "${no_umi}" > "${sample_id}.log"
    mv captured_args.txt "${sample_id}.args"
  )
}

run_case remap false false
if grep -qx -- '--skip_remap' "${test_dir}/remap.args"; then
  echo "FAIL: false skip_remap emitted --skip_remap" >&2
  exit 1
fi
grep -qx -- '--no_umi' "${test_dir}/remap.args"
grep -qx -- 'False' "${test_dir}/remap.args"
grep -qx -- 'Number of barcodes: 2' "${test_dir}/remap.log"

run_case skip true true
grep -qx -- '--skip_remap' "${test_dir}/skip.args"
grep -qx -- 'True' "${test_dir}/skip.args"
grep -qx -- '--no_umi' "${test_dir}/skip.args"

echo "PASS: Souporcell boolean arguments"
