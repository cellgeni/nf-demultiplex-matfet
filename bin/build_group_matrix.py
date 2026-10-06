#!/usr/bin/env python3
"""Harmonise per-library Souporcell/Vartrix allele matrices into group matrices.

Run by the BUILD_GROUP_MATRIX process; written and verified during validation
(Gate 1 of the validation plan; see docs/pipeline_design.md). It implements the group-matrix
technical contract: union of covered variants in reference-panel order,
concatenated cells in samplesheet order, and explicit provenance/QC metadata.

Standard library only, single streaming pass per library, so it is safe for
large real inputs under LSF.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field

MM_HEADER = "%%MatrixMarket matrix coordinate real general"

# Fallback contig order, used only when no panel header is available.
DEFAULT_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


class ContractError(Exception):
    """A pass/fail violation of the group-matrix contract."""


@dataclass
class Library:
    group_id: str
    sample_id: str
    sample_role: str
    expected_genotypes: int
    souporcell_dir: str
    keys: list[tuple[str, int, str, str]] = field(default_factory=list)
    vcf_records: list[str] = field(default_factory=list)
    barcodes: list[str] = field(default_factory=list)
    n_variants: int = 0
    n_cells: int = 0
    n_entries: int = 0
    alt_sum: int = 0
    ref_sum: int = 0
    first_col: int = 0
    last_col: int = 0

    def path(self, name: str) -> str:
        return os.path.join(self.souporcell_dir, name)


# --------------------------------------------------------------------------
# input reading
# --------------------------------------------------------------------------

def read_groups_csv(path: str, group_id: str) -> list[Library]:
    required = {"group_id", "sample_id", "sample_role", "expected_genotypes", "souporcell_dir"}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ContractError(f"{path}: missing columns {sorted(missing)}")
        rows = [row for row in reader if row["group_id"].strip() == group_id]

    if not rows:
        raise ContractError(f"{path}: no libraries for group_id {group_id!r}")

    libraries: list[Library] = []
    seen: set[str] = set()
    expected = {row["expected_genotypes"].strip() for row in rows}
    if len(expected) != 1:
        raise ContractError(f"group {group_id}: conflicting expected_genotypes {sorted(expected)}")

    for row in rows:
        sample_id = row["sample_id"].strip()
        if sample_id in seen:
            raise ContractError(f"group {group_id}: duplicate sample_id {sample_id!r}")
        seen.add(sample_id)
        role = row["sample_role"].strip()
        if role not in {"mixture", "maternal_reference"}:
            raise ContractError(f"{sample_id}: unknown sample_role {role!r}")
        libraries.append(
            Library(
                group_id=group_id,
                sample_id=sample_id,
                sample_role=role,
                expected_genotypes=int(row["expected_genotypes"]),
                souporcell_dir=os.path.abspath(os.path.expanduser(row["souporcell_dir"].strip())),
            )
        )
    return libraries


def read_contig_order(panel_vcf: str | None) -> list[str]:
    """Contig order from the reference panel header (header read only)."""
    if not panel_vcf:
        return list(DEFAULT_CONTIGS)
    contigs: list[str] = []
    with open(panel_vcf) as handle:
        for line in handle:
            if line.startswith("##contig=<ID="):
                contigs.append(line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">\n").strip())
            elif not line.startswith("##"):
                break
    if not contigs:
        raise ContractError(f"{panel_vcf}: no ##contig header lines found")
    return contigs


def read_vcf(library: Library) -> None:
    """Read the library's covered VCF, validating the variant-key contract."""
    path = library.path("common_variants_covered.vcf")
    keys: list[tuple[str, int, str, str]] = []
    records: list[str] = []
    seen_keys: set[tuple[str, int, str, str]] = set()
    seen_positions: dict[tuple[str, int], tuple[str, str]] = {}

    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ContractError(f"{path}:{lineno}: fewer than 8 VCF columns")
            chrom, pos_text, ref, alt = fields[0], fields[1], fields[3], fields[4]
            if "," in alt:
                raise ContractError(f"{path}:{lineno}: multiallelic record {chrom}:{pos_text} ALT={alt}")
            try:
                pos = int(pos_text)
            except ValueError as exc:
                raise ContractError(f"{path}:{lineno}: non-integer POS {pos_text!r}") from exc
            key = (chrom, pos, ref.upper(), alt.upper())
            if key in seen_keys:
                raise ContractError(f"{path}:{lineno}: duplicate variant key {key}")
            previous = seen_positions.get((chrom, pos))
            if previous is not None and previous != (key[2], key[3]):
                raise ContractError(
                    f"{path}:{lineno}: conflicting alleles at {chrom}:{pos} "
                    f"({previous} vs {(key[2], key[3])})"
                )
            seen_keys.add(key)
            seen_positions[(chrom, pos)] = (key[2], key[3])
            keys.append(key)
            # store the 8 fixed VCF columns as one line to keep large groups light
            records.append("\t".join(fields[:8]))

    library.keys = keys
    library.vcf_records = records
    library.n_variants = len(keys)


def read_barcodes(library: Library) -> None:
    path = library.path("barcodes.tsv")
    with open(path) as handle:
        barcodes = [line.strip() for line in handle if line.strip()]
    if len(set(barcodes)) != len(barcodes):
        raise ContractError(f"{path}: duplicate barcodes within library")
    library.barcodes = barcodes
    library.n_cells = len(barcodes)


def read_mtx_header(path: str) -> tuple[int, int, int]:
    with open(path) as handle:
        first = handle.readline()
        if not first.startswith("%%MatrixMarket matrix coordinate"):
            raise ContractError(f"{path}: unexpected Matrix Market banner {first.strip()!r}")
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                raise ContractError(f"{path}: malformed dimension line {line.strip()!r}")
            return int(parts[0]), int(parts[1]), int(parts[2])
    raise ContractError(f"{path}: no dimension line")


def open_mtx_body(path: str):
    """Open a Matrix Market file positioned on its first entry line."""
    handle = open(path)
    handle.readline()
    for line in handle:
        if line.startswith("%"):
            continue
        parts = line.split()
        if len(parts) != 3:
            handle.close()
            raise ContractError(f"{path}: malformed dimension line {line.strip()!r}")
        return handle, (int(parts[0]), int(parts[1]), int(parts[2]))
    handle.close()
    raise ContractError(f"{path}: no dimension line")


def parse_count(token: str, path: str, lineno: int) -> int:
    """Vartrix declares 'real' but allele counts must be non-negative integers."""
    try:
        value = int(token)
    except ValueError:
        try:
            as_float = float(token)
        except ValueError as exc:
            raise ContractError(f"{path}:{lineno}: non-numeric value {token!r}") from exc
        if as_float != int(as_float):
            raise ContractError(f"{path}:{lineno}: non-integer allele count {token!r}")
        value = int(as_float)
    if value < 0:
        raise ContractError(f"{path}:{lineno}: negative allele count {token!r}")
    return value


# --------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------

def build_union(libraries: list[Library], contigs: list[str]) -> list[tuple[str, int, str, str]]:
    rank = {name: index for index, name in enumerate(contigs)}
    union: dict[tuple[str, int, str, str], None] = {}
    position_alleles: dict[tuple[str, int], tuple[str, str]] = {}

    for library in libraries:
        for key in library.keys:
            chrom, pos, ref, alt = key
            if chrom not in rank:
                raise ContractError(
                    f"{library.sample_id}: contig {chrom!r} is absent from the reference panel order"
                )
            previous = position_alleles.get((chrom, pos))
            if previous is not None and previous != (ref, alt):
                raise ContractError(
                    f"group {library.group_id}: conflicting alleles at {chrom}:{pos} "
                    f"({previous} vs {(ref, alt)}) between libraries"
                )
            position_alleles[(chrom, pos)] = (ref, alt)
            union[key] = None

    return sorted(union, key=lambda key: (rank[key[0]], key[1], key[2], key[3]))


def write_group_vcf(path: str, union: list[tuple[str, int, str, str]],
                    libraries: list[Library], contigs: list[str]) -> None:
    records: dict[tuple[str, int, str, str], str] = {}
    for library in libraries:
        for key, fields in zip(library.keys, library.vcf_records):
            records.setdefault(key, fields)

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d")
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write(f"##fileDate={stamp}\n")
        out.write("##source=build_group_matrix.py\n")
        for contig in contigs:
            out.write(f"##contig=<ID={contig}>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for key in union:
            out.write(records[key] + "\n")


def stream_library(library: Library, rowmap: list[int], col_offset: int,
                   alt_out, ref_out) -> None:
    """Append one library's ALT/REF entries, remapped onto the group axes."""
    alt_path = library.path("alt.mtx")
    ref_path = library.path("ref.mtx")
    alt_in, alt_dims = open_mtx_body(alt_path)
    ref_in, ref_dims = open_mtx_body(ref_path)

    try:
        if alt_dims != ref_dims:
            raise ContractError(
                f"{library.sample_id}: ALT dimensions {alt_dims} differ from REF {ref_dims}"
            )
        n_rows, n_cols, n_entries = alt_dims
        if n_rows != library.n_variants:
            raise ContractError(
                f"{library.sample_id}: matrix rows {n_rows} != VCF records {library.n_variants}"
            )
        if n_cols != library.n_cells:
            raise ContractError(
                f"{library.sample_id}: matrix columns {n_cols} != barcodes {library.n_cells}"
            )

        seen = 0
        alt_total = 0
        ref_total = 0
        for lineno, (alt_line, ref_line) in enumerate(zip(alt_in, ref_in), 1):
            alt_fields = alt_line.split()
            ref_fields = ref_line.split()
            if len(alt_fields) != 3 or len(ref_fields) != 3:
                raise ContractError(f"{alt_path}/{ref_path}: malformed entry at body line {lineno}")
            if alt_fields[0] != ref_fields[0] or alt_fields[1] != ref_fields[1]:
                raise ContractError(
                    f"{library.sample_id}: ALT/REF coordinate mismatch at body line {lineno} "
                    f"({alt_fields[:2]} vs {ref_fields[:2]})"
                )
            row = int(alt_fields[0])
            col = int(alt_fields[1])
            if not 1 <= row <= n_rows:
                raise ContractError(f"{alt_path}:{lineno}: row {row} outside 1..{n_rows}")
            if not 1 <= col <= n_cols:
                raise ContractError(f"{alt_path}:{lineno}: column {col} outside 1..{n_cols}")
            alt_total += parse_count(alt_fields[2], alt_path, lineno)
            ref_total += parse_count(ref_fields[2], ref_path, lineno)

            group_row = rowmap[row - 1]
            group_col = col + col_offset
            # Explicit zeros are part of the Vartrix convention and are preserved.
            alt_out.write(f"{group_row} {group_col} {alt_fields[2]}\n")
            ref_out.write(f"{group_row} {group_col} {ref_fields[2]}\n")
            seen += 1

        if alt_in.readline() or ref_in.readline():
            raise ContractError(f"{library.sample_id}: ALT and REF have different entry counts")
        if seen != n_entries:
            raise ContractError(
                f"{library.sample_id}: declared {n_entries} entries but read {seen}"
            )
    finally:
        alt_in.close()
        ref_in.close()

    library.n_entries = n_entries
    library.alt_sum = alt_total
    library.ref_sum = ref_total


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bgzip_and_index(vcf_path: str, bgzip_cmd: list[str], tabix_cmd: list[str]) -> bool:
    gz_path = vcf_path + ".gz"
    try:
        with open(gz_path, "wb") as out:
            subprocess.run(bgzip_cmd + ["-c", vcf_path], stdout=out, check=True)
        subprocess.run(tabix_cmd + ["-p", "vcf", "-f", gz_path], check=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        if os.path.exists(gz_path):
            os.remove(gz_path)
        print(f"WARNING: could not create {gz_path}: {exc}", file=sys.stderr)
        print(
            "WARNING: the bgzip copy is only for tools that need an indexed VCF; "
            "Souporcell discovery reads variants.vcf directly.",
            file=sys.stderr,
        )
        return False
    return True


def build_group(groups_csv: str, group_id: str, outdir: str, panel_vcf: str | None,
                bgzip_cmd: list[str], tabix_cmd: list[str], do_hash: bool) -> dict:
    libraries = read_groups_csv(groups_csv, group_id)
    contigs = read_contig_order(panel_vcf)

    for library in libraries:
        for name in ("alt.mtx", "ref.mtx", "barcodes.tsv", "common_variants_covered.vcf"):
            if not os.path.isfile(library.path(name)):
                raise ContractError(f"{library.sample_id}: missing input {library.path(name)}")
        read_vcf(library)
        read_barcodes(library)

    union = build_union(libraries, contigs)
    union_index = {key: index + 1 for index, key in enumerate(union)}

    group_dir = os.path.join(outdir, group_id)
    os.makedirs(group_dir, exist_ok=True)

    total_cells = sum(library.n_cells for library in libraries)
    total_entries = sum(read_mtx_header(library.path("alt.mtx"))[2] for library in libraries)

    alt_path = os.path.join(group_dir, "alt.mtx")
    ref_path = os.path.join(group_dir, "ref.mtx")
    with open(alt_path, "w") as alt_out, open(ref_path, "w") as ref_out:
        for out in (alt_out, ref_out):
            out.write(MM_HEADER + "\n")
            out.write(f"% written by build_group_matrix.py for group {group_id}\n")
            out.write(f"{len(union)} {total_cells} {total_entries}\n")

        col_offset = 0
        for library in libraries:
            rowmap = [union_index[key] for key in library.keys]
            library.first_col = col_offset + 1
            library.last_col = col_offset + library.n_cells
            stream_library(library, rowmap, col_offset, alt_out, ref_out)
            col_offset += library.n_cells

    observed_entries = sum(library.n_entries for library in libraries)
    if observed_entries != total_entries:
        raise ContractError(
            f"group {group_id}: entry count changed during the build "
            f"({total_entries} declared, {observed_entries} written)"
        )

    vcf_path = os.path.join(group_dir, "variants.vcf")
    write_group_vcf(vcf_path, union, libraries, contigs)

    group_barcodes: list[str] = []
    with open(os.path.join(group_dir, "barcode_provenance.tsv"), "w") as prov:
        prov.write("group_barcode\tsample_id\traw_barcode\tgroup_column\tsource_column\n")
        column = 0
        for library in libraries:
            for source_column, barcode in enumerate(library.barcodes, 1):
                column += 1
                group_barcode = f"{library.sample_id}::{barcode}"
                group_barcodes.append(group_barcode)
                prov.write(
                    f"{group_barcode}\t{library.sample_id}\t{barcode}\t{column}\t{source_column}\n"
                )

    if len(set(group_barcodes)) != len(group_barcodes):
        raise ContractError(f"group {group_id}: group barcodes are not unique")

    with open(os.path.join(group_dir, "barcodes.tsv"), "w") as out:
        for barcode in group_barcodes:
            out.write(barcode + "\n")

    with open(os.path.join(group_dir, "library_offsets.tsv"), "w") as out:
        out.write("sample_id\tsample_role\tfirst_group_column\tlast_group_column\tn_cells\n")
        for library in libraries:
            out.write(
                f"{library.sample_id}\t{library.sample_role}\t{library.first_col}\t"
                f"{library.last_col}\t{library.n_cells}\n"
            )

    with open(os.path.join(group_dir, "matrix_qc.tsv"), "w") as out:
        out.write("scope\tsample_id\tn_variants\tn_cells\tn_entries\talt_sum\tref_sum\tstatus\n")
        for library in libraries:
            out.write(
                f"library\t{library.sample_id}\t{library.n_variants}\t{library.n_cells}\t"
                f"{library.n_entries}\t{library.alt_sum}\t{library.ref_sum}\tpass\n"
            )
        out.write(
            f"group\t{group_id}\t{len(union)}\t{total_cells}\t{total_entries}\t"
            f"{sum(l.alt_sum for l in libraries)}\t{sum(l.ref_sum for l in libraries)}\tpass\n"
        )

    gz_ok = bgzip_and_index(vcf_path, bgzip_cmd, tabix_cmd)

    manifest = {
        "tool": "build_group_matrix.py",
        "tool_version": "0.1.0",
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "group_id": group_id,
        "expected_genotypes": libraries[0].expected_genotypes,
        "parameters": {
            "groups_csv": os.path.abspath(groups_csv),
            "panel_vcf": os.path.abspath(panel_vcf) if panel_vcf else None,
            "contig_order_source": "panel_header" if panel_vcf else "default",
            "hashed_inputs": do_hash,
            "bgzip_index_written": gz_ok,
        },
        "group": {
            "n_variants": len(union),
            "n_cells": total_cells,
            "n_entries": total_entries,
            "alt_sum": sum(library.alt_sum for library in libraries),
            "ref_sum": sum(library.ref_sum for library in libraries),
        },
        "libraries": [],
        "outputs": {},
    }

    for library in libraries:
        entry = {
            "sample_id": library.sample_id,
            "sample_role": library.sample_role,
            "souporcell_dir": library.souporcell_dir,
            "n_variants": library.n_variants,
            "n_cells": library.n_cells,
            "n_entries": library.n_entries,
            "alt_sum": library.alt_sum,
            "ref_sum": library.ref_sum,
            "first_group_column": library.first_col,
            "last_group_column": library.last_col,
            "inputs": {},
        }
        for name in ("alt.mtx", "ref.mtx", "barcodes.tsv", "common_variants_covered.vcf"):
            path = library.path(name)
            info = {
                "path": path,
                "size": os.path.getsize(path),
                "mtime_utc": datetime.datetime.fromtimestamp(
                    os.path.getmtime(path), datetime.timezone.utc
                ).isoformat(timespec="seconds"),
            }
            if do_hash:
                info["sha256"] = sha256(path)
            entry["inputs"][name] = info
        manifest["libraries"].append(entry)

    for name in sorted(os.listdir(group_dir)):
        if name == "manifest.json":
            continue
        path = os.path.join(group_dir, name)
        info = {"size": os.path.getsize(path)}
        if do_hash:
            info["sha256"] = sha256(path)
        manifest["outputs"][name] = info

    with open(os.path.join(group_dir, "manifest.json"), "w") as out:
        json.dump(manifest, out, indent=2, sort_keys=True)
        out.write("\n")

    return manifest


def index_only(group_id: str, outdir: str, bgzip_cmd: list[str], tabix_cmd: list[str],
               do_hash: bool) -> dict:
    """Add variants.vcf.gz/.tbi to an existing build and refresh its manifest.

    Used when the build ran on a host without a container runtime: the matrices
    are already proven, so only the bgzip derivative and the manifest's record
    of it need to change.
    """
    group_dir = os.path.join(outdir, group_id)
    manifest_path = os.path.join(group_dir, "manifest.json")
    vcf_path = os.path.join(group_dir, "variants.vcf")
    for path in (manifest_path, vcf_path):
        if not os.path.isfile(path):
            raise ContractError(f"{path}: not an existing group build")

    with open(manifest_path) as handle:
        manifest = json.load(handle)

    if not bgzip_and_index(vcf_path, bgzip_cmd, tabix_cmd):
        raise ContractError(f"{group_dir}: bgzip/tabix failed; manifest left unchanged")

    manifest["parameters"]["bgzip_index_written"] = True
    manifest["indexed_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds"
    )
    for name in ("variants.vcf.gz", "variants.vcf.gz.tbi"):
        path = os.path.join(group_dir, name)
        info = {"size": os.path.getsize(path)}
        if do_hash:
            info["sha256"] = sha256(path)
        manifest["outputs"][name] = info

    with open(manifest_path, "w") as out:
        json.dump(manifest, out, indent=2, sort_keys=True)
        out.write("\n")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--groups", required=True, help="groups samplesheet CSV")
    parser.add_argument("--group-id", required=True, help="group_id to build")
    parser.add_argument("--outdir", required=True, help="output root; <outdir>/<group_id> is written")
    parser.add_argument("--panel-vcf", default=None,
                        help="reference panel VCF; only its ##contig header is read, to fix variant order")
    parser.add_argument("--singularity-image", default=None,
                        help="image used to run bgzip/tabix when they are not on PATH")
    parser.add_argument("--singularity-bind", default="/lustre,/nfs",
                        help="paths bound into the image; the defaults are where this project's "
                             "inputs and outputs live")
    parser.add_argument("--index-only", action="store_true",
                        help="do not rebuild; only write variants.vcf.gz/.tbi for an existing "
                             "build and refresh its manifest")
    parser.add_argument("--no-hash", action="store_true",
                        help="skip sha256 of inputs/outputs (faster; manifest then records size and mtime only)")
    args = parser.parse_args(argv)

    if args.singularity_image:
        prefix = ["singularity", "exec", "-B", args.singularity_bind, args.singularity_image]
        bgzip_cmd, tabix_cmd = prefix + ["bgzip"], prefix + ["tabix"]
    else:
        bgzip_cmd = [shutil.which("bgzip") or "bgzip"]
        tabix_cmd = [shutil.which("tabix") or "tabix"]

    if args.index_only:
        try:
            index_only(args.group_id, args.outdir, bgzip_cmd, tabix_cmd, not args.no_hash)
        except ContractError as exc:
            print(f"CONTRACT FAILURE: {exc}", file=sys.stderr)
            return 2
        print(f"indexed {args.group_id}: variants.vcf.gz and .tbi written")
        return 0

    try:
        manifest = build_group(
            groups_csv=args.groups,
            group_id=args.group_id,
            outdir=args.outdir,
            panel_vcf=args.panel_vcf,
            bgzip_cmd=bgzip_cmd,
            tabix_cmd=tabix_cmd,
            do_hash=not args.no_hash,
        )
    except ContractError as exc:
        print(f"CONTRACT FAILURE: {exc}", file=sys.stderr)
        return 2

    group = manifest["group"]
    print(
        f"built group {args.group_id}: {group['n_variants']} variants x {group['n_cells']} cells, "
        f"{group['n_entries']} entries from {len(manifest['libraries'])} libraries"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
