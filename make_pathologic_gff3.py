#!/usr/bin/env python3
"""
make_pathologic_gff3.py

Creates a Pathway Tools / PathoLogic-friendly GFF3 by:
  1) reading genome FASTA (contig IDs + lengths),
  2) reading structural GFF3,
  3) reading eggNOG-mapper emapper.annotations,
  4) adding functional attributes to mRNA/CDS (product, Ontology_term, Dbxref=EC:...),
  5) adding required ##sequence-region directives (one per contig).

Example:
  python3 make_pathologic_gff3.py \
    --genome Anastatica_hierochuntica.fasta \
    --gff Anastatica_hierochuntica.gff \
    --eggnog Anastatica_eggnog_mapper.emapper.annotations \
    --out Anastatica_hierochuntica.pwt.gff3 \
    --filter primary-dot1 \
    --annotate mRNA \
    --sequence-region \
    --embed-fasta \
    --missing mrna_without_eggnog.txt \
    --report pwt.report.tsv
"""

import argparse
import sys
import re
import urllib.parse
from collections import defaultdict

ID_RE = re.compile(r'(?:^|;)ID=([^;]+)')
PARENT_RE = re.compile(r'(?:^|;)Parent=([^;]+)')

def die(msg: str, code: int = 1) -> None:
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

def url_encode_product(s: str) -> str:
    # strict GFF3-safe encoding for spaces/special chars
    return urllib.parse.quote(s, safe="()[],:._+-/")

def clean_field(s: str) -> str:
    s = s.strip()
    return "" if s in ("", "-", "--") else s

def parse_fasta_lengths(genome_fasta: str):
    """
    Returns:
      contig_order: list of contig IDs in FASTA order
      contig_len: dict {contig_id: length}
    """
    contig_order = []
    contig_len = {}
    cur = None
    length = 0

    with open(genome_fasta, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                # finalize previous
                if cur is not None:
                    contig_len[cur] = length
                hdr = line[1:].strip()
                if not hdr:
                    continue
                cur = hdr.split()[0]
                contig_order.append(cur)
                length = 0
            else:
                if cur is None:
                    continue
                length += len(line.strip())

    if cur is not None:
        contig_len[cur] = length

    if not contig_order:
        die(f"No FASTA headers found in {genome_fasta}")
    return contig_order, contig_len

def load_eggnog_map(eggnog_path: str, key_mode: str):
    """
    key_mode:
      - exact
      - strip-isoform  (Ah01T00010.2 -> Ah01T00010)
    """
    m = {}
    with open(eggnog_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue

            query = cols[0]
            desc = clean_field(cols[7])   # Description
            pref = clean_field(cols[8])   # Preferred_name
            gos  = clean_field(cols[9])   # GOs
            ecs  = clean_field(cols[10])  # EC

            product = clean_field(pref if pref else desc)

            key = query
            if key_mode == "strip-isoform":
                key = re.sub(r"\.\d+$", "", key)

            m[key] = {
                "product": url_encode_product(product) if product else "",
                "go": gos.replace(" ", ""),
                "ec": ecs.replace(" ", ""),
                "raw_query": query,
            }

    if not m:
        die(f"No usable rows parsed from eggNOG file: {eggnog_path}")
    return m

def has_attr(attrs: str, key: str) -> bool:
    return re.search(r'(?:^|;)' + re.escape(key) + r'=', attrs) is not None

def set_attr(attrs: str, key: str, value: str, overwrite: bool) -> str:
    if not value:
        return attrs
    if has_attr(attrs, key):
        if not overwrite:
            return attrs
        attrs = re.sub(r'((?:^|;)' + re.escape(key) + r'=)[^;]*', r'\1' + value, attrs)
        return attrs
    sep = "" if (attrs == "" or attrs.endswith(";")) else ";"
    return attrs + sep + f"{key}={value}"

def choose_primary_transcripts(gff_path: str, mode: str):
    """
    mode:
      - none
      - primary-dot1
    """
    if mode == "none":
        return set()

    keep = set()
    with open(gff_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9:
                continue
            if p[2] != "mRNA":
                continue
            mid = ID_RE.search(p[8])
            if not mid:
                continue
            tx = mid.group(1)
            if mode == "primary-dot1":
                if tx.endswith(".1"):
                    keep.add(tx)
            else:
                die(f"Unknown filter mode: {mode}")

    return keep

def main():
    ap = argparse.ArgumentParser(description="Create PathoLogic-friendly GFF3 by adding eggNOG GO/EC/product + sequence-region.")
    ap.add_argument("--genome", required=True, help="Genome FASTA")
    ap.add_argument("--gff", required=True, help="Input GFF3")
    ap.add_argument("--eggnog", required=True, help="eggNOG emapper.annotations")
    ap.add_argument("--out", required=True, help="Output GFF3")
    ap.add_argument("--embed-fasta", action="store_true", help="Append ##FASTA and genome sequence at end of GFF3")
    ap.add_argument("--overwrite", action="store_true", help="Overwrite existing product/Ontology_term/Dbxref")
    ap.add_argument("--annotate", choices=["mRNA", "CDS", "both"], default="mRNA",
                    help="Where to write functional attributes (default: mRNA)")
    ap.add_argument("--filter", choices=["none", "primary-dot1"], default="none",
                    help="Filter transcripts (default: none)")
    ap.add_argument("--warn-missing-contigs", action="store_true",
                    help="Warn if any GFF3 seqid not present in genome FASTA")
    ap.add_argument("--report", default=None, help="Write TSV metrics report")
    ap.add_argument("--missing", default=None, help="Write unmatched mRNA IDs (no eggNOG match) one-per-line")
    ap.add_argument("--id-fallback", choices=["none", "strip-isoform"], default="none",
                    help="Fallback for eggNOG query IDs (default: none)")
    ap.add_argument("--sequence-region", action="store_true",
                    help="Add ##sequence-region <seqid> 1 <len> directives (required by PathoLogic)")
    args = ap.parse_args()

    contig_order, contig_len = parse_fasta_lengths(args.genome)

    egg_key_mode = "strip-isoform" if args.id_fallback == "strip-isoform" else "exact"
    egg = load_eggnog_map(args.eggnog, key_mode=egg_key_mode)

    keep_tx = choose_primary_transcripts(args.gff, args.filter)
    stats = defaultdict(int)
    unmatched_mrnas = []

    def egg_lookup(tid: str):
        rec = egg.get(tid)
        if rec:
            return rec
        if args.id_fallback == "strip-isoform":
            base = re.sub(r"\.\d+$", "", tid)
            return egg.get(base)
        return None

    # We insert ##sequence-region lines once, right after the initial GFF3 directives/comments
    inserted_seq_region = False
    saw_gff_version = False

    with open(args.gff, "r", encoding="utf-8") as fin, open(args.out, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                fout.write(line)
                continue

            if line.startswith("##gff-version"):
                saw_gff_version = True
                fout.write(line)
                continue

            # Write other directives/comments, but before the first feature line,
            # inject ##sequence-region if requested and not yet inserted.
            if line.startswith("#"):
                fout.write(line)
                continue

            # First non-comment feature line reached
            if not inserted_seq_region:
                # Ensure gff-version exists at top if missing
                if not saw_gff_version:
                    fout.write("##gff-version 3\n")
                    saw_gff_version = True

                if args.sequence_region:
                    for cid in contig_order:
                        L = contig_len.get(cid, 0)
                        if L > 0:
                            fout.write(f"##sequence-region {cid} 1 {L}\n")
                inserted_seq_region = True

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                fout.write(line)
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = parts

            if args.warn_missing_contigs and seqid not in contig_len:
                stats["warn_missing_contig_lines"] += 1

            # filtering
            if args.filter != "none":
                if ftype == "gene":
                    pass
                elif ftype == "mRNA":
                    mid = ID_RE.search(attrs)
                    if not mid:
                        stats["dropped_mRNA_no_id"] += 1
                        continue
                    tx = mid.group(1)
                    if tx not in keep_tx:
                        stats["dropped_mRNA_not_primary"] += 1
                        continue
                else:
                    par = PARENT_RE.search(attrs)
                    if not par:
                        stats["dropped_child_no_parent"] += 1
                        continue
                    if par.group(1) not in keep_tx:
                        stats["dropped_child_not_primary"] += 1
                        continue

            # annotate
            if ftype in ("mRNA", "CDS"):
                do_mrna = (args.annotate in ("mRNA", "both") and ftype == "mRNA")
                do_cds  = (args.annotate in ("CDS", "both") and ftype == "CDS")

                if do_mrna or do_cds:
                    tid = None
                    if ftype == "mRNA":
                        mid = ID_RE.search(attrs)
                        tid = mid.group(1) if mid else None
                    else:
                        par = PARENT_RE.search(attrs)
                        tid = par.group(1) if par else None

                    rec = egg_lookup(tid) if tid else None
                    if rec:
                        before = attrs

                        attrs = set_attr(attrs, "product", rec["product"], args.overwrite)
                        if rec["go"]:
                            attrs = set_attr(attrs, "Ontology_term", rec["go"], args.overwrite)
                        if rec["ec"]:
                            ecs = ",".join([f"EC:{x}" for x in rec["ec"].split(",") if x])
                            if ecs:
                                attrs = set_attr(attrs, "Dbxref", ecs, args.overwrite)

                        if attrs != before:
                            stats[f"annotated_{ftype}"] += 1

                        parts[8] = attrs
                        line = "\t".join(parts) + "\n"
                    else:
                        stats[f"no_eggnog_match_{ftype}"] += 1
                        if ftype == "mRNA" and tid:
                            unmatched_mrnas.append(tid)

            fout.write(line)

        # If the file had no features (only comments), still insert directives if asked
        if not inserted_seq_region:
            if not saw_gff_version:
                fout.write("##gff-version 3\n")
            if args.sequence_region:
                for cid in contig_order:
                    L = contig_len.get(cid, 0)
                    if L > 0:
                        fout.write(f"##sequence-region {cid} 1 {L}\n")

        if args.embed_fasta:
            fout.write("##FASTA\n")
            with open(args.genome, "r", encoding="utf-8") as gf:
                for gl in gf:
                    fout.write(gl)

    if args.missing:
        with open(args.missing, "w", encoding="utf-8") as m:
            for tid in sorted(set(unmatched_mrnas)):
                m.write(tid + "\n")

    if args.report:
        with open(args.report, "w", encoding="utf-8") as r:
            r.write("metric\tvalue\n")
            for k in sorted(stats.keys()):
                r.write(f"{k}\t{stats[k]}\n")

    if args.warn_missing_contigs and stats["warn_missing_contig_lines"] > 0:
        print(f"WARNING: {stats['warn_missing_contig_lines']} feature lines had seqid not in genome FASTA", file=sys.stderr)

    print(
        f"Wrote {args.out}\n"
        f"Annotated mRNA: {stats.get('annotated_mRNA',0)}\n"
        f"Annotated CDS: {stats.get('annotated_CDS',0)}\n"
        f"mRNA without eggNOG: {stats.get('no_eggnog_match_mRNA',0)}",
        file=sys.stderr
    )

if __name__ == "__main__":
    main()
