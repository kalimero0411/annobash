#!/usr/bin/env python3
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
    return urllib.parse.quote(s, safe="()[],:._+-/")

def clean_field(s: str) -> str:
    s = s.strip()
    return "" if s in ("", "-", "--") else s

def parse_fasta_lengths(genome_fasta: str):
    contig_order = []
    contig_len = {}
    cur = None
    length = 0
    with open(genome_fasta, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
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
    m = {}
    with open(eggnog_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            query = cols[0]
            desc = clean_field(cols[7])
            pref = clean_field(cols[8])
            gos  = clean_field(cols[9])
            ecs  = clean_field(cols[10])

            product = clean_field(pref if pref else desc)
            key = query
            if key_mode == "strip-isoform":
                key = re.sub(r"\.\d+$", "", key)

            m[key] = {
                "product": url_encode_product(product) if product else "",
                "go": gos.replace(" ", ""),
                "ec": ecs.replace(" ", ""),
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
    if mode == "none":
        return set()
    keep = set()
    with open(gff_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "mRNA":
                continue
            mid = ID_RE.search(cols[8])
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default=None, help="Genome FASTA (optional unless --sequence-region).")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--eggnog", required=True)
    ap.add_argument("--out", required=True)

    ap.add_argument("--embed-fasta", action="store_true", help="Append ##FASTA and genome sequence (requires --genome).")
    ap.add_argument("--sequence-region", action="store_true",
                    help="Add ##sequence-region directives AND 'region' features (requires --genome + --taxon-id).")
    ap.add_argument("--taxon-id", default=None,
                    help="NCBI taxon ID (digits). Required if --sequence-region. Will be written as Dbxref=taxon:<id> on region features.")

    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--annotate", choices=["mRNA", "CDS", "both"], default="mRNA")
    ap.add_argument("--filter", choices=["none", "primary-dot1"], default="none")
    ap.add_argument("--id-fallback", choices=["none", "strip-isoform"], default="none")
    ap.add_argument("--report", default=None)
    ap.add_argument("--missing", default=None)
    args = ap.parse_args()

    # Enforce requirements
    if args.sequence_region and not args.genome:
        die("--genome is mandatory when using --sequence-region (region end must come from FASTA).")
    if args.sequence_region and not args.taxon_id:
        die("--taxon-id is mandatory when using --sequence-region (mpwt expects Dbxref=taxon:<id> on region).")
    if args.taxon_id and not re.fullmatch(r"\d+", args.taxon_id):
        die("--taxon-id must be numeric (NCBI taxon id).")

    contig_order, contig_len = ([], {})
    if args.genome:
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

    def write_sequence_region_and_region_features(fout):
        tax = args.taxon_id  # guaranteed present when called
        for sid in contig_order:
            L = int(contig_len.get(sid, 0))
            if L > 0:
                fout.write(f"##sequence-region {sid} 1 {L}\n")
        for sid in contig_order:
            L = int(contig_len.get(sid, 0))
            if L > 0:
                rid = f"region:{sid}"
                # mpwt requirement: ;Dbxref=taxon:taxonid;
                attrs = f"ID={rid};Name={sid};Dbxref=taxon:{tax}"
                fout.write(f"{sid}\tpathologic\tregion\t1\t{L}\t.\t+\t.\t{attrs}\n")

    saw_gff_version = False
    inserted_regions = False

    with open(args.gff, "r", encoding="utf-8") as fin, open(args.out, "w", encoding="utf-8") as fout:
        for line in fin:
            if not line.strip():
                fout.write(line)
                continue
            if line.startswith("##gff-version"):
                saw_gff_version = True
                fout.write(line)
                continue
            if line.startswith("#"):
                fout.write(line)
                continue

            if not inserted_regions:
                if not saw_gff_version:
                    fout.write("##gff-version 3\n")
                    saw_gff_version = True
                if args.sequence_region:
                    write_sequence_region_and_region_features(fout)
                inserted_regions = True

            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9:
                fout.write(line)
                continue

            seqid, source, ftype, start, end, score, strand, phase, attrs = cols

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
                        cols[8] = attrs
                    else:
                        stats[f"no_eggnog_match_{ftype}"] += 1
                        if ftype == "mRNA" and tid:
                            unmatched_mrnas.append(tid)

            fout.write("\t".join(cols) + "\n")

        if not inserted_regions:
            if not saw_gff_version:
                fout.write("##gff-version 3\n")
            if args.sequence_region:
                write_sequence_region_and_region_features(fout)

        if args.embed_fasta:
            if not args.genome:
                print("WARNING: --embed-fasta requested but --genome not provided; skipping FASTA embedding.", file=sys.stderr)
            else:
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

    print(f"Wrote {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
