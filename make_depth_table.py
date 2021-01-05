import os
import sys
import argparse

import pysam

BAM_DIR = "/cephfs/covid/bham/nicholsz/artifacts/elan2/staging/alignment"

parser = argparse.ArgumentParser()
parser.add_argument("--bestls", required=True)
parser.add_argument("--query", "-q", required=False, default=None)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("--long", action="store_true")
group.add_argument("--wide", action="store_true")
args = parser.parse_args()

# Load best
best_bams = {}
with open(args.bestls) as bestls_fh:
    for line in bestls_fh:
        if line[0] == '#' or line[0] == '[':
            continue
        fields = line.strip().split('\t')
        central_sample_id, run_name, suffix = fields[1].split('.', 2)
        best_bams[central_sample_id] = run_name

seen_cogs = set([])
for bam in os.scandir(BAM_DIR):
    if not bam.is_file() or not bam.name.endswith("bam"):
        continue

    central_sample_id, run_name, suffix = bam.name.split('.', 2)
    if central_sample_id not in best_bams:
        # Likely a low quality sequence that has no best sequence
        #sys.stderr.write("[WARN] %s has unknown central_sample_id\n" % bam.name)
        pass
    else:
        if best_bams[central_sample_id] != run_name:
            continue
    seen_cogs.add(central_sample_id)

    if args.query:
        if args.query + '.' not in bam.name:
            continue

    bam_fh = pysam.AlignmentFile(bam.path)

    # Assume the reference is valid, because its come through swell/elan
    sys.stderr.write("[NOTE] %s\n" % bam.name)
    if len(bam_fh.references) > 1:
        sys.stderr.write("[WARN] %s has %d references.\n" % (bam.name, len(bam_fh.references)))
    ref = bam_fh.references[0]

    # Init counters
    # Let's do it with a dict to start with
    depths = {}
    for i0 in range(bam_fh.lengths[0]):
        depths[i0+1] = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0,
            'N': 0,
            '-': 0,
        }

    # Count
    # NOTE We'll use a read iterator for speed as the pileup usually turns out slower
    # pysam count_coverage works in the same way but does not support counting indels currently
    for read in bam_fh.fetch(contig=ref):
        for qry_pos0, ref_pos0 in read.get_aligned_pairs():
            if qry_pos0 is None and ref_pos0 is None:
                # Ignore clipping
                continue

            if ref_pos0 is None:
                # Ignore insertions as we're only SNP/del aware for the var table
                continue

            ref_pos1 = ref_pos0 + 1

            if not qry_pos0:
                # Deletion (no query against ref)
                depths[ref_pos1]['-'] += 1
            else:
                # Match or substitution
                depths[ref_pos1][read.seq[qry_pos0]] += 1

    # Print
    if args.long:
        for i0 in range(bam_fh.lengths[0]):
            i1 = i0 + 1
            n_a = depths[i1]['A']
            n_c = depths[i1]['C']
            n_g = depths[i1]['G']
            n_t = depths[i1]['T']
            n_n = depths[i1]['N']
            n_del = depths[i1]['-']

            n_all = sum(depths[i1].values())
            n_acgt = n_a + n_c + n_g + n_t
            print(','.join([str(x) for x in [
                central_sample_id,
                run_name,
                i1,
                n_all,
                n_acgt,
                n_a,
                n_c,
                n_g,
                n_t,
                n_n,
                n_del,
            ]]))

    elif args.wide:

        n_all = []
        n_a = []
        n_c = []
        n_g = []
        n_t = []
        n_n = []
        n_del = []

        for i0 in range(bam_fh.lengths[0]):
            i1 = i0 + 1
            n_a.append( depths[i1]['A'] )
            n_c.append( depths[i1]['C'] )
            n_g.append( depths[i1]['G'] )
            n_t.append( depths[i1]['T'] )
            n_n.append( depths[i1]['N'] )
            n_del.append( depths[i1]['-'] )
            n_all.append( sum(depths[i1].values()) )

        print(','.join([str(x) for x in [
            central_sample_id,
            run_name,
            ":".join([str(x) for x in n_all]),
            ":".join([str(x) for x in n_a]),
            ":".join([str(x) for x in n_c]),
            ":".join([str(x) for x in n_g]),
            ":".join([str(x) for x in n_t]),
            ":".join([str(x) for x in n_n]),
            ":".join([str(x) for x in n_del]),
        ]]))

    bam_fh.close()

missing_cogs = set(best_bams.keys()) - seen_cogs
sys.stderr.write("[WARN] %d best sequences could not be matched to a BAM\n" % len(missing_cogs))

