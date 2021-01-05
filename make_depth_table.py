import os
import sys
import argparse

import pysam

BAM_DIR = "/cephfs/covid/bham/nicholsz/artifacts/elan2/staging/alignment"
DEPTH_DIR = "/cephfs/covid/bham/nicholsz/artifacts/elan2/staging/depth"

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

    total_depths = {} # use a dict to safeguard against non-consecutive pos in depth file
    with open(os.path.join(DEPTH_DIR, bam.name+'.depth')) as depth_fh:
        for line in depth_fh:
            ref, pos, tot_depth = line.strip().split('\t')
            total_depths[int(pos)] = int(tot_depth)

    # Assume the reference is valid, because its come through swell/elan
    sys.stderr.write("[NOTE] %s\n" % bam.name)
    if len(bam_fh.references) > 1:
        sys.stderr.write("[WARN] %s has %d references.\n" % (bam.name, len(bam_fh.references)))
    ref = bam_fh.references[0]

    depths = bam_fh.count_coverage(
            contig=ref,
            quality_threshold=0,
            read_callback="all"
    )
    n_all = [total_depths.get(i, 0) for i in range(bam_fh.lengths[0])]

    if args.long:
        for i in range(bam_fh.lengths[0]):
            n_a = depths[0][i]
            n_c = depths[1][i]
            n_g = depths[2][i]
            n_t = depths[3][i]
            n_all = total_depths.get(i, 0)
            n_valid = n_a + n_c + n_g + n_t
            print(','.join([str(x) for x in [
                central_sample_id,
                run_name,
                i+1,
                n_all,
                n_valid,
                n_a,
                n_c,
                n_g,
                n_t,
            ]]))

    if args.wide:
        print(','.join([str(x) for x in [
            central_sample_id,
            run_name,
            ":".join([str(x) for x in n_all]),
            ":".join([str(x) for x in depths[0]]),
            ":".join([str(x) for x in depths[1]]),
            ":".join([str(x) for x in depths[2]]),
            ":".join([str(x) for x in depths[3]]),
        ]]))

    bam_fh.close()

missing_cogs = set(best_bams.keys()) - seen_cogs
sys.stderr.write("[WARN] %d best sequences could not be matched to a BAM\n" % len(missing_cogs))

