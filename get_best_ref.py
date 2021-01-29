import os
import sys
import csv
import argparse

from readfq import readfq # thanks heng

parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
parser.add_argument("--metrics", required=True)
parser.add_argument("--latest", required=False)
parser.add_argument("--out-ls", required=True)
args = parser.parse_args()

# Check paths
if not os.path.isfile(args.fasta):
    sys.stderr.write("[FAIL] Could not open FASTA %s.\n" % args.fasta)
    sys.exit(1)

if not os.path.isfile(args.metrics):
    sys.stderr.write("[FAIL] Could not open Ocarina metrics output %s.\n" % args.metrics)
    sys.exit(1)

# Keep track of the last best refs so we can flag if the best ref for a sample
# has changed today
previous_best = {}
if args.latest:
    if not os.path.isfile(args.latest):
        sys.stderr.write("[FAIL] Could not previous best ref list %s.\n" % args.latest)
        sys.exit(1)

    with open(args.latest) as latest_fh:
        for line in latest_fh:
            if line[0] == '[' or line[0] == '#':
                continue
            fields = line.strip().split('\t')
            previous_best[ fields[0] ] = fields[1]

best_qc = {}
n_len_discarded = 0

# Use the metrics from Majora via Ocarina to determine the best sequence
# for each central_sample_id by tracking the FASTA with the fewest N sites.
with open(args.metrics) as ocarina_out_fh:
    for row in csv.DictReader(ocarina_out_fh, delimiter='\t'):
        fasta_path = os.path.basename(row["fasta_path"])
        central_sample_id = row["central_sample_id"]
        num_bases = int(row["num_bases"])
        pc_acgt = float(row["pc_acgt"])
        pc_masked = float(row["pc_masked"])

        # Calculate absolute masked bases
        num_masked = num_bases * (pc_masked/100.0)

        # Remove genomes shorter than 29 Kbp
        if num_bases < 29000:
            n_len_discarded += 1
            continue

        if central_sample_id not in best_qc:
            # If this is the first genome, assume it is the best
            best_qc[central_sample_id] = [num_masked, row["published_name"], fasta_path]
        else:
            # If the best has more masked sites than the current QC
            # swap to the better run
            if best_qc[central_sample_id][0] > num_masked:
                best_qc[central_sample_id] = [num_masked, row["published_name"], fasta_path]
            elif best_qc[central_sample_id][0] == num_masked:
                # Use run_name lexo to break tie
                if row["run_name"] > best_qc[central_sample_id][1].split(':')[1]:
                    best_qc[central_sample_id] = [num_masked, row["published_name"], fasta_path]

sys.stderr.write("[NOTE] %s best sequences found. Non-best sequences discarded (for length: %d). Writing FASTA.\n" % (len(best_qc), n_len_discarded))

# Emit all (central_sample_id, best FASTA filename) pairs to stderr
best_published_names = set([])
with open(args.out_ls, 'w') as out_ls_fh:
    for central_sample_id in best_qc:
        status = 1 # assume new
        if best_qc[central_sample_id][2] == previous_best.get(central_sample_id):
            status = 0 # unless new best ref matches last best ref

        out_ls_fh.write('\t'.join([
            central_sample_id,
            best_qc[central_sample_id][2],
            best_qc[central_sample_id][1],
            str(status),
        ]) + '\n')

        # Add this sequence's PAG to the best_published_names set
        best_published_names.add( best_qc[central_sample_id][1] )

# Iterate the matched FASTA and print out sequences that have a name in the best_published_names set
seen_best_published_names = set([])
with open(args.fasta) as latest_fasta_fh:
    for name, seq, qual in readfq(latest_fasta_fh):
        # Apparently I write the names out wrong so that's good
        curr_pag = name.split('|')[0].replace('COGUK', 'COG-UK')
        central_sample_id = curr_pag.split('/')[1]
        if curr_pag in best_published_names:
            sys.stdout.write('>%s\n%s\n' % (central_sample_id, seq))
            seen_best_published_names.add(curr_pag)
sys.stderr.write("[NOTE] %s best sequences written.\n" % len(seen_best_published_names))
sys.stderr.write("[NOTE] %s best sequences missing.\n" % (len(best_published_names) - len(seen_best_published_names)))

if len(seen_best_published_names) != len(best_published_names):
    missing = best_published_names - seen_best_published_names
    for pag in missing:
        sys.stderr.write("[WARN] Best sequence found for %s but not written\n" % pag)
