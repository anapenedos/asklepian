import os
import sys
import csv
import argparse

from readfq import readfq # thanks heng

parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
parser.add_argument("--elandir", required=True)
args = parser.parse_args()

# Check paths
if not os.path.isfile(args.fasta):
    sys.stderr.write("[FAIL] Could not open FASTA %s.\n" % args.fasta)
    sys.exit(1)

if not os.path.isdir(args.elandir):
    sys.stderr.write("[FAIL] Could not open ELANDIR %s.\n" % args.elandir)
    sys.exit(1)

# Open each QC file in the Elan QC directory and determine the best sequence
# for each central_sample_id by tracking the FASTA with the fewest N sites.
best_qc = {}
n_len_discarded = 0
n_suppressed_discarded = 0

# Grab a list of FASTA in the published/latest directory
# It's quite slow to hit them individually with os.path.exists, so let's just
# turn them into a set
available_fasta_fp = set([])
ELANFA_DIR = os.path.join(args.elandir, 'latest', 'fasta')
for path in os.scandir(ELANFA_DIR):
    available_fasta_fp.add(path.name)
sys.stderr.write("[NOTE] %d FASTA counted in latest/fasta\n" % len(available_fasta_fp))

ELANQC_DIR = os.path.join(args.elandir, 'latest', 'qc')
for path in os.scandir(ELANQC_DIR):
    if path.is_file() and path.name.endswith('.qc'):
        central_sample_id, run_name, suffix = path.name.split('.', 2)
        with open(os.path.join(ELANQC_DIR, path)) as qc_fh:
            for row in csv.DictReader(qc_fh, delimiter='\t'):
                num_bases = int(row["num_bases"])
                pc_acgt = float(row["pc_acgt"])
                pc_masked = float(row["pc_masked"])

                # Calculate absolute masked bases
                num_masked = num_bases * (pc_masked/100.0)

                # Remove genomes shorter than 29 Kbp
                if num_bases < 29000:
                    n_len_discarded += 1
                    continue

                # Skip FASTA that have been suppressed (missing from FASTA dir)
                if not row["fasta_path"] in available_fasta_fp:
                    n_suppressed_discarded += 1
                    continue

                if central_sample_id not in best_qc:
                    # If this is the first genome, assume it is the best
                    best_qc[central_sample_id] = [num_masked, row["fasta_path"]]
                else:
                    # If the best has more masked sites than the current QC
                    # swap to the better run
                    if best_qc[central_sample_id][0] > num_masked:
                        best_qc[central_sample_id] = [num_masked, row["fasta_path"]]

sys.stderr.write("[NOTE] %s best sequences found. Non-best sequences discarded (for length: %d, suppressed %d). Writing FASTA.\n" % (len(best_qc), n_len_discarded, n_suppressed_discarded))

# Emit all (central_sample_id, best FASTA filename) pairs to stderr
for central_sample_id in best_qc:
    sys.stderr.write('\t'.join([
        central_sample_id,
        best_qc[central_sample_id][1]
    ]) + '\n')

# Emit all (central_sample_id, best FASTA) pairs, as a new FASTA on stdout
for central_sample_id in best_qc:
    with open(os.path.join(ELANFA_DIR, best_qc[central_sample_id][1])) as fasta_fh:
        for name, seq, qual in readfq(fasta_fh):
            sys.stdout.write('>%s\n%s\n' % (central_sample_id, seq))
            break # Write only the first sequence (the only one handled by Elan)

