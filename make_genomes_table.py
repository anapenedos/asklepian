import os
import sys
import csv
import argparse

from readfq import readfq # cheers heng

parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
parser.add_argument("--meta", required=True)
args = parser.parse_args()

# Check files exist
for fpt, fp in ("FASTA", args.fasta), ("META", args.meta):
    if not os.path.isfile(fp):
        sys.stderr.write("[FAIL] Could not open %s %s.\n" % (fpt, fp))
        sys.exit(1)
    else:
        sys.stderr.write("[NOTE] %s: %s\n" % (fpt, fp))


# Grab and load the metadata table and extract the sample_date
# NOTE sample_date defined as collection_date else received_date
collection_or_received_date = {}
with open(args.meta) as metadata_fh:
    for row in csv.DictReader(metadata_fh, delimiter='\t'):
        cogid = row["central_sample_id"]
        sample_date = row["collection_date"]

        # Try the received date if collection date is invalid
        if not sample_date or sample_date == "None" or len(sample_date) == 0:
            sample_date = row["received_date"]

        # Give up with an error if impossibly, the sample_date could not be assigned...
        if not sample_date or sample_date == "None" or len(sample_date) == 0:
            sys.stderr.write("[FAIL] No sample date for %s\n" % cogid)
            sys.exit(2)
        else:
            if cogid not in collection_or_received_date:
                collection_or_received_date[cogid] = sample_date
            else:
                if collection_or_received_date[cogid] != sample_date:
                    sys.stderr.write("[FAIL] Sample date mismatch for %s\n" % cogid)
                    sys.exit(2)

sys.stderr.write("[NOTE] %d samples with sample_date loaded\n" % len(collection_or_received_date))


# Load the FASTA, lookup and emit the sample_date and genome sequence
print(','.join([
    "COG-ID",
    "Sample_date",
    "Sequence",
]))
with open(args.fasta) as all_fh:
    for name, seq, qual in readfq(all_fh):
        central_sample_id = name

        print(','.join([
            central_sample_id,
            collection_or_received_date[central_sample_id],
            seq
        ]))

