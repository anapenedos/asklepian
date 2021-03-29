import os
import sys
import csv
import argparse

from readfq import readfq # cheers heng

parser = argparse.ArgumentParser()
parser.add_argument("--fasta", required=True)
parser.add_argument("--meta", required=True)
parser.add_argument("--best-ls", required=True)
args = parser.parse_args()

# Check files exist
for fpt, fp in ("FASTA", args.fasta), ("META", args.meta), ("BEST-LS", args.best_ls):
    if not os.path.isfile(fp):
        sys.stderr.write("[FAIL] Could not open %s %s.\n" % (fpt, fp))
        sys.exit(1)
    else:
        sys.stderr.write("[NOTE] %s: %s\n" % (fpt, fp))


# Map cog to PAG name
best_pags = {}
with open(args.best_ls) as best_fh:
    for line in best_fh:
        cogid, climb_fn, pag_name, new = line.strip().split('\t')
        best_pags[cogid] = pag_name
sys.stderr.write("[NOTE] %d best PAGs loaded\n" % len(best_pags))

# Grab and load the metadata table and extract the required metadata
# NOTE sample_date defined as collection_date else received_date
parsed_metadata = {}
seen_pags = set([])
with open(args.meta) as metadata_fh:
    for row in csv.DictReader(metadata_fh, delimiter='\t'):
        cogid = row["central_sample_id"]
        pag_name = row["published_name"]

        if cogid not in best_pags:
            # Ignore cogs without a best PAG, they will have been omitted
            # e.g. by get_best_ref for being too short
            continue

        if best_pags[cogid] == pag_name:
            seen_pags.add(pag_name)
        else:
            continue

        sample_date = row["collection_date"]

        # Try the received date if collection date is invalid
        if not sample_date or sample_date == "None" or len(sample_date) == 0:
            sample_date = row["received_date"]

        # Give up with an error if impossibly, the sample_date could not be assigned...
        if not sample_date or sample_date == "None" or len(sample_date) == 0:
            sys.stderr.write("[FAIL] No sample date for %s\n" % cogid)
            sys.exit(2)

        parsed_metadata[cogid] = {
            "adm1": row["adm1"],
            "collection_pillar": row["collection_pillar"],
            "collection_or_received_date": sample_date,
            "published_date": row["published_date"],
        }

sys.stderr.write("[NOTE] %d samples with metadata loaded\n" % len(parsed_metadata))
if len(seen_pags) != len(best_pags):
    missing = set(best_pags.values()) - seen_pags
    for pag in missing:
        sys.stderr.write("[WARN] Best PAG found for %s but not matched to metadata\n" % pag)
    sys.exit(3)

# Load the FASTA, lookup and emit the sample_date and genome sequence
print(','.join([
    "COG-ID",
    "Sample_date",
    "Adm1",
    "Pillar",
    "Published_date",
    "Sequence",
]))
with open(args.fasta) as all_fh:
    for name, seq, qual in readfq(all_fh):
        central_sample_id = name

        print(','.join([
            central_sample_id,
            parsed_metadata[central_sample_id]["collection_or_received_date"],
            parsed_metadata[central_sample_id]["adm1"],
            parsed_metadata[central_sample_id]["collection_pillar"],
            parsed_metadata[central_sample_id]["published_date"],
            seq,
        ]))

