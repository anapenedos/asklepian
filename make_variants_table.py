import os
import sys
import argparse

from readfq import readfq # cheers heng

parser = argparse.ArgumentParser()
parser.add_argument("--ref", required=True)
parser.add_argument("--msa", required=True)
args = parser.parse_args()

# Check files exist
for fpt, fp in ("REF", args.ref), ("MSA", args.msa):
    if not os.path.isfile(fp):
        sys.stderr.write("[FAIL] Could not open %s %s.\n" % (fpt, fp))
        sys.exit(1)
    else:
        sys.stderr.write("[NOTE] %s: %s\n" % (fpt, fp))

# Load the ref and assign it to ref_seq
with open(args.ref) as canon_fh:
    for name, seq, qual in readfq(canon_fh):
        break
    if not name:
        sys.stderr.write("[FAIL] Could not read sequence from reference.\n")
        sys.exit(2)
    else:
        ref_seq = seq

# Open the MSA, iterate over each sequence and walk the genome to find
# diagreements with the loaded reference
# NOTE This particular MSA does not handle insertions
print(','.join([
    "COG-ID",
    "Position",
    "Reference_Base",
    "Alternate_Base",
    "Is_Indel"
]))
with open(args.msa) as all_fh:
    for name, seq, qual in readfq(all_fh):
        central_sample_id = name.split('/')[1]

        query_on_ref_pos = 0
        current_deletion_len = 0

        for qbase in seq:
            if qbase == '-':
                # Extend the length of the current deletion
                current_deletion_len += 1
            else:
                if current_deletion_len > 0:
                    # We've come to the end of a deletion, output it
                    print(','.join([
                        central_sample_id,
                        #str( "%d-%d" % ((query_on_ref_pos-current_deletion_len)+1, query_on_ref_pos) ),
                        str((query_on_ref_pos-current_deletion_len)+1),
                        "",
                        "%dD" % current_deletion_len,
                        "1",
                    ]))
                    current_deletion_len = 0

            # Now deletions are handled, check for single nucleotide variants
            # NOTE This includes missing data such as N
            # NOTE This algorithm does not consider INS against ref
            if qbase != ref_seq[query_on_ref_pos]:
                if current_deletion_len == 0:
                    # SNV detected and we aren't in an active DEL
                    print(','.join([
                        central_sample_id,
                        str(query_on_ref_pos+1),
                        ref_seq[query_on_ref_pos],
                        qbase,
                        "0",
                    ]))

            # Advance pointer (this is overkill here but a useful starting point
            # for a future algo walking the ref for insertions)
            query_on_ref_pos += 1

        if current_deletion_len > 0:
            # Output the last deletion, if there is one 
            # (this is almost always going to be garbage but we include it for completeness)
            print(','.join([
                central_sample_id,
                #str( "%d-%d" % ((query_on_ref_pos-current_deletion_len)+1, query_on_ref_pos) ),
                str((query_on_ref_pos-current_deletion_len)+1),
                "",
                "%dD" % current_deletion_len,
                "1",
            ]))

