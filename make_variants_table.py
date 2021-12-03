import os
import sys
import argparse
from readfq import readfq # cheers heng
from dataclasses import dataclass
from typing import Union
import pandas as pd


def check_exist(ref, msa):
    # Check files exist
    for fpt, fp in ("REF", ref), ("MSA", msa):
        if not os.path.isfile(fp):
            raise FileNotFoundError(f"[FAIL] Could not open {fpt} {fp}.")
        else:
            sys.stderr.write("[NOTE] %s: %s\n" % (fpt, fp))


def load_ref_seq(ref):
    # Load the ref and assign it to ref_seq
    with open(ref) as canon_fh:
        for name, seq, qual in readfq(canon_fh):
            break
        if not name:
            raise ValueError("[FAIL] Could not read sequence from reference.")
        else:
            return seq


def get_last_base_call_in_seq(sequence):
    """
    Given a sequence, returns the position of the last base in the sequence
    that is not a deletion ('-').
    """
    last_base_pos = None
    for pos_from_end, base in enumerate(sequence[::-1]):
        if base != '-':
            last_base_pos = len(sequence) - pos_from_end
            break
    return last_base_pos


@dataclass
class SeqComparisonState:
    """
    Class to keep track of the comparison between two sequences.
    The comparison of sequences:
             1 2 3 4 |5 6 7 8 9 10 11 12 13 14 15 16 17
        ref: A T G C |G G C T G A  A  T  T  A  A |G  G
        seq: A C - C |- - - T G A  A  C  T  -  - |-  -
    would add these lines to the output object provided (StringIO or list):
        | Pos | Ref | Seq | Is_indel |
        |:---:|:---:|:---:|:--------:|
        |  2  |  T  |  C  |    0     |
        |  3  |     | 1D  |    1     |
        |  5  |     | 3D  |    1     |
        | 12  |  T  |  C  |    0     |
        | 14  |  A  |  N  |    0     |
        | 15  |  A  |  N  |    0     |
    Pipes in ref and seq mark the start and end positions within which variant
    calls are made. Calls are not made outside that range unless a base has
    been called in seq prior to that position due to low quality sequence.
    Beginning and end of range can be defined with `analyses_start` and
    `analyses_end` attributes.
    """
    seq_name: str
    seq_end: int
    output: Union[str, list]
    curr_pos: int = 0
    analyses_start: int = 256
    analyses_end: int = 29675
    del_len: int = 0
    seq_started: bool = False
    in_variant_call_range: bool = False

    def _emit(self, pos, ref_base, seq_base, is_indel):
        """Write relevant variant call to output given."""
        if isinstance(self.output, str):
            self.output += (
                ','.join(
                    [self.seq_name, str(pos), ref_base, seq_base, str(is_indel)
                     ]))
            self.output += '\n'
        elif isinstance(self.output, list):
            self.output.append(
                [self.seq_name, pos, ref_base, seq_base, is_indel])

    def process_pair(self, ref_base, seq_base):
        """
        Compare the nucleotide in the reference and the sequence being
        compared, keeping track of deletions and differences (not insertions).
        """
        self.curr_pos += 1
        # Handle beginning and end of sequence differently.
        if self.analyses_start <= self.curr_pos <= self.analyses_end:
            self.in_variant_call_range = True
        else:
            self.in_variant_call_range = False

        # Handle deletions.
        if seq_base == '-':
            # if in the range of positions in the sequence for which a base has
            # been called, track deletion length regardless of whether in the
            # variant calling range
            if self.seq_started and self.curr_pos < self.seq_end:
                self.del_len += 1
                return
            # if out of base calls position range, if in the variant call
            # range, emit an 'N' so that variant caller knows this is a low QC
            # seq rather than a real deletion
            elif self.in_variant_call_range:
                self._emit(self.curr_pos, ref_base, 'N', 0)
                return
            # if no base has been called yet and we are out of the variant call
            # range, skip position
            else:
                return
        # get length of previous deletion if base found
        elif self.del_len != 0 and self.seq_started:
            self._emit(
                self.curr_pos - self.del_len, '', f'{self.del_len}D', 1)
            self.del_len = 0
        # Mark sequence started when the first base is called
        else:
            if not self.seq_started:
                self.seq_started = True

        # Handle end-of-sequence marker
        if seq_base is None:
            # We only want the end-of-deletion handling above.
            return

        # Report differences
        if seq_base != ref_base:
            self._emit(self.curr_pos, ref_base, seq_base, 0)


def process_seq(
        central_sample_id, seq, ref_seq, output='',
        first_analysed_nt=256, last_analysed_nt=29675):
    """
    Return variants table for a sequence.

    Parameters
    ----------
    msa : str or pathlib.Path
        Path to a multiple sequence alignment file.
    ref_seq_fp : str or pathlib.Path
        Path to a reference sequence to compare wach sequence in the MSA to.
    output : str or list, default ''
        If '' is passed, returns a csv-like string to be passed to stdout.
        If [] is passed, returns a list of lists to be passed to a dataframe
        (each sub-list is a row).
    first_analysed_nt : int, default 256
        The position of the first base used for variant calling (genome termini
        ignored due to low QC base calls).
    last_analysed_nt : int, default 29675
        Like `first_analysed_nt`, but the last position used for variant
        calling.

    Returns
    -------
    str or list
        A csv-like str object or a list of lists to be passed to the `data`
        keyword argument of a pandas dataframe.
    """
    # create a sequence comparison object
    comparator = SeqComparisonState(
        central_sample_id,
        get_last_base_call_in_seq(seq),
        output,
        analyses_start=first_analysed_nt,
        analyses_end=last_analysed_nt)
    # compare seq to ref
    for ref_base, seq_base in zip(ref_seq, seq):
        comparator.process_pair(ref_base, seq_base)
    # ensure last deletion if any is processed
    comparator.process_pair(None, None)
    return comparator.output


def process_msa_to_cmd_line(
        msa, ref_seq_fp, first_analysed_nt=256, last_analysed_nt=29675):
    """
    Compare sequences in a MSA to a reference sequence and send the variants
    table to a pandas dataframe.

    Parameters
    ----------
    msa : str or pathlib.Path
        Path to a multiple sequence alignment file.
    ref_seq_fp : str or pathlib.Path
        Path to a reference sequence to compare wach sequence in the MSA to.
    first_analysed_nt : int, default 256
        The position of the first base used for variant calling (genome termini
        ignored due to low QC base calls).
    last_analysed_nt : int, default 29675
        Like `first_analysed_nt`, but the last position used for variant
        calling.

    pd.DataFrame
        Pandas dataframe containing the variants table for the MSA.
    """
    try:
        ref_seq = load_ref_seq(ref_seq_fp)
    except ValueError as e:
        sys.stderr.write(f'{e}\n')
        sys.exit(2)

    column_names = ['COG-ID', 'Position', 'Reference_Base', 'Alternate_Base',
                    'Is_Indel']
    sys.stdout.write(','.join(column_names))
    sys.stdout.write('\n')
    with open(msa) as all_fh:
        for central_sample_id, seq, _ in readfq(all_fh):
            output = process_seq(
                central_sample_id, seq, ref_seq, output='',
                first_analysed_nt=first_analysed_nt,
                last_analysed_nt=last_analysed_nt)
            sys.stdout.write(output)


def process_msa_to_dataframe(
        msa, ref_seq_fp, first_analysed_nt=256, last_analysed_nt=29675):
    """
    Compare sequences in a MSA to a reference sequence and send the variants
    table to stdout in csv format.

    Parameters
    ----------
    msa : str or pathlib.Path
        Path to a multiple sequence alignment file.
    ref_seq_fp : str or pathlib.Path
        Path to a reference sequence to compare wach sequence in the MSA to.
    first_analysed_nt : int, default 256
        The position of the first base used for variant calling (genome termini
        ignored due to low QC base calls).
    last_analysed_nt : int, default 29675
        Like `first_analysed_nt`, but the last position used for variant
        calling.
    """
    ref_seq = load_ref_seq(ref_seq_fp)

    column_names = ['COG-ID', 'Position', 'Reference_Base', 'Alternate_Base',
                    'Is_Indel']
    data = []
    with open(msa) as all_fh:
        for central_sample_id, seq, _ in readfq(all_fh):
            data.extend(process_seq(
                central_sample_id, seq, ref_seq, output=[],
                first_analysed_nt=first_analysed_nt,
                last_analysed_nt=last_analysed_nt))
    return pd.DataFrame(columns=column_names, data=data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True)
    parser.add_argument("--msa", required=True)
    args = parser.parse_args()
    try:
        check_exist(args.ref, args.msa)
    except FileNotFoundError as e:
        sys.stderr.write(f'{e}\n')
        sys.exit(1)
    process_msa_to_cmd_line(args.msa, args.ref)
