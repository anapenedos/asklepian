#!/usr/bin/bash

# Activate env
eval "$(conda shell.bash hook)"
conda activate samstudio8
source ~/.ocarina
set -euo pipefail

# Init outdir
cd $ASKLEPIAN_DIR
DATESTAMP=${ELAN_DATE} # set by mqtt wrapper
OUTDIR="$ASKLEPIAN_OUTDIR/$DATESTAMP"
mkdir -p $OUTDIR

# Get best ref for each central_sample_id
# TODO There is no need to do this de novo but we do so here for simplicity. We would do well to replace this as it will likely become quite slow, especially with IO perf on login node.
#      We can work toward some future system where Majora keeps a "symlink" to the best QC pag for a PAG Group.
python get_best_ref.py --fasta $COG_PUBLISHED_DIR/elan.latest.consensus.matched.fasta --elandir $COG_PUBLISHED_DIR > $OUTDIR/best_refs.paired.fasta 2> $OUTDIR/best_refs.paired.ls

################################################################################
# This stanza emulates the key first rules from phylopipe 1_preprocess_uk
# that are required to generate a naive insertion-unaware MSA.
# We perform them here independently while we work to stabilise phylopipe.
# See https://github.com/COG-UK/grapevine/blob/master/rules/1_preprocess_uk.smk

# uk_minimap2_to_reference
minimap2 -t 24 -a -x asm5 $WUHAN_FP $OUTDIR/best_refs.paired.fasta > $OUTDIR/output.sam 2> $OUTDIR/mm2.log

# uk_full_untrimmed_alignment
datafunk sam_2_fasta -s $OUTDIR/output.sam -r $WUHAN_FP -o $OUTDIR/naive_msa.fasta 2> $OURDIR/dfunk.log

# The 'naive' MSA here is not filtered on person-level identifiers, neither does
# it consider any filtering or masking. We use it in the knowledge that some
# of the variants may be dirty. Outside of its use for the full variant table here,
# it should not be considered an equivalent drop-in for the phylopipe MSA.
################################################################################

# Make and push genome table
python make_genomes_table.py --fasta $OUTDIR/best_refs.paired.fasta --meta $COG_PUBLISHED_DIR/majora.latest.metadata.matched.tsv > $OUTDIR/genome_table_$DATESTAMP.csv
python upload_azure.py -c genomics -f $OUTDIR/genome_table_$DATESTAMP.csv

# Make and push variant table
python make_variants_table.py --ref $WUHAN_FP --msa $OUTDIR/naive_msa.fasta > $OUTDIR/variant_table_$DATESTAMP.csv
python upload_azure.py -c genomics -f $OUTDIR/variant_table_$DATESTAMP.csv

# Clean up
rm $OUTDIR/best_refs.paired.fasta
rm $OUTDIR/output.sam
#rm $OUTDIR/naive_msa.fasta
rm $OUTDIR/genome_table_$DATESTAMP.csv
#rm $OUTDIR/variant_table_$DATESTAMP.csv
