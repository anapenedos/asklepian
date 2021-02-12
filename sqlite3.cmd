# NOTE We use COG_ID not COG-ID because why on earth is it a hyphen
CREATE TABLE variants("COG_ID" TEXT, "Position" INTEGER, "Reference_Base" TEXT, "Alternate_Base" TEXT, "Is_Indel" INTEGER);
.mode csv
.import /cephfs/covid/bham/results/variants/latest/naive_variant_table.csv variants
CREATE INDEX idx_variants_position ON variants(Position);
CREATE INDEX idx_variants_cog_id ON variants("COG_ID");

# Remove the header because it is faster than tailing the file and hitting disk again
DELETE FROM variants WHERE cog_id="COG-ID";

# Output counts
#SELECT COUNT(*) from variants;
#SELECT COUNT(DISTINCT "COG-ID") from variants;
