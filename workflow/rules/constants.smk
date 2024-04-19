# Versions
YAK_VERSION = "0.1"
KRAKEN2_VERSION="2.1.3"
SEQTK_VERSION="1.4"

FILTER_Z=config.get("z_filter", -2)
KRAKEN2_DB=config["kraken2_db"]

# Config parameters
ILLUMINA_ENDEDNESS=config.get("illumina","pair")
