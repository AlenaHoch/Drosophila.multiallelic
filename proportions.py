import gffutils
from Bio import SeqIO

# Input file paths (replace these with your actual file paths)
gff_file = "genes_annotations_612.gff3"
fasta_file = "dmel-all-chromosome-r6.12.fasta"

# Specify chromosomes to analyze (e.g., 2L, 2R, 3L, 3R)
chromosomes_of_interest = ["2L", "2R", "3L", "3R"]

# Create a database from the GFF file for efficient querying
db = gffutils.create_db(gff_file, dbfn="temp.db", force=True, keep_order=True)

# Initialize variables to store aggregated values
total_length_combined = 0
coding_length_combined = 0

# Calculate total sequence length for chromosomes of interest and aggregate it
for record in SeqIO.parse(fasta_file, "fasta"):
    if record.id in chromosomes_of_interest:
        total_length_combined += len(record.seq)

# Calculate gene length for all chromosomes of interest and aggregate it
for feature in db.features_of_type("gene"):
    if feature.chrom in chromosomes_of_interest:
        coding_length_combined += feature.end - feature.start + 1

# Calculate proportion of DNA in genes across all selected chromosomes
if total_length_combined > 0:  # Avoid division by zero
    proportion_coding = (coding_length_combined / total_length_combined) * 100
    
    print(f"   Total Length: {total_length_combined} bp")
    print(f"   Coding Length: {coding_length_combined} bp")
    print(f"   Proportion Coding DNA: {proportion_coding:.2f}%\n")
else:
    print("Total length is zero; cannot calculate proportions.")


