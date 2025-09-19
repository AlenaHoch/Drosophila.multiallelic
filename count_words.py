import pandas as pd

# load .txt file
with open("genes_" + pop + ".txt", "r") as file:
    words = file.readlines()

# make list out of the gene names
words = [word.strip() for word in words]

# count frequency of the words
df = pd.DataFrame(words, columns=["word"])  # create dataframe with column word
df = df["word"].value_counts().reset_index()  # count frequencies
df.columns = ["genename", "frequency"]  # make header

# calculate the amount of loci found in genes
n_loci_in_genes = df["frequency"].sum()

# load again your loci to count the total amount of loci
inp = pd.read_csv("linvilla_1_10_allelefreq.tsv", sep="\t")
n_loci_total = len(inp)/15

# calculate the proportion of loci found in genes
proportion = (n_loci_in_genes / n_loci_total) * 100

print(f"Total amount of loci: {n_loci_total}")
print(f"Proportion of loci in genes: {proportion:.2f}%\n")

df.to_csv("genes_" + pop + ".tsv", sep='\t', index=False) 