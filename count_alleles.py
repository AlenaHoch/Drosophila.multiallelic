
import pandas as pd

# change your number of timepoints here
n_timepoints = 15

#import tsv-file (change path here)
# that is the tsv file you originally got from the data_conversion script and also used as input for the F-Matrix script
inp = pd.read_csv("linvilla_1_10_allelefreq.tsv", sep="\t", usecols=["n_alleles"])

# count frequencies
counts = inp["n_alleles"].value_counts()

# make new dataframe
result_df = counts.reset_index()

result_df.columns = ['n_alleles', 'frequency']

# divide by the number of timepoints
result_df['frequency'] = result_df['frequency'] / n_timepoints

result_df.to_csv("n_alleles.tsv", sep="\t", index=False)