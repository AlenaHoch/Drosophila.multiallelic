import pandas as pd

#change the name of your population here:
pop="linvilla"
#change the number of time points you have here:
n_timepoints = 15

inp = pd.read_csv("linvilla_1_10_allelefreq.tsv", sep="\t")

# only every 15th row, since we have every position 15 times in the .tsv-file
inpf = inp.iloc[::n_timepoints]

# only columns chr and pos
df = inpf.iloc[:, 1:3] 

# new column with pos+1
df['end'] = df.iloc[:, 1] + 1

# save as .bed
df.to_csv(pop +'.bed', sep='\t', index=False, header=False)  # no header for .bed