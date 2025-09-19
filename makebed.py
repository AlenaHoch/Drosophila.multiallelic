import pandas as pd

#change the name of your population here:
pop="linvilla"

inp = pd.read_csv("linvilla_1_10_allelefreq.tsv", sep="\t")

# Behalte nur jede fünfzehnte Zeile
inpf = inp.iloc[::15]

# Nur die ersten beiden Spalten behalten
df = inpf.iloc[:, 1:3]  # Auswahl der ersten zwei Spalten mit iloc

# Neue dritte Spalte hinzufügen: Wert der zweiten Spalte (Index 1) + 1
df['C'] = df.iloc[:, 1] + 1

# Als BED-Datei speichern (tabulatorgetrennt)
df.to_csv(pop +'.bed', sep='\t', index=False, header=False)  # Kein Header für BED-Format