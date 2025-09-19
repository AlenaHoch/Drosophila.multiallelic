import pandas as pd

# TXT-Datei einlesen
with open("genes_" + pop + ".txt", "r") as file:
    words = file.readlines()

# Entferne Leerzeichen/Zeilenumbrüche und erstelle eine Liste der Wörter
words = [word.strip() for word in words]

# Zähle die Häufigkeiten der Wörter mit Pandas
df = pd.DataFrame(words, columns=["word"])  # Erstelle einen DataFrame mit einer Spalte "Word"
df = df["word"].value_counts().reset_index()  # Zähle Vorkommen jedes Wortes und erstelle neuen DataFrame
df.columns = ["word", "frequency"]  # Benenne die Spalten um

df.to_csv("genes_" + pop + ".tsv", sep='\t', index=False) 