# Check posizione
# (cpdb_env) Mac:Pitone nome_utente$ pwd
# /Users/nome_utente/Desktop/Pitone

# in bash => conda activate cpdb_env 
# cd /Users/nome_pc/Desktop/Piton
# ls ### osserva se compare il tuo script in '.py'
# avvia: (cpdb_env) Mac:Piton nome_pc$ python Git4students2.py

# nel py deve essere presente
import pandas as pd
import os
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# 1. Definisci i percorsi (possono stare fuori dal blocco main)
cpdb_file_path = 'cellphonedb.zip'
meta_file_path = 'cellphonedb_meta_2026_2.txt'
counts_file_path = 'cellphonedb_counts_2026_final.txt'
out_path = 'results/analisi_2026_2'

# 2. Protezione per il multiprocessing (FONDAMENTALE! EVITI CHE IL CODICE VADA AVANTI IN LOOP!)
if __name__ == '__main__':
    
    # Crea la cartella di output se non esiste
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # Debug dei Barcode
    print("Verifica coerenza file...")
    counts_test = pd.read_csv(counts_file_path, sep='\t', index_col=0)
    meta_test = pd.read_csv(meta_file_path, sep='\t', header=None)

    print(f"Esempio Barcode Counts: {counts_test.columns[0]}")
    print(f"Esempio Barcode Meta: {meta_test.iloc[0, 0]}")

    missing = set(meta_test[0]) - set(counts_test.columns)
    if missing:
        print(f"ERRORE: Ci sono {len(missing)} cellule nel Meta che non sono nei Counts!")
        print(f"Esempio di cellula mancante: {list(missing)[0]}")
    else:
        print("Allineamento Barcode: OK!")

    # Lancio ufficiale di CellphoneDB
    print("Inizio analisi statistica CellphoneDB (potrebbe richiedere tempo)...")
    
    cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path = cpdb_file_path,
        meta_file_path = meta_file_path,
        counts_file_path = counts_file_path,
        counts_data = 'hgnc_symbol',
        threshold = 0.1,
        iterations = 1000,
        threads = 4, # Specifica i thread qui
        output_path = out_path
    )

    # Risultati finali
    print("Analisi completata!")
    print("Primi risultati dei p-values:")
    print(cpdb_results['pvalues'].head())
# controlla la barra di caricamento 
# Fonte: https://github.com/ventolab/CellphoneDB/tree/master

    # Per i plots 
    # Python: in bash => pip install ktplotspy
    # R => library(ktplots)
  
