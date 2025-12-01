# check posizione
# (cpdb_env) Mac:Pitone nome_utente$ pwd
# /Users/nome_utente/Desktop/Pitone

import os

# =========================================================================
# ⚠️ CONFIGURA QUESTO: Usa un percorso generico o la posizione relativa
# =========================================================================
# Sostituisci '/Users/nome_utente/Desktop/Pitone' con un percorso valido per il tuo ambiente
NOME_UTENTE = 'nome_utente_generico'
BASE_PATH = '/PATH_AL_TUO_PROGETTO/Pitone' 
# Se esegui lo script dalla cartella 'Pitone', puoi usare:
# BASE_PATH = os.getcwd() 

# 1. DEFINIZIONE DEI PERCORSI E DELLE VARIABILI (FATTO QUI)
CPDB_DB_PATH = os.path.join(BASE_PATH, 'cellphonedb_v5.0.0.zip') 
META_FILE = os.path.join(BASE_PATH, 'test_meta.txt')
COUNTS_FILE = os.path.join(BASE_PATH, 'test.h5ad') 
DEGS_FILE = os.path.join(BASE_PATH, 'test_degs.txt')
ACTIVE_TFS_FILE = os.path.join(BASE_PATH, 'test_active_tfs.txt')
MICROENVS_FILE = os.path.join(BASE_PATH, 'test_microenviroments.txt')
OUTPUT_DIR = os.path.join(BASE_PATH, 'Risultati_Salvati_Finali') 

import os
import sys

# Importa il metodo di analisi corretto
from cellphonedb.src.core.methods import cpdb_degs_analysis_method
from cellphonedb.utils import db_utils

# ⚠️ DEFINISCE ESATTAMENTE LA CARTELLA DOVE VERRÀ SALVATO IL FILE
# Usa BASE_PATH per la directory di destinazione
TARGET_DIR = BASE_PATH
VERSIONE_DB = 'v5.0.0' 

print(f"Tentativo di scaricamento di CellphoneDB {VERSIONE_DB} in: {TARGET_DIR}")

# Usa la sintassi corretta: target_dir (non cpdb_target_dir)
# Questo creerà il file cellphonedb_v5.0.0.zip nella cartella 'Pitone'
db_utils.download_database(target_dir=TARGET_DIR, cpdb_version=VERSIONE_DB)

print(f"Scaricamento completato. Controlla la cartella {os.path.basename(BASE_PATH)}.")

# ---
# ⚠️ METODO DEGs (Differentially Expressed Genes)
# ---

# 1. CONFIGURAZIONE DEI PERCORSI ASSOLUTI
# BASE_PATH è già definito sopra

# CPDB_DB_PATH è ora impostato su un nome generico, non specifico della versione
CPDB_DB_PATH = os.path.join(BASE_PATH, 'cellphonedb.zip')
META_FILE = os.path.join(BASE_PATH, 'test_meta.txt')
COUNTS_FILE = os.path.join(BASE_PATH, 'test.h5ad')
DEGS_FILE = os.path.join(BASE_PATH, 'test_degs.txt')
ACTIVE_TFS_FILE = os.path.join(BASE_PATH, 'test_active_tfs.txt')
MICROENVS_FILE = os.path.join(BASE_PATH, 'test_microenviroments.txt')

OUTPUT_DIR = os.path.join(BASE_PATH, 'Risultati_Finali_Assoluti')

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("\n--- Avvio dell'analisi CellphoneDB (Metodo DEGs) ---")
print(f"I risultati verranno salvati in: {OUTPUT_DIR}")

# 2. CHIAMATA ALLA FUNZIONE DI ANALISI
try:
    cpdb_results = cpdb_degs_analysis_method.call(
        cpdb_file_path = CPDB_DB_PATH,
        meta_file_path = META_FILE,
        counts_file_path = COUNTS_FILE,
        
        degs_file_path = DEGS_FILE,
        
        # ⚠️ PARAMETRI MODIFICATI PER IL DEBUG ⚠️
        counts_data = 'ensembl', #'gene_name', 
        threshold = 0.00,          
        
        active_tfs_file_path = ACTIVE_TFS_FILE,
        microenvs_file_path = MICROENVS_FILE,
        score_interactions = True,
        
        output_path = OUTPUT_DIR
    )
    print("--- Analisi COMPLETATA CON SUCCESSO! ---")

except Exception as e:
    # Se fallisce qui, significa che il file non esiste fisicamente o è corrotto.
    print(f"ERRORE CRITICO DURANTE L'ANALISI: {e}")

# ---
# ⚠️ METODO STATISTICO (Statistical Analysis)
# ---

# Importa il metodo di analisi statistica (il più rigoroso)
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# -----------------------------------------------------------
# ⚠️ 1. CONFIGURAZIONE DEI PERCORSI ASSOLUTI (MODIFICA QUESTI!)
# -----------------------------------------------------------

# La cartella che contiene TUTTI i tuoi file di input e il database (cellphonedb.zip)
# BASE_PATH è già definito sopra

# USARE I PERCORSI ASSOLUTI PER I TUOI FILE REALI
CPDB_DB_PATH = os.path.join(BASE_PATH, 'cellphonedb.zip')  
META_FILE = os.path.join(BASE_PATH, 'mia_meta.txt')        # ⬅️ IL TUO FILE META
COUNTS_FILE = os.path.join(BASE_PATH, 'miei_counts.h5ad')  # ⬅️ IL TUO FILE DI CONTEGGI 

# Cartella di output 
OUTPUT_DIR = os.path.join(BASE_PATH, 'Risultati_Analisi_Statistica_Reali')


# -----------------------------------------------------------
# ⚠️ 1. CONFIGURAZIONE PER I DATI DI ESEMPIO
# -----------------------------------------------------------

# USA I NOMI DEI FILE DI ESEMPIO CHE ESISTONO
CPDB_DB_PATH = os.path.join(BASE_PATH, 'cellphonedb.zip')
META_FILE = os.path.join(BASE_PATH, 'test_meta.txt')
COUNTS_FILE = os.path.join(BASE_PATH, 'test.h5ad') 
OUTPUT_DIR = os.path.join(BASE_PATH, 'Risultati_Analisi_Statistica_Test')

os.makedirs(OUTPUT_DIR, exist_ok=True)

print("\n--- Avvio dell'analisi STATISTICA CellphoneDB (sui dati di esempio) ---")

# -----------------------------------------------------------
# ⚠️ 2. CHIAMATA ALLA FUNZIONE STATISTICA
# -----------------------------------------------------------

try:
    cpdb_results = cpdb_statistical_analysis_method.call(
        cpdb_file_path = CPDB_DB_PATH,
        meta_file_path = META_FILE,
        counts_file_path = COUNTS_FILE,
        
        # PARAMETRI CHE HANNO FUNZIONATO CON I DATI DI ESEMPIO
        counts_data = 'ensembl',  # Formato che ha funzionato
        threshold = 0.00,         
        
        # Parametri statistici
        iterations = 1000,        # Numero di permutazioni 
        threads = 4,
        score_interactions = True,
        
        output_path = OUTPUT_DIR
    )
    print("--- Analisi STATISTICA COMPLETATA CON SUCCESSO! ---")

except Exception as e:
    print(f"ERRORE CRITICO: {e}", file=sys.stderr)