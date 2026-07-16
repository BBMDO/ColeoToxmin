import pandas as pd
from pathlib import Path

STRUCT_FILE = Path("02_features/struct/candidates_struct_features.parquet")

if not STRUCT_FILE.exists():
    raise SystemExit(f"[ERRO] Arquivo não encontrado: {STRUCT_FILE}")

df = pd.read_parquet(STRUCT_FILE)

print("\n==============================")
print("VALIDAÇÃO 03h STRUCT FEATURES")
print("==============================\n")

print(f"N candidatos: {len(df)}\n")

required_cols = [
    "dssp_n_res",
    "sec_helix_frac",
    "sec_sheet_frac",
    "sec_coil_frac",
]

missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise SystemExit(f"[ERRO] Colunas DSSP ausentes: {missing}")

print("✔ Colunas DSSP presentes")

# ---- Teste 1: tudo zerado?
zero_frac = (df["dssp_n_res"].fillna(0) == 0).mean()
print(f"Frac. dssp_n_res == 0: {zero_frac:.2%}")

if zero_frac > 0.8:
    print("🚨 ALERTA: >80% dos candidatos sem DSSP válido")
else:
    print("✔ DSSP parece ter sido calculado")

# ---- Teste 2: estatísticas básicas
print("\nEstatísticas dssp_n_res:")
print(df["dssp_n_res"].describe())

# ---- Teste 3: frações de secondary structure
print("\nMédias das frações:")
print(df[["sec_helix_frac","sec_sheet_frac","sec_coil_frac"]].mean())

# ---- Teste 4: soma das frações (~1?)
df["sec_sum"] = (
    df["sec_helix_frac"].fillna(0)
    + df["sec_sheet_frac"].fillna(0)
    + df["sec_coil_frac"].fillna(0)
)

print("\nResumo soma das frações:")
print(df["sec_sum"].describe())

bad_sum = ((df["sec_sum"] > 1.2) | (df["sec_sum"] < 0.8)).mean()
print(f"Frac. com soma fora de ~1 (0.8–1.2): {bad_sum:.2%}")

print("\n==============================\n")
