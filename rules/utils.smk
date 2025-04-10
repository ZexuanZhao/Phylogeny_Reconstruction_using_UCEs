def get_asm(wildcards):
    asm = sample_sheet.loc[(wildcards), ["asm"]].dropna()
    return {"asm": asm.asm}