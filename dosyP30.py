
import tkinter as tk
from tkinter import ttk, messagebox
import math
import re
import csv
import os

# Constantes
k_B = 1.380649e-23
gamma_H = 2.675e8

VISCOSITES = {
    "H2O": 0.0010,
    "DMSO": 0.0020,
    "CDCl3": 0.0006,
    "MeOD": 0.00055,
    "Acetone-d6": 0.0003,
    "DMF": 0.00092,
    "Acetic": 0.0012,
}

G_MAX_SPEC = {"400": 0.30, "500": 0.35}
DELTA_BOUNDS = (500e-6, 2500e-6)
CSV_FILE = r"U:\CBD_PARTAGE\JMS\DOSY_EXPERIMENTS\experiences_dosy.csv"

# Fonctions calcul
def masse_moleculaire(formule):
    masses = {"H": 1.008, "C": 12.01, "N": 14.01, "O": 16.00, "Cl": 35.45, "Br": 79.90, "I": 126.90}
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formule)
    mw = 0.0
    for elem, count in tokens:
        n = int(count) if count else 1
        mw += masses.get(elem, 0.0) * n
    return mw

def r_h_lineaire(MW, k_lin=3.205e-12):
    return k_lin * (MW ** 0.5)

def r_h_polymere(DP, l_monomere=0.25e-9, nu=0.588):
    return l_monomere * (DP ** nu)

def D_stokes_einstein(r_h, T_K, eta):
    return (k_B * T_K) / (6.0 * math.pi * eta * r_h)

def b_tanner(g, delta, Delta):
    return (gamma_H**2) * (g**2) * (delta**2) * (Delta - delta/3.0)

def attenuation_ratio(D, g_min, g_max, delta, Delta):
    return math.exp(-D * (b_tanner(g_max, delta, Delta) - b_tanner(g_min, delta, Delta)))

def clamp_delta(delta, Delta):
    return min(max(delta, DELTA_BOUNDS[0]), DELTA_BOUNDS[1])

def solve_delta(D, g_max, Delta, ratio_min=0.02, ratio_max=0.05):
    g_min = 0.02 * g_max
    g_max_val = 0.98 * g_max
    ratio_target = (ratio_min + ratio_max) / 2.0
    lo, hi = DELTA_BOUNDS
    for _ in range(60):
        mid = (lo + hi) / 2.0
        R = attenuation_ratio(D, g_min, g_max_val, mid, Delta)
        if abs(R - ratio_target) < 1e-6:
            delta = mid
            break
        if R > ratio_target:
            lo = mid
        else:
            hi = mid
    else:
        delta = (lo + hi) / 2.0
    return clamp_delta(delta, Delta), attenuation_ratio(D, g_min, g_max_val, delta, Delta)

# UI
root = tk.Tk()
root.title("Calcul p30 DOSY RMN")

mode_var = tk.StringVar(value="non_poly")
spectro_var = tk.StringVar()
solvant_var = tk.StringVar()
d20_var = tk.StringVar(value="0.1")
temp_var = tk.StringVar(value="298")
ratio_min_var = tk.StringVar(value="0.02")
ratio_max_var = tk.StringVar(value="0.05")

formule_var = tk.StringVar()
mw_var = tk.StringVar()
formule_mono_var = tk.StringVar()
mw_mono_var = tk.StringVar()
dp_var = tk.StringVar()

result_text = tk.Text(root, height=10, width=60)

last_result = None

# Fonctions UI
def switch_mode():
    for widget in frame_dynamic.winfo_children():
        widget.destroy()
    if mode_var.get() == "non_poly":
        ttk.Label(frame_dynamic, text="Formule brute:").grid(row=0, column=0)
        ttk.Entry(frame_dynamic, textvariable=formule_var, width=15).grid(row=0, column=1)
        ttk.Label(frame_dynamic, text="OU MW:").grid(row=0, column=2)
        ttk.Entry(frame_dynamic, textvariable=mw_var, width=10).grid(row=0, column=3)
    else:
        ttk.Label(frame_dynamic, text="Formule monomère:").grid(row=0, column=0)
        ttk.Entry(frame_dynamic, textvariable=formule_mono_var, width=15).grid(row=0, column=1)
        ttk.Label(frame_dynamic, text="OU MW monomère:").grid(row=0, column=2)
        ttk.Entry(frame_dynamic, textvariable=mw_mono_var, width=10).grid(row=0, column=3)
        ttk.Label(frame_dynamic, text="DP:").grid(row=1, column=0)
        ttk.Entry(frame_dynamic, textvariable=dp_var, width=10).grid(row=1, column=1)

def reset_fields():
    formule_var.set(""); mw_var.set("")
    formule_mono_var.set(""); mw_mono_var.set(""); dp_var.set("")
    spectro_var.set(""); solvant_var.set("")
    d20_var.set("0.1"); temp_var.set("298")
    ratio_min_var.set("0.02"); ratio_max_var.set("0.05")
    result_text.delete("1.0", tk.END)

def calculer():
    try:
        spectro = spectro_var.get().strip()
        solvant = solvant_var.get().strip()
        d20_s = float(d20_var.get()); T_K = float(temp_var.get())
        eta = VISCOSITES.get(solvant); g_spec_max = G_MAX_SPEC.get(spectro)
        if eta is None or g_spec_max is None:
            messagebox.showerror("Erreur", "Spectro ou solvant invalide."); return
        ratio_min = float(ratio_min_var.get()); ratio_max = float(ratio_max_var.get())
        result_text.delete("1.0", tk.END)

        global last_result
        last_result = {"spectro": spectro, "solvant": solvant, "d20": d20_s, "temp": T_K, "DP": "", "MW": 0.0}

        if mode_var.get() == "non_poly":
            MW = 0.0
            if mw_var.get().strip() != "":
                MW = float(mw_var.get())
            elif formule_var.get().strip() != "":
                MW = masse_moleculaire(formule_var.get().strip())
            else:
                messagebox.showerror("Erreur", "Entrer formule ou MW."); return

            r_lin = r_h_lineaire(MW)
            D_lin = D_stokes_einstein(r_lin, T_K, eta)
            delta_lin, _ = solve_delta(D_lin, g_spec_max, d20_s, ratio_min, ratio_max)

            result_text.insert(tk.END, f"MW: {MW:.2f} g/mol\n")
            result_text.insert(tk.END, f"p30 calculé (linéaire): {delta_lin*1e6:.1f} µs\n")

            last_result["MW"] = MW

        else:
            DP = int(dp_var.get()) if dp_var.get().strip() != "" else 0
            if DP <= 0:
                messagebox.showerror("Erreur", "DP invalide."); return
            MW_mono = 0.0
            if mw_mono_var.get().strip() != "":
                MW_mono = float(mw_mono_var.get())
            elif formule_mono_var.get().strip() != "":
                MW_mono = masse_moleculaire(formule_mono_var.get().strip())
            else:
                messagebox.showerror("Erreur", "Entrer formule monomère ou MW."); return

            MW_total = MW_mono * DP
            r_h = r_h_polymere(DP)
            D = D_stokes_einstein(r_h, T_K, eta)
            delta, _ = solve_delta(D, g_spec_max, d20_s, ratio_min, ratio_max)

            result_text.insert(tk.END, f"MW total: {MW_total:.2f} g/mol\n")
            result_text.insert(tk.END, f"p30 calculé (polymère): {delta*1e6:.1f} µs\n")

            last_result["MW"] = MW_mono
            last_result["DP"] = DP

    except Exception as e:
        messagebox.showerror("Erreur", str(e))

def export_csv():
    if not last_result:
        messagebox.showerror("Erreur", "Aucun calcul effectué."); return
    win = tk.Toplevel(root)
    win.title("Entrer p30 optimisé expérimental")
    ttk.Label(win, text="p30 optimisé (µs):").pack(pady=5)
    p30_var = tk.StringVar()
    ttk.Entry(win, textvariable=p30_var).pack(pady=5)
    def save_exp():
        try:
            p30_exp = float(p30_var.get())
            file_exists = os.path.isfile(CSV_FILE)
            with open(CSV_FILE, "a", newline="") as f:
                writer = csv.writer(f)
                if not file_exists:
                    writer.writerow(["spectrometer", "solvent", "molecular_weight", "d20_s", "p30_us", "température_K", "DP"])
                writer.writerow([last_result["spectro"], last_result["solvant"], last_result["MW"], last_result["d20"], p30_exp, last_result["temp"], last_result.get("DP", "")])
            messagebox.showinfo("Export", f"Valeur expérimentale ajoutée à {CSV_FILE}")
            win.destroy()
        except:
            messagebox.showerror("Erreur", "Valeur invalide.")
    ttk.Button(win, text="Ajouter", command=save_exp).pack(pady=10)

# Layout
frm_top = tk.Frame(root); frm_top.pack(pady=5)
ttk.Label(frm_top, text="Mode:").pack(side=tk.LEFT)
ttk.Radiobutton(frm_top, text="Non-Polymère", variable=mode_var, value="non_poly", command=switch_mode).pack(side=tk.LEFT)
ttk.Radiobutton(frm_top, text="Polymère", variable=mode_var, value="poly", command=switch_mode).pack(side=tk.LEFT)

frm_common = tk.Frame(root); frm_common.pack(pady=5)
ttk.Label(frm_common, text="Spectro:").grid(row=0, column=0)
ttk.Combobox(frm_common, textvariable=spectro_var, values=list(G_MAX_SPEC.keys()), width=8).grid(row=0, column=1)
ttk.Label(frm_common, text="Solvant:").grid(row=0, column=2)
ttk.Combobox(frm_common, textvariable=solvant_var, values=list(VISCOSITES.keys()), width=15).grid(row=0, column=3)
ttk.Label(frm_common, text="Δ (s):").grid(row=1, column=0)
ttk.Entry(frm_common, textvariable=d20_var, width=10).grid(row=1, column=1)
ttk.Label(frm_common, text="Température (K):").grid(row=1, column=2)
ttk.Entry(frm_common, textvariable=temp_var, width=10).grid(row=1, column=3)

frame_dynamic = tk.Frame(root); frame_dynamic.pack(pady=5)
switch_mode()

frm_btn = tk.Frame(root); frm_btn.pack(pady=10)
ttk.Button(frm_btn, text="Calculer", command=calculer).pack(side=tk.LEFT, padx=5)
ttk.Button(frm_btn, text="Reset", command=reset_fields).pack(side=tk.LEFT, padx=5)
ttk.Button(frm_btn, text="Exporter calcul", command=export_csv).pack(side=tk.LEFT, padx=5)

result_text.pack(pady=10)

root.mainloop()