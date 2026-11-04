"""
Script principal del experimento obligatorio (seccion 9 del PDF).

Ejecuta:
 1. Verificacion del ejemplo diagonal (seccion 3.1 / 4 / 5).
 2. Construccion de la senal suavizada (n=128, s=0.03).
 3. Para eta en {1e-4, 1e-3, 1e-2}: compara pseudoinversa, SVD truncada
    (varios k) y Tikhonov (malla log de lambda/sigma1^2).
 4. Guarda figuras en ../figures y un resumen numerico en ../data.

Correr desde la carpeta src/:
    python experimento.py
Requiere: numpy, matplotlib (pip install numpy matplotlib).
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from reconstruccion import (
    svd_compacta, pseudoinversa_solucion, svd_truncada_solucion,
    tikhonov_solucion, tikhonov_por_minimos_cuadrados_aumentado,
    error_relativo, residuo, construir_senal, construir_suavizado,
    generar_ruido, ejemplo_diagonal, regla_discrepancia,
)

FIG_DIR = "../figures"
DATA_DIR = "../data"
SEMILLA = 20261027  # fecha en que se fijo el material de referencia (reproducible)

resumen = {}

# ---------------------------------------------------------------------
# 1. Ejemplo diagonal (verificacion de formulas contra el PDF)
# ---------------------------------------------------------------------
ej = ejemplo_diagonal()
err_exacto = error_relativo(ej["x_exacto"], ej["x_true"])
err_trunc1 = error_relativo(ej["x_trunc1"], ej["x_true"])
err_tik = error_relativo(ej["x_tik"], ej["x_true"])
res_exacto = residuo(ej["A"], ej["x_exacto"], ej["b"])
res_trunc1 = residuo(ej["A"], ej["x_trunc1"], ej["b"])
res_tik = residuo(ej["A"], ej["x_tik"], ej["b"])

resumen["ejemplo_diagonal"] = {
    "x_exacto": ej["x_exacto"].tolist(),
    "x_trunc1": ej["x_trunc1"].tolist(),
    "x_tik_lam1e-4": ej["x_tik"].tolist(),
    "error_abs_exacto": float(np.linalg.norm(ej["x_exacto"] - ej["x_true"])),
    "error_abs_trunc1": float(np.linalg.norm(ej["x_trunc1"] - ej["x_true"])),
    "error_abs_tik": float(np.linalg.norm(ej["x_tik"] - ej["x_true"])),
    "residuo_exacto": float(res_exacto),
    "residuo_trunc1": float(res_trunc1),
    "residuo_tik": float(res_tik),
}
print("=== Ejemplo diagonal ===")
print(f"x_exacto = {ej['x_exacto']}  (esperado (1, 11))")
print(f"error absoluto exacto = {np.linalg.norm(ej['x_exacto']-ej['x_true']):.4f} (esperado 10)")
print(f"x_trunc(k=1) = {ej['x_trunc1']}  (esperado (1, 0))")
print(f"residuo_trunc1 = {res_trunc1:.4f} (esperado 0.011)")
print(f"x_tikhonov(lam=1e-4) = {ej['x_tik']}  (esperado approx (0.99990, 0.10891))")
print(f"error_abs_tik = {np.linalg.norm(ej['x_tik']-ej['x_true']):.4f} (esperado approx 0.8911)\n")

# Figura 1 del PDF (barras Verdad / Ajuste exacto / Tikhonov)
fig, ax = plt.subplots(figsize=(6, 4))
etiquetas = ["Componente 1", "Componente 2"]
ancho = 0.25
x_pos = np.arange(2)
ax.bar(x_pos - ancho, ej["x_true"], ancho, label="Verdad")
ax.bar(x_pos, ej["x_exacto"], ancho, label="Ajuste exacto")
ax.bar(x_pos + ancho, ej["x_tik"], ancho, label="Tikhonov")
ax.set_xticks(x_pos)
ax.set_xticklabels(etiquetas)
ax.set_ylabel("Valor reconstruido")
ax.set_title("Reproduccion: residuo cero no garantiza buena reconstruccion")
ax.legend()
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/01_ejemplo_diagonal.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------
# 2. Senal suavizada (seccion 8)
# ---------------------------------------------------------------------
n = 128
s = 0.03
t, x_true = construir_senal(n=n)
A = construir_suavizado(t, s=s)
simetrica = np.allclose(A, A.T)
print(f"=== Senal suavizada n={n}, s={s} ===")
print(f"A es simetrica: {simetrica} (se espera False)")

datos_limpios = A @ x_true

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(t, x_true, label="Senal verdadera $x_{true}$")
ax.plot(t, datos_limpios, label="Datos limpios $Ax_{true}$")
ax.set_xlabel("t")
ax.set_ylabel("amplitud")
ax.set_title(f"Senal con region suave y pulso (n={n}, s={s})")
ax.legend()
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/02_senal_y_datos_limpios.png", dpi=150)
plt.close(fig)

# espectro singular (comun a todos los niveles de ruido, depende solo de A)
U, S, V, r = svd_compacta(A)
coef_b_limpio = np.abs(U.T @ datos_limpios)

fig, axs = plt.subplots(1, 2, figsize=(10, 4))
axs[0].semilogy(np.arange(1, len(S) + 1), S, "o-", ms=3)
axs[0].set_xlabel("i")
axs[0].set_ylabel("$\\sigma_i$")
axs[0].set_title("Espectro singular de A")
axs[1].semilogy(np.arange(1, len(S) + 1), coef_b_limpio, "o-", ms=3, color="C1")
axs[1].set_xlabel("i")
axs[1].set_ylabel("$|u_i^T b|$ (datos limpios)")
axs[1].set_title("Coeficientes espectrales de los datos")
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/03_espectro_singular.png", dpi=150)
plt.close(fig)

resumen["senal"] = {
    "n": n, "s": s, "A_simetrica": bool(simetrica), "rango_numerico": int(r),
    "sigma_max": float(S[0]), "sigma_min": float(S[-1]),
    "kappa2_A": float(S[0] / S[-1]),
}

# ---------------------------------------------------------------------
# 3. Comparacion de metodos para distintos niveles de ruido eta
# ---------------------------------------------------------------------
etas = [1e-4, 1e-3, 1e-2]
k_lista = [5, 10, 20, 30, 40, 60, 80, 100, 128]
lam_rel = np.logspace(-10, 0, 40)  # lambda / sigma_1^2
lambdas = lam_rel * S[0] ** 2

resumen["experimentos"] = {}

fig_sig, ax_sig = plt.subplots(figsize=(7, 4.5))
ax_sig.plot(t, x_true, "k--", lw=1.5, label="Verdad")

fig_lcurve, ax_lcurve = plt.subplots(figsize=(6, 5))

colores = {1e-4: "C0", 1e-3: "C1", 1e-2: "C2"}

for eta in etas:
    e = generar_ruido(A, x_true, eta, seed=SEMILLA)
    b = A @ x_true + e
    delta = np.linalg.norm(e)

    # pseudoinversa (sin regularizar)
    x_pinv, _ = pseudoinversa_solucion(A, b)
    err_pinv = error_relativo(x_pinv, x_true)
    res_pinv = residuo(A, x_pinv, b)

    # SVD truncada: barrer k
    errores_trunc = []
    residuos_trunc = []
    for k in k_lista:
        xk, _ = svd_truncada_solucion(A, b, k)
        errores_trunc.append(error_relativo(xk, x_true))
        residuos_trunc.append(residuo(A, xk, b))
    k_oraculo = k_lista[int(np.argmin(errores_trunc))]
    err_trunc_oraculo = min(errores_trunc)

    # Tikhonov: barrer malla de lambda
    errores_tik = []
    residuos_tik = []
    normas_tik = []
    for lam in lambdas:
        xl, _ = tikhonov_solucion(A, b, lam)
        errores_tik.append(error_relativo(xl, x_true))
        residuos_tik.append(residuo(A, xl, b))
        normas_tik.append(np.linalg.norm(xl))
    idx_oraculo = int(np.argmin(errores_tik))
    lam_oraculo = lambdas[idx_oraculo]
    err_tik_oraculo = errores_tik[idx_oraculo]

    # regla de discrepancia (sin usar x_true), tau=1.5
    lambdas_desc = lambdas[::-1]  # de mayor a menor
    lam_discrepancia = regla_discrepancia(A, b, lambdas_desc, delta, tau=1.5)
    if lam_discrepancia is not None:
        x_desc, _ = tikhonov_solucion(A, b, lam_discrepancia)
        err_discrepancia = error_relativo(x_desc, x_true)
    else:
        err_discrepancia = None

    # contraste: formula espectral vs sistema aumentado (QR), en lambda oraculo
    x_qr = tikhonov_por_minimos_cuadrados_aumentado(A, b, lam_oraculo)
    x_spec, _ = tikhonov_solucion(A, b, lam_oraculo)
    diff_metodos = float(np.linalg.norm(x_qr - x_spec))

    resumen["experimentos"][str(eta)] = {
        "delta_norma_ruido": float(delta),
        "pseudoinversa": {"error_relativo": float(err_pinv), "residuo": float(res_pinv)},
        "svd_truncada": {
            "k_oraculo": int(k_oraculo),
            "error_oraculo": float(err_trunc_oraculo),
            "k_vs_error": list(zip([int(k) for k in k_lista], [float(e) for e in errores_trunc])),
        },
        "tikhonov": {
            "lambda_oraculo": float(lam_oraculo),
            "lambda_rel_oraculo": float(lam_oraculo / S[0] ** 2),
            "error_oraculo": float(err_tik_oraculo),
            "lambda_discrepancia": None if lam_discrepancia is None else float(lam_discrepancia),
            "error_discrepancia": err_discrepancia,
        },
        "verificacion_qr_vs_formula_espectral": diff_metodos,
    }

    print(f"\n=== eta = {eta} (||e||={delta:.4g}) ===")
    print(f"pseudoinversa: error_rel={err_pinv:.4g}, residuo={res_pinv:.4g}")
    print(f"SVD truncada oraculo: k={k_oraculo}, error_rel={err_trunc_oraculo:.4g}")
    print(f"Tikhonov oraculo: lambda/sigma1^2={lam_oraculo/S[0]**2:.3g}, error_rel={err_tik_oraculo:.4g}")
    if lam_discrepancia is not None:
        print(f"Tikhonov discrepancia (tau=1.5): lambda/sigma1^2={lam_discrepancia/S[0]**2:.3g}, error_rel={err_discrepancia:.4g}")
    else:
        print("Tikhonov discrepancia: no se encontro lambda en la malla que cumpla el criterio")
    print(f"|formula espectral - sistema aumentado QR| = {diff_metodos:.2e} (verificacion)")

    # curva de reconstruccion para este eta (con el lambda oraculo)
    x_mejor, _ = tikhonov_solucion(A, b, lam_oraculo)
    ax_sig.plot(t, x_mejor, color=colores[eta], alpha=0.85,
                label=f"Tikhonov oraculo, $\\eta$={eta:g}")

    ax_lcurve.loglog(residuos_tik, normas_tik, "-", color=colores[eta],
                      label=f"$\\eta$={eta:g}")
    ax_lcurve.loglog([res_pinv], [np.linalg.norm(x_pinv)], "x", color=colores[eta])

ax_sig.set_xlabel("t")
ax_sig.set_ylabel("amplitud")
ax_sig.set_title("Reconstruccion con Tikhonov (lambda oraculo) para distintos niveles de ruido")
ax_sig.legend(fontsize=8)
fig_sig.tight_layout()
fig_sig.savefig(f"{FIG_DIR}/04_reconstrucciones_tikhonov.png", dpi=150)
plt.close(fig_sig)

ax_lcurve.set_xlabel("$\\|Ax_\\lambda-b\\|$ (residuo)")
ax_lcurve.set_ylabel("$\\|x_\\lambda\\|$ (norma de la solucion)")
ax_lcurve.set_title("Curva L para distintos niveles de ruido")
ax_lcurve.legend()
fig_lcurve.tight_layout()
fig_lcurve.savefig(f"{FIG_DIR}/05_curva_L.png", dpi=150)
plt.close(fig_lcurve)

# ---------------------------------------------------------------------
# Figura adicional: error vs lambda (curva de error, caso eta=1e-2) con
# el minimo oraculo marcado, para mostrar el filtro espectral en accion.
# ---------------------------------------------------------------------
eta_demo = 1e-2
e_demo = generar_ruido(A, x_true, eta_demo, seed=SEMILLA)
b_demo = A @ x_true + e_demo
errores_demo = []
for lam in lambdas:
    xl, _ = tikhonov_solucion(A, b_demo, lam)
    errores_demo.append(error_relativo(xl, x_true))
idx_min = int(np.argmin(errores_demo))

fig, ax = plt.subplots(figsize=(6.5, 4.5))
ax.semilogx(lam_rel, errores_demo, "-o", ms=3)
ax.axvline(lam_rel[idx_min], color="red", ls="--",
           label=f"minimo oraculo ($\\lambda/\\sigma_1^2$={lam_rel[idx_min]:.2e})")
ax.set_xlabel("$\\lambda/\\sigma_1^2$")
ax.set_ylabel("error relativo $\\|x_\\lambda-x_{true}\\|/\\|x_{true}\\|$")
ax.set_title(f"Curva de error oraculo vs. $\\lambda$ ($\\eta$={eta_demo:g})")
ax.legend()
fig.tight_layout()
fig.savefig(f"{FIG_DIR}/06_curva_error_oraculo.png", dpi=150)
plt.close(fig)

# ---------------------------------------------------------------------
# Guardar resumen numerico
# ---------------------------------------------------------------------
with open(f"{DATA_DIR}/resumen_resultados.json", "w") as f:
    json.dump(resumen, f, indent=2)

print("\nFiguras guardadas en figures/, resumen numerico en data/resumen_resultados.json")
