"""
Reconstruccion de senales con ruido: pseudoinversa, SVD truncada y
regularizacion de Tikhonov.

Grupo: Carlos Ledo, Matias Desmonteix, Paula Romero
Curso: Algebra Lineal Aplicada (UP 2026-2) - Prof. Oswaldo Velasquez Castanon

Implementa las formulas del material de referencia "Reconstruccion de
senales con ruido" (seccion 2, 4 y 5 del PDF) a partir de la SVD compacta
de A, sin usar funciones de "caja negra" para la parte espectral (se usa
numpy.linalg.svd como la herramienta de clase, pero el ensamblado de la
pseudoinversa, el truncamiento y Tikhonov se hacen a mano con las formulas
deducidas en el documento).
"""

import numpy as np


def svd_compacta(A, tol=None):
    """SVD compacta de A (rango numerico r), devuelve U_r, S_r (vector), V_r.

    tol: tolerancia de rango numerico. Si es None se usa la regla estandar
    de numpy/LAPACK: tol = max(m, n) * eps * sigma_max.
    """
    U, S, Vt = np.linalg.svd(A, full_matrices=False)
    if tol is None:
        m, n = A.shape
        eps = np.finfo(S.dtype).eps
        tol = max(m, n) * eps * S[0]
    r = int(np.sum(S > tol))
    return U[:, :r], S[:r], Vt[:r, :].T, r


def pseudoinversa_solucion(A, b, tol=None):
    """x+ = A+ b = sum_i (u_i^T b / sigma_i) v_i   (Proposicion 2.1 del PDF)."""
    U, S, V, r = svd_compacta(A, tol=tol)
    coef = (U.T @ b) / S
    x = V @ coef
    return x, {"U": U, "S": S, "V": V, "r": r, "coef": coef}


def svd_truncada_solucion(A, b, k, tol=None):
    """x_k = sum_{i=1}^k (u_i^T b / sigma_i) v_i, k <= r (seccion 4)."""
    U, S, V, r = svd_compacta(A, tol=tol)
    k = min(k, r)
    coef = (U[:, :k].T @ b) / S[:k]
    x = V[:, :k] @ coef
    return x, {"U": U, "S": S, "V": V, "r": r, "k": k}


def tikhonov_solucion(A, b, lam, tol=None):
    """x_lambda = sum_i sigma_i/(sigma_i^2+lambda) (u_i^T b) v_i (Teorema 5.1)."""
    U, S, V, r = svd_compacta(A, tol=tol)
    factores = S / (S**2 + lam)
    coef = factores * (U.T @ b)
    x = V @ coef
    f_i = (S**2) / (S**2 + lam)  # factores de filtro f_i(lambda), para diagnostico
    return x, {"U": U, "S": S, "V": V, "r": r, "f": f_i}


def tikhonov_por_minimos_cuadrados_aumentado(A, b, lam):
    """Resuelve el mismo problema por el sistema aumentado [A; sqrt(lam) I] x = [b; 0]
    via QR (np.linalg.lstsq), como contraste con la formula espectral (seccion 9)."""
    n = A.shape[1]
    A_aum = np.vstack([A, np.sqrt(lam) * np.eye(n)])
    b_aum = np.concatenate([b, np.zeros(n)])
    x, *_ = np.linalg.lstsq(A_aum, b_aum, rcond=None)
    return x


def error_relativo(x_hat, x_true):
    return np.linalg.norm(x_hat - x_true) / np.linalg.norm(x_true)


def residuo(A, x_hat, b):
    return np.linalg.norm(A @ x_hat - b)


def construir_senal(n=128, centro=0.3, ancho=80.0, pulso_ini=0.6, pulso_fin=0.75,
                     pulso_altura=0.6):
    """Senal x_j (seccion 8): gaussiana + pulso rectangular, t_j = j/(n-1)."""
    t = np.arange(n) / (n - 1)
    gauss = np.exp(-ancho * (t - centro) ** 2)
    pulso = pulso_altura * ((t >= pulso_ini) & (t <= pulso_fin)).astype(float)
    x = gauss + pulso
    return t, x


def construir_suavizado(t, s=0.03):
    """Matriz de suavizado gaussiano discreto normalizado por filas (seccion 8).
    No simetrica en general: no se asume A = A^T."""
    n = len(t)
    Ti, Tj = np.meshgrid(t, t, indexing="ij")
    A_tilde = np.exp(-((Ti - Tj) ** 2) / (2 * s**2))
    filas = A_tilde.sum(axis=1, keepdims=True)
    A = A_tilde / filas
    return A


def generar_ruido(A, x_true, eta, seed):
    """e = eta * ||A x_true||_2 * z/||z||_2, z ~ N(0,I) con semilla fija."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal(len(x_true))
    Ax = A @ x_true
    e = eta * np.linalg.norm(Ax) * z / np.linalg.norm(z)
    return e


def ejemplo_diagonal():
    """Ejemplo de la seccion 3.1 del PDF, para verificar las formulas."""
    A = np.array([[1.0, 0.0], [0.0, 1e-3]])
    x_true = np.array([1.0, 1.0])
    b = np.array([1.0, 0.011])
    x_exacto = np.linalg.solve(A, b)  # ajuste exacto (residuo 0)
    x_trunc1, _ = svd_truncada_solucion(A, b, k=1)
    x_tik, _ = tikhonov_solucion(A, b, lam=1e-4)
    return {
        "A": A, "x_true": x_true, "b": b,
        "x_exacto": x_exacto,
        "x_trunc1": x_trunc1,
        "x_tik": x_tik,
    }


def regla_discrepancia(A, b, lambdas, delta, tau=1.5):
    """Principio de discrepancia (seccion 7): escoge el mayor lambda de la
    malla cuyo residuo no exceda tau*delta.

    La malla se recorre de mayor a menor (se ordena aqui, sin depender del
    orden en que llegue) y se devuelve el PRIMER lambda que cumple el
    criterio. Como el residuo crece con lambda, ese es el mayor que cumple.

    Correccion (4 nov 2026): la version anterior no se detenia en la primera
    coincidencia y seguia sobrescribiendo 'elegido', por lo que devolvia el
    MENOR lambda que cumple (siempre el extremo inferior de la malla)."""
    for lam in sorted(lambdas, reverse=True):
        x, _ = tikhonov_solucion(A, b, lam)
        res = residuo(A, x, b)
        if res <= tau * delta:
            return lam
    return None
