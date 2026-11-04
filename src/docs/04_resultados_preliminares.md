# Resultados preliminares e interpretación inicial

Todos los números provienen de `src/experimento.py` (semilla fija
`20261027`) y están guardados íntegros en `data/resumen_resultados.json`.
Las figuras están en `figures/`.

## 1. Verificación: ejemplo diagonal

Reprodujimos el ejemplo exacto del PDF del profesor (`A=diag(1,10⁻³)`) como
prueba de que nuestra implementación de pseudoinversa, SVD truncada y
Tikhonov es correcta antes de aplicarla a un problema nuevo:

| Método | $\hat x$ | error absoluto | residuo |
|---|---|---|---|
| Ajuste exacto | $(1,\ 11)$ | $10$ | $0$ |
| SVD truncada $k=1$ | $(1,\ 0)$ | $1$ | $0.011$ |
| Tikhonov $\lambda=10^{-4}$ | $(0.9999,\ 0.1089)$ | $0.8911$ | $0.0109$ |

Coincide exactamente con el PDF (figura `01_ejemplo_diagonal.png`).

## 2. La matriz de observación

$A$ (suavizado gaussiano, $n=128$, $s=0.03$) resultó **no simétrica**
(verificado numéricamente) y severamente mal condicionada:

- Rango numérico: $89$ de $128$.
- $\sigma_1\approx1.001$, $\sigma_{89}\approx4.2\times10^{-14}$.
- $\kappa_2(A)\approx2.4\times10^{13}$.

Esto confirma que la pseudoinversa sin regularizar es inviable: dividir
entre valores singulares del orden de $10^{-14}$ amplifica el ruido por un
factor equivalente.

## 3. Comparación de métodos por nivel de ruido

| $\eta$ | pseudoinversa (error rel.) | SVD truncada óptima (k, error rel.) | Tikhonov óptimo ($\lambda/\sigma_1^2$, error rel.) |
|---|---|---|---|
| $10^{-4}$ | $8.3\times10^{7}$ | $k=40$, $0.093$ | $1.3\times10^{-6}$, $0.093$ |
| $10^{-3}$ | $8.3\times10^{8}$ | $k=40$, $0.104$ | $1.3\times10^{-5}$, $0.096$ |
| $10^{-2}$ | $8.3\times10^{9}$ | $k=30$, $0.123$ | $4.6\times10^{-4}$, $0.111$ |

("óptimo" = oráculo, es decir el que minimiza el error contra $x_{true}$
conocido; no es un criterio aplicable sin la verdad, ver
`02_desarrollo_matematico.md`, §5.)

**Interpretación:**
- La pseudoinversa sin regularizar es inútil (error relativo $>10^7$ en
  todos los casos): confirma en un problema de tamaño real lo que el
  ejemplo diagonal muestra en miniatura.
- SVD truncada y Tikhonov dan errores del mismo orden (~9–12%), con
  Tikhonov consistentemente un poco mejor que el mejor truncamiento: el
  filtro suave aprovecha parcialmente las direcciones intermedias que el
  corte abrupto descarta por completo.
- El $k$ óptimo (30–40 de 89 direcciones con $\sigma_i>0$ numérico) baja al
  subir el ruido, como predice la teoría: con más ruido conviene confiar en
  menos direcciones.
- Verificamos la fórmula espectral de Tikhonov contra un sistema aumentado
  resuelto por mínimos cuadrados (QR): la diferencia entre ambos métodos es
  del orden de $10^{-12}$–$10^{-13}$ en norma, es decir, coinciden hasta el
  error de redondeo (columna "verificación" en el JSON).

## 4. Principio de discrepancia (corregido)

Con $\tau=1.5$ y la misma malla de 40 puntos logarítmicos
($\lambda/\sigma_1^2$ entre $10^{-10}$ y $1$), el principio de discrepancia
(`02_desarrollo_matematico.md`, §5) elige un $\lambda$ cercano al óptimo oráculo
sin usar $x_{true}$:

| $\eta$ | Tikhonov oráculo ($\lambda/\sigma_1^2$, error rel.) | Discrepancia ($\lambda/\sigma_1^2$, error rel.) | razón de errores |
|---|---|---|---|
| $10^{-4}$ | $1.3\times10^{-6}$, $0.093$ | $1.3\times10^{-5}$, $0.095$ | $1.03$ |
| $10^{-3}$ | $1.3\times10^{-5}$, $0.096$ | $4.6\times10^{-4}$, $0.109$ | $1.13$ |
| $10^{-2}$ | $4.6\times10^{-4}$, $0.111$ | $8.9\times10^{-3}$, $0.122$ | $1.10$ |

**Interpretación:** el error queda entre un 3 % y un 13 % por encima del
óptimo oráculo. El $\lambda$ elegido es sistemáticamente mayor que el
oráculo (sobre-regulariza un poco), lo esperable con $\tau>1$: el criterio
acepta un residuo algo mayor que el nivel del ruido. Además, el $\lambda$
elegido crece con $\eta$, como debe ser.

**Corrección respecto de la versión anterior de este avance.** En la
versión subida el 2 de noviembre, el principio de discrepancia parecía
fallar: elegía siempre el $\lambda$ más pequeño de la malla
($\lambda/\sigma_1^2\approx10^{-10}$), con errores relativos de $0.71$,
$7.1$ y $71$, es decir, entre $7.7\times$ y $640\times$ el óptimo oráculo.
Atribuimos ese fallo al mal condicionamiento de $A$ y a la resolución de la
malla, pero la causa real era un error de programación en
`regla_discrepancia`: el bucle recorría la malla de mayor a menor y no se
detenía en la primera coincidencia, de modo que devolvía el **menor**
$\lambda$ que cumple el criterio en vez del mayor. La pista era que el
$\lambda$ elegido era idéntico para los tres niveles de ruido. Corregido el
bucle (y sin cambiar la malla ni $\tau$), el criterio funciona como describe
la teoría.

## 5. Verificaciones aplicadas

- Fórmulas de pseudoinversa, SVD truncada y Tikhonov reproducen
  exactamente el ejemplo numérico del propio material del curso.
- $A\ne A^T$ confirmado numéricamente (no se asumió, se verificó).
- Solución espectral de Tikhonov vs. sistema aumentado (QR): coinciden
  hasta error de redondeo, para los tres niveles de ruido.
- Error relativo calculado solo donde se conoce $x_{true}$ (síntesis
  controlada); residuo y norma de la solución, que sí se pueden calcular sin
  conocerla, se reportan por separado (curva L, figura `05_curva_L.png`).
