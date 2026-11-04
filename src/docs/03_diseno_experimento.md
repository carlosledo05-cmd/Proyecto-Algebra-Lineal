# Diseño del experimento

Este documento describe cómo se construye y se evalúa el experimento
obligatorio de [1, sec. 9]: la señal, la matriz de observación, el ruido, los
métodos con sus parámetros, las métricas y las verificaciones. Cada decisión
se justifica con el material del profesor o con un resultado de
`02_desarrollo_matematico.md`. El código está en `src/` y la §7 indica qué
función implementa cada paso.

## 1. Señal verdadera

Usamos $n=128$ puntos equiespaciados en $[0,1]$, $t_j=j/(n-1)$ con
$j=0,\dots,127$, y la señal propuesta en [1, sec. 8]:
$$x_j=\exp\!\left[-80\,(t_j-0.3)^2\right]+0.6\cdot\mathbb 1_{[0.6,\,0.75]}(t_j).$$

- **Región suave:** una gaussiana centrada en $0.3$, con desviación
  estándar $1/\sqrt{160}\approx0.079$.
- **Pulso:** altura $0.6$ en $[0.6,\,0.75]$, que contiene 19 puntos de la
  malla ($j=77,\dots,95$).
- $\|x_{true}\|\approx4.96$.

*Por qué esta señal.* Las dos partes ponen a prueba cosas distintas. La
gaussiana es suave y está bien descrita por los primeros vectores singulares,
que corresponden a frecuencias bajas. El pulso tiene dos saltos, cuyo
contenido de alta frecuencia es justamente lo que el suavizado destruye. Así
se ve en una misma reconstrucción dónde funcionan los métodos y dónde chocan
con el límite de lo recuperable. Es además la señal que propone [1].

## 2. Matriz de observación

### 2.1 Construcción

El suavizado gaussiano discreto normalizado por filas es [1, sec. 8]
$$\tilde a_{ij}=\exp\!\left(-\frac{(t_i-t_j)^2}{2s^2}\right),\qquad a_{ij}=\frac{\tilde a_{ij}}{d_i},\qquad d_i=\sum_{\ell=1}^n\tilde a_{i\ell},\qquad s=0.03.$$
En forma matricial, $A=D^{-1}\tilde A$ con $D=\mathrm{diag}(d_1,\dots,d_n)$.
Los bordes usan exactamente esta regla, sin convolución periódica
[1, sec. 8].

*Escala del suavizado.* El paso de la malla es $h=1/127\approx0.0079$, así
que $s\approx3.8\,h$: cada dato observado es un promedio ponderado de unos
diez puntos vecinos a cada lado (hasta $3s$). Usamos $s=0.03$ porque es el
valor que fija el experimento obligatorio [1, sec. 9].

### 2.2 Propiedades

**A conserva las constantes.** Si $\mathbf 1$ es el vector de unos,
$(A\mathbf 1)_i=\sum_j\tilde a_{ij}/d_i=1$, es decir, $A\mathbf 1=\mathbf 1$.
Numéricamente, $\max_i|(A\mathbf 1)_i-1|\approx3\times10^{-16}$.

**$\sigma_1\ge1$.** Como $\|A\mathbf 1\|/\|\mathbf 1\|=1$ y
$\sigma_1=\max_{x\ne0}\|Ax\|/\|x\|$, se cumple $\sigma_1\ge1$. Numéricamente,
$\sigma_1\approx1.001$.

**A no es simétrica.** $\tilde A$ es simétrica, pero
$$A=A^T\iff D^{-1}\tilde A=\tilde AD^{-1}\iff\frac{\tilde a_{ij}}{d_i}=\frac{\tilde a_{ij}}{d_j}\ \ \forall i,j\iff d_i=d_j\ \ \forall i,j,$$
donde la última equivalencia usa que $\tilde a_{ij}>0$. Las sumas de fila no
son iguales:

- En el interior, la suma aproxima una integral gaussiana:
  $d_i\approx\frac1h\int_{-\infty}^{\infty}e^{-u^2/(2s^2)}\,du=\frac{s\sqrt{2\pi}}{h}\approx9.55$.
- En un borde solo queda la mitad de la campana más el término central:
  $d_1\approx(9.55+1)/2\approx5.28$.

Los valores calculados coinciden con estas aproximaciones (9.550 y 5.275),
y $\|A-A^T\|_F\approx0.24$. Por eso, como advierte [1, sec. 8], no se supone
$A=A^T$. En consecuencia, los vectores singulares izquierdos y derechos son
distintos ($u_i\ne v_i$) y hay que usar la SVD, no una descomposición en
valores propios.

### 2.3 Espectro y rango numérico

Para la convolución continua con una gaussiana de desviación $s$, la
frecuencia $\omega$ se atenúa por el factor $e^{-s^2\omega^2/2}$. Esa
atenuación crece muy rápido con la frecuencia, y explica por qué los valores
singulares de $A$ decaen tan rápido (en escala logarítmica, la curva de
`figures/03_espectro_singular.png` se curva hacia abajo). Con la tolerancia
de `02_desarrollo_matematico.md` §1.5:

| Cantidad | Valor |
| --- | --- |
| $\sigma_1$ | $\approx1.001$ |
| Tolerancia de rango $\tau_{\text{rango}}=\max(m,n)\,\varepsilon\,\sigma_1$ | $\approx2.85\times10^{-14}$ |
| Rango numérico $r$ | 89 de 128 (se descartan 39 direcciones) |
| $\sigma_{89}$ y $\sigma_{90}$ | $\approx4.21\times10^{-14}$ y $\approx2.10\times10^{-14}$ |
| Condicionamiento efectivo $\kappa_r=\sigma_1/\sigma_{89}$ | $\approx2.38\times10^{13}$ |
| Parte de $x_{true}$ en las direcciones descartadas | $\approx5.4\,\%$ de $\|x_{true}\|$ |

La última fila es el piso de error de todos los métodos del estudio
(`02_desarrollo_matematico.md` §1.5).

## 3. Ruido

### 3.1 Generación

Se genera $z\in\mathbb R^m$ con entradas normales estándar independientes,
usando el generador `numpy.random.default_rng` con semilla fija
`20261027`, y se escala [1, sec. 8]:
$$e=\eta\,\|Ax_{true}\|\,\frac{z}{\|z\|}.$$
Entonces $\|e\|/\|Ax_{true}\|=\eta\,\|z\|/\|z\|=\eta$ exactamente: el ruido
relativo es el nominal, sin variación aleatoria en su tamaño.

*Distribución de la dirección.* La distribución normal estándar en
$\mathbb R^m$ es invariante por rotaciones, así que $z/\|z\|$ es uniforme en
la esfera unitaria. Este es el hecho que usa `02_desarrollo_matematico.md`
§2.5 para obtener $\mathbb E[(u_i^Te)^2]=\delta^2/m$.

### 3.2 Niveles de ruido

Usamos los tres niveles que fija [1, sec. 9]. Con
$\|Ax_{true}\|\approx4.685$:

| $\eta$ | $\delta=\|e\|$ | $\delta/\sqrt m$ (nivel típico por coeficiente) |
| --- | --- | --- |
| $10^{-4}$ | $4.685\times10^{-4}$ | $4.1\times10^{-5}$ |
| $10^{-3}$ | $4.685\times10^{-3}$ | $4.1\times10^{-4}$ |
| $10^{-2}$ | $4.685\times10^{-2}$ | $4.1\times10^{-3}$ |

### 3.3 Limitación declarada

El código genera el mismo $z$ para los tres niveles, porque usa la misma
semilla en cada llamada. Por eso el ruido de un nivel es exactamente el del
anterior multiplicado por 10, y cada nivel tiene una sola realización. Los
resultados describen esa realización, no un promedio. Repetir el experimento
con varias semillas y reportar media y dispersión queda pendiente.

## 4. Métodos y parámetros

### 4.1 Pseudoinversa numérica

$x^\dagger=V_r\Sigma_r^{-1}U_r^Tb$ con $r=89$, fijado por
$\tau_{\text{rango}}$ (`02_desarrollo_matematico.md` §1.5). No tiene
parámetro libre.

### 4.2 SVD truncada

$x_k$ para $k$ en la lista $\{5,10,20,30,40,60,80,100,128\}$. Los valores
mayores que $r$ se recortan a $r$ y repiten la pseudoinversa numérica. Como
la lista es gruesa para ubicar el óptimo, el siguiente cambio del código es
recorrer todos los $k$ de 1 a $r$ (`02_desarrollo_matematico.md` §3.5).

### 4.3 Tikhonov: malla de λ

Usamos 40 valores logarítmicamente equiespaciados [1, sec. 9]:
$$\lambda_j=\sigma_1^2\cdot10^{-10+10j/39},\qquad j=0,\dots,39,$$
es decir, $\lambda/\sigma_1^2$ entre $10^{-10}$ y $1$, con razón
$q=10^{10/39}\approx1.80$ entre puntos consecutivos.

*Por qué relativa a $\sigma_1^2$.* Si se cambian las unidades de los datos,
$A\to cA$ y $b\to cb$, entonces
$\|cAx-cb\|^2+c^2\lambda\|x\|^2=c^2\left(\|Ax-b\|^2+\lambda\|x\|^2\right)$: la
solución no cambia si $\lambda$ se escala por $c^2$. El parámetro con sentido
es $\lambda/\sigma_1^2$, que no depende de las unidades. Lo comprobamos
numéricamente con $c=7$ (diferencia $\approx10^{-13}$).

*Por qué ese rango.* En el extremo inferior, $\sqrt\lambda=10^{-5}\sigma_1$:
las 48 direcciones más fuertes conservan más del 99 % de su peso y solo se
apagan las de $\sigma_i$ muy pequeño, así que la regularización es débil. En
el superior, $\lambda=\sigma_1^2$ y hasta la dirección más
fuerte queda a la mitad ($f_1=1/2$): regularización excesiva. Entre ambos
extremos quedan todos los valores relevantes. En el experimento, el
$\lambda$ oráculo y el de la discrepancia caen en los índices 16 a 31 de la
malla, lejos de los bordes, así que la malla no recorta el óptimo.

*Por qué 40 puntos.* Con $q\approx1.80$, la regla de discrepancia queda a
menos de un paso del valor exacto (`02_desarrollo_matematico.md` §5.2), y la
curva de error frente a $\lambda$ es suave
(`figures/06_curva_error_oraculo.png`), así que esa resolución basta para
este avance.

### 4.4 Selección del parámetro

- **Oráculo:** el $k$ o el $\lambda$ de la malla con menor error frente a
  $x_{true}$. Solo sirve como referencia.
- **Principio de discrepancia:** el mayor $\lambda_j$ con
  $\|Ax_{\lambda_j}-b\|\le\tau\delta$, con $\tau=1.5$ y $\delta=\|e\|$ conocido
  exactamente, lo que es una idealización del experimento sintético
  (`02_desarrollo_matematico.md` §5.2).
- **Curva L:** graficada; la localización de la esquina queda pendiente.

## 5. Métricas

Para cada reconstrucción $\hat x$ medimos [1, sec. 9]:

| Métrica | Fórmula | ¿Usa $x_{true}$? |
| --- | --- | --- |
| Error relativo | $\|\hat x-x_{true}\|/\|x_{true}\|$ | Sí: solo para evaluar |
| Residuo | $\|A\hat x-b\|$ | No |
| Norma de la solución | $\|\hat x\|$ | No |

El residuo y la norma son lo único que una regla de selección puede usar; el
error relativo es lo que permite juzgarla después.

## 6. Verificaciones

| # | Verificación | Qué comprueba | Resultado esperado | Estado |
| --- | --- | --- | --- | --- |
| 1 | Ejemplo diagonal de [1, sec. 3.1] | Fórmulas de los tres métodos | $(1,11)$, $(1,0)$ y $(0.99990,\,0.10891)$, con errores 10, 1 y 0.8911 | Implementada (se imprime) |
| 2 | Fórmula espectral de Tikhonov frente a sistema aumentado | Implementación de Tikhonov | Diferencia del orden del redondeo | Implementada con `lstsq`; pasar a QR, pendiente |
| 3 | $A\ne A^T$ | Que no se suponga simetría | $\|A-A^T\|_F\approx0.24$ | Implementada |
| 4 | $A\mathbf 1=\mathbf 1$ | Construcción de la matriz | Diferencia $\approx10^{-16}$ | Pendiente |
| 5 | $\|e\|/\|Ax_{true}\|=\eta$ | Escala del ruido | Igualdad hasta el redondeo | Pendiente |
| 6 | Identidad sesgo más ruido de la SVD truncada | Implementación de $x_k$ | Igualdad hasta el redondeo | Pendiente |
| 7 | $\|b_\perp\|<\tau\delta<\|b\|$ | Que el umbral de la discrepancia sea alcanzable | Se cumple para los tres $\eta$ | Pendiente |

Hoy las verificaciones 1 a 3 imprimen el valor obtenido junto al esperado.
Queda pendiente convertir todas en comprobaciones (`assert`) que detengan el
script si fallan.

## 7. Implementación

| Paso | Función (`src/reconstruccion.py`) | Fórmula |
| --- | --- | --- |
| SVD compacta y rango numérico | `svd_compacta` | `02_desarrollo_matematico.md` §1.5 |
| Pseudoinversa numérica | `pseudoinversa_solucion` | §1.3 |
| SVD truncada | `svd_truncada_solucion` | §3.1 |
| Tikhonov (fórmula espectral) | `tikhonov_solucion` | §4.3 |
| Tikhonov (sistema aumentado) | `tikhonov_por_minimos_cuadrados_aumentado` | §4.8 |
| Principio de discrepancia | `regla_discrepancia` | §5.2 |
| Señal, matriz y ruido | `construir_senal`, `construir_suavizado`, `generar_ruido` | §§1–3 de este documento |
| Ejemplo diagonal | `ejemplo_diagonal` | §2.4, §3.4 y §4.5 |

`src/experimento.py` ejecuta las verificaciones, recorre los tres niveles de
ruido y los tres métodos, y guarda las seis figuras en `figures/` y todos los
números en `data/resumen_resultados.json`.

## 8. Cumplimiento de lo que pide el material del profesor

| Requisito de [1] | Estado |
| --- | --- |
| Reproducir el ejemplo diagonal y comprobar fórmulas y residuos (sec. 9.1) | Cumplido |
| Señal y matriz con $s=0.03$; mostrar señal, datos limpios y **datos perturbados** (sec. 9.2) | Falta graficar los datos perturbados |
| $\eta\in\{10^{-4},10^{-3},10^{-2}\}$; registrar semilla y **tolerancias** (sec. 9.3) | Semilla registrada; la tolerancia está en este documento pero no en el JSON |
| Comparar pseudoinversa numérica, SVD truncada y Tikhonov en malla relativa a $\sigma_1^2$ (sec. 9.4) | Cumplido; falta el barrido completo de $k$ |
| Medir error, **residuo y norma**; mostrar espectro y $\lvert u_i^Tb\rvert$ (sec. 9.5) | Faltan residuo y norma en las tablas; los coeficientes se grafican con datos limpios, no con los observados |
| Regla basada en el ruido, comparada con el oráculo (sec. 9.6) | Cumplido (principio de discrepancia) |
| Contraste con el sistema aumentado (sec. 9) | Cumplido; el método es `lstsq`, no QR |
| Repositorio con generador, matriz, semilla, tolerancias, métodos, criterio y tablas (sec. 11) | Cumplido salvo la tolerancia en los resultados |

## Referencias

La numeración [n] remite a la lista única de fuentes, en formato IEEE, del
`README.md`.
