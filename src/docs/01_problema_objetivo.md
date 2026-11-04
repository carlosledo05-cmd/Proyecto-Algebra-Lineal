# Problema, objetivo y alcance

Este documento delimita el problema, los objetivos y el alcance del avance.
El desarrollo matemático está en `02_desarrollo_matematico.md`, el diseño del
experimento en `03_diseno_experimento.md` y los resultados en
`04_resultados_preliminares.md`. La notación es la de
`02_desarrollo_matematico.md` y sigue la del material del profesor [1].

## 1. Problema

Una señal desconocida $x_{true}\in\mathbb R^n$ se observa a través de un
proceso lineal que la suaviza y la contamina con ruido:
$$b = Ax_{true} + e.$$
El suavizado gaussiano atenúa casi por completo las componentes de alta
frecuencia. En la base de vectores singulares, la información de $x_{true}$ en
la dirección $v_i$ llega a los datos multiplicada por $\sigma_i$, y muchos
$\sigma_i$ son menores que el ruido. Por eso un vector que reproduce $b$ con
gran precisión puede estar muy lejos de $x_{true}$ [1, sec. 1]. El problema es
reconstruir información útil de $x_{true}$ sin amplificar excesivamente lo que
no sabemos, que es $e$. La teoría de problemas inversos discretos [2], [3]
ofrece el marco para este análisis.

## 2. Naturaleza del estudio

Es un estudio numérico comparativo con datos sintéticos controlados.
Construimos $x_{true}$, $A$ y $e$, de modo que el error de cada reconstrucción
puede medirse exactamente.

- **Lo que variamos:** el método (pseudoinversa numérica, SVD truncada,
  Tikhonov), su parámetro ($k$ o $\lambda$) y el nivel de ruido $\eta$.
- **Lo que mantenemos fijo:** el tamaño $n$, el ancho del suavizado $s$, la
  señal y la semilla del ruido.
- **Lo que medimos:** el error relativo $\|\hat x-x_{true}\|/\|x_{true}\|$, el
  residuo $\|A\hat x-b\|$ y la norma $\|\hat x\|$ de cada reconstrucción
  $\hat x$.

Conocer $x_{true}$ obliga a una distinción que mantenemos en todo el trabajo:
el error frente a $x_{true}$ sirve para evaluar un método, nunca para elegir
su parámetro [1, sec. 7].

## 3. Preguntas

1. ¿Por qué la solución de mínimos cuadrados de mínima norma, que ajusta $b$
   casi exactamente, puede estar lejos de $x_{true}$?
2. ¿Cuánto reducen ese error la SVD truncada y la regularización de Tikhonov,
   y qué pierden a cambio?
3. ¿Puede elegirse el parámetro sin conocer $x_{true}$, y cuánto se pierde
   frente al mejor parámetro posible?

## 4. Objetivos

**Objetivo general.** Comparar la pseudoinversa numérica, la SVD truncada y
la regularización de Tikhonov en la reconstrucción de una señal suavizada con
ruido, cuantificando su error frente a la señal verdadera y evaluando si una
regla basada solo en el nivel de ruido elige un parámetro cercano al óptimo.

**Objetivos específicos y cómo se verifican:**

1. Deducir las soluciones de los tres métodos a partir de la SVD y
   comprobarlas con el ejemplo diagonal de [1, sec. 3.1]. Se cumple si los
   valores obtenidos coinciden con los publicados.
2. Implementar los tres métodos y contrastar la fórmula espectral de Tikhonov
   con el sistema aumentado resuelto por un solucionador estable
   [1, sec. 9]. Se cumple si ambas soluciones difieren solo en el orden del
   error de redondeo.
3. Para cada $\eta\in\{10^{-4},10^{-3},10^{-2}\}$, medir el error relativo, el
   residuo y la norma de cada método, y determinar el $k$ y el $\lambda$ que
   minimizan el error (selección oráculo). Se cumple al reportar esas
   cantidades en tablas.
4. Aplicar el principio de discrepancia sin consultar $x_{true}$ y reportar,
   para cada $\eta$, la razón entre su error y el error oráculo.

## 5. Alcance

- Señal unidimensional de $n=128$ puntos, suavizado gaussiano con $s=0.03$ y
  tres niveles de ruido. Usamos una sola realización del ruido por nivel
  (semilla fija), así que los resultados describen esa realización y no un
  promedio.
- Métodos: pseudoinversa numérica, SVD truncada para varios $k$ y Tikhonov
  para 40 valores de $\lambda/\sigma_1^2$ entre $10^{-10}$ y $1$
  [1, sec. 9].
- Selección del parámetro: principio de discrepancia (implementado) y curva L
  (graficada; la localización de su esquina queda pendiente).
- Quedan fuera, como posibles ampliaciones, la reconstrucción de imágenes, la
  penalización por diferencias y la validación cruzada generalizada
  [1, sec. 10].
- No buscamos recuperar lo que la observación destruye: la parte de
  $x_{true}$ en las direcciones que $A$ anula no puede recuperarse con ningún
  método (`02_desarrollo_matematico.md`, §1.3 y §1.5).

## 6. Conexión con el curso

El trabajo usa mínimos cuadrados y ecuaciones normales, los subespacios
fundamentales de una matriz y las proyecciones ortogonales, la SVD, el número
de condición, las matrices simétricas definidas positivas y las normas
vectoriales y matriciales. La SVD se usa como herramienta estudiada en clase,
sin reimplementarla [1, Resumen].

## 7. Cambios respecto de la presentación inicial

> **[Completar por el grupo]** Indicar en dos o tres líneas qué
> retroalimentación recibimos en la presentación de organización y qué
> cambiamos por ella en el problema, el objetivo o el alcance. La consigna de
> la entrega parcial pide incorporarla.

## Referencias

La numeración [n] remite a la lista única de fuentes, en formato IEEE, del
`README.md`.
