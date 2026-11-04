# Guion del video final (≈10 minutos)

Tema: Reconstrucción de señales con ruido — pseudoinversa, SVD y Tikhonov.
Reparto de tiempos según la estructura sugerida por la pauta de la entrega
(1′ problema · 3′ explicación espectral/Tikhonov · 4′ reconstrucciones · 2′
conclusiones). **La asignación de quién habla en cada bloque es una
propuesta inicial — ajustarla según de quién sea más cómoda cada parte antes
de grabar el video final.**

| # | Bloque | Tiempo | Contenido | A cargo (propuesto) |
|---|---|---|---|---|
| 1 | Problema | 0:00–1:00 | Qué es un problema inverso ($b=Ax_{true}+e$), por qué suavizar + ruido hace que un buen ajuste no sea una buena reconstrucción (gancho: el ejemplo diagonal, residuo 0 pero error 10) | Carlos Ledo |
| 2 | Teoría — pseudoinversa y amplificación | 1:00–2:30 | SVD compacta, fórmula $x^\dagger=\sum u_i^Tb/\sigma_i\, v_i$, por qué $\sigma_i$ pequeño amplifica el ruido, número de condición | Matías Desmonteix |
| 3 | Teoría — SVD truncada y Tikhonov | 2:30–4:00 | Truncamiento (corte abrupto) vs. Tikhonov (filtro suave $f_i(\lambda)$), unicidad de $x_\lambda$, trade-off sesgo/varianza | Paula Romero |
| 4 | Experimento y datos | 4:00–5:30 | Construcción de la señal (gaussiana + pulso), matriz de suavizado ($s=0.03$, no simétrica), niveles de ruido $\eta$, verificación del ejemplo diagonal contra el PDF | Carlos Ledo |
| 5 | Resultados: comparación de métodos | 5:30–7:30 | Mostrar figuras `02`–`04`: espectro singular, reconstrucciones con Tikhonov oráculo para los tres $\eta$, tabla de errores relativos (pseudoinversa inservible, SVD truncada ≈ Tikhonov, Tikhonov ligeramente mejor) | Matías Desmonteix |
| 6 | Resultados: selección del parámetro | 7:30–8:30 | Curva L (`05`) y curva de error oráculo (`06`); principio de discrepancia: sin usar $x_{true}$ elige un $\lambda$ cuyo error queda entre 3 % y 13 % sobre el óptimo oráculo (tabla de `docs/04_resultados_preliminares.md` §4). Contar brevemente el error de programación que encontramos y corregimos (la búsqueda devolvía el menor $\lambda$ en vez del mayor), como ejemplo de por qué verificamos los resultados | Paula Romero |
| 7 | Conclusiones y pendientes | 8:30–10:00 | Qué se aprendió (residuo bajo ≠ buena reconstrucción; regularización necesaria cuando $\kappa_2$ es grande), qué falta para la entrega final (posible ampliación a imágenes o a la penalización por diferencias), agradecimientos | Carlos Ledo, Matías Desmonteix y Paula Romero (cierre conjunto) |

## Notas de producción
- Grabar los bloques 4–6 compartiendo pantalla con las figuras de
  `figures/` ya generadas (no es necesario correr el código en vivo).
- Bloque 1 y 7 pueden grabarse a cámara, sin pantalla compartida.
- Pendiente antes de grabar la versión final: cronometrar un ensayo
  completo y ajustar los tiempos si algún bloque se alarga.
