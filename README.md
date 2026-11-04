# Reconstrucción de señales con ruido: pseudoinversa, SVD y Tikhonov

**Álgebra Lineal Aplicada — UP 2026-2** · Prof. Oswaldo Velásquez Castañón
**Grupo:** Carlos Ledo Correa, Matías Desmonteix, Paula Romero
**Etiqueta de esta entrega:** avance-parcial

## Estado del avance (entrega parcial)

Esta entrega cubre lo pedido para el avance parcial: problema y objetivo
delimitados, desarrollo matemático completo con justificación, una
implementación funcionando con resultados numéricos y gráficos, y
verificaciones cruzadas contra el ejemplo del material del profesor. **No**
es el proyecto terminado: faltan las ampliaciones (ver `cronograma.md`),
en particular comparar el principio de discrepancia con la curva L y
redactar las conclusiones finales.

**Corrección del 4 de noviembre.** La versión anterior de este avance
reportaba que el principio de discrepancia fallaba. La causa era un error de
programación en `regla_discrepancia` (devolvía el menor $\lambda$ que cumple
el criterio en vez del mayor), ya corregido. Con la corrección, el criterio
queda a 3–13 % del óptimo oráculo. Detalle en
`docs/04_resultados_preliminares.md` §4.

Lo que SÍ está hecho y puede revisarse:
- El problema, los objetivos y el alcance (`docs/01_problema_objetivo.md`).
- El desarrollo matemático completo (`docs/02_desarrollo_matematico.md`),
  siguiendo la notación y las demostraciones del PDF del profesor.
- La implementación en Python de pseudoinversa, SVD truncada y Tikhonov
  (`src/reconstruccion.py`), verificada contra el ejemplo diagonal del PDF.
- El experimento obligatorio completo: señal sintética de 128 puntos,
  suavizado gaussiano, tres niveles de ruido, comparación de los tres
  métodos, curva L y curva de error oráculo (`src/experimento.py`).
- Resultados numéricos (`data/resumen_resultados.json`) y 6 figuras
  (`figures/`), con interpretación inicial (`docs/04_resultados_preliminares.md`).

## Estructura del repositorio

```
├── README.md                       (este archivo)
├── docs/
│   ├── 01_problema_objetivo.md       problema, objetivos, alcance
│   ├── 02_desarrollo_matematico.md   teoría y demostraciones
│   ├── 03_diseno_experimento.md      señal, matriz, ruido, parámetros
│   ├── 04_resultados_preliminares.md tablas, verificaciones, interpretación
│   └── contribuciones.md             registro de aportes
├── src/
│   ├── reconstruccion.py           funciones: pseudoinversa, SVD truncada, Tikhonov
│   └── experimento.py              script que corre todo y genera figuras/datos
├── figures/                        6 figuras PNG generadas por experimento.py
├── data/
│   └── resumen_resultados.json     todos los números reportados, reproducibles
├── guion_video.md                  guion del video final (10 min)
└── cronograma.md                   actividades pendientes y fechas
```

## Cómo ejecutar

Requisitos: Python 3, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib
cd src
python experimento.py
```

Esto imprime en consola la verificación del ejemplo diagonal y las métricas
por nivel de ruido, y regenera las 6 figuras en `../figures/` y el resumen
en `../data/resumen_resultados.json`. El experimento usa una semilla fija
(`20261027`) para ser reproducible exactamente.

Para solo revisar (sin ejecutar nada), basta con leer en este orden:
1. `docs/01_problema_objetivo.md` (qué se estudia y para qué)
2. `docs/02_desarrollo_matematico.md` (teoría)
3. `docs/03_diseno_experimento.md` (cómo se construye el experimento)
4. `src/reconstruccion.py` (implementación de las fórmulas)
5. `docs/04_resultados_preliminares.md` + `figures/` (qué salió y qué significa)

## Fuentes (formato IEEE)

[1] O. Velásquez Castañón, "Reconstrucción de señales con ruido: pseudoinversa,
valores singulares y regularización de Tikhonov," Material de referencia
para proyectos, Universidad del Pacífico, Álgebra Lineal Aplicada 2026-2,
oct. 2026.

[2] P. C. Hansen, *Rank-Deficient and Discrete Ill-Posed Problems: Numerical
Aspects of Linear Inversion*. Philadelphia, PA, USA: SIAM, 1998.
doi:10.1137/1.9780898719697.

[3] P. C. Hansen, *Discrete Inverse Problems: Insight and Algorithms*.
Philadelphia, PA, USA: SIAM, 2010. doi:10.1137/1.9780898718836.

[4] NumPy Developers, "numpy.linalg.matrix_rank," NumPy Reference. [En
línea]. Disponible:
https://numpy.org/doc/stable/reference/generated/numpy.linalg.matrix_rank.html.
Consultado: 4 nov. 2026.

## Uso de IA

- **Herramienta y versión:** Claude (Anthropic), modelo Claude Sonnet 5, usado
  dentro del entorno de trabajo del integrante Carlos Ledo. Las correcciones
  del 4 de noviembre (error del principio de discrepancia y actualización de
  la documentación) se hicieron con Claude Opus 5.5, también a pedido de
  Carlos Ledo.
- **Finalidad:** redactar el desarrollo matemático a partir del PDF del
  profesor, escribir el código de `reconstruccion.py` y `experimento.py`,
  generar las figuras, y organizar el repositorio, el guion de video y el
  cronograma.
- **Instrucciones relevantes dadas:** reproducir fielmente las fórmulas del
  material del profesor; verificar todo contra el ejemplo numérico del PDF
  antes de aplicarlo al experimento nuevo; documentar honestamente cualquier
  resultado inesperado en vez de ocultarlo (así se detectó y corrigió el
  error del principio de discrepancia).
- **Partes asistidas:** desarrollo matemático (docs/), implementación
  numérica (src/), generación de figuras, README, guion de video y
  cronograma.
- **Verificaciones numéricas (incluidas en el código):** el ejemplo diagonal
  se contrasta número por número contra el PDF del profesor (coincide
  exactamente); la solución de Tikhonov se calcula por dos caminos
  independientes (fórmula espectral vs. sistema aumentado resuelto por QR),
  con diferencia del orden de $10^{-12}$. Los resultados corregidos del
  principio de discrepancia se obtuvieron con una reimplementación
  independiente del experimento, que reproduce exactamente los demás
  números de `data/resumen_resultados.json`.
- **Revisión manual por el grupo:** la revisión de que las fórmulas de
  `docs/02_desarrollo_matematico.md` correspondan a las del PDF original está
  asignada por secciones y se registra en `docs/contribuciones.md`, donde
  cada integrante anota qué verificó. Cada integrante es responsable de la
  exactitud del contenido que declaró haber revisado.
