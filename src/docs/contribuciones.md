# Registro de contribuciones

**Nota de transparencia:** el borrador de este avance (código, desarrollo
matemático, README, guion y cronograma) fue redactado con asistencia de IA
(Claude, ver "Uso de IA" en `README.md`) a pedido de Carlos Ledo. Las
casillas de "revisión y verificación" deben completarlas Matías y Paula
cuando revisen su parte asignada, antes de subir el avance — son ellos
quienes dan fe de la exactitud de lo que verificaron, como exige la
consigna ("Todos son responsables de la exactitud del contenido").

| Integrante | Secciones a cargo | Revisión y verificación realizada |
|---|---|---|
| **Carlos Ledo** | Coordinación general del avance; estructura del repo; problema/objetivo/alcance; experimento (construcción de la señal y el ruido) | Pidió y coordinó la generación del avance; pendiente: verificar manualmente la sección de problema/objetivo contra la consigna del docente |
| **Matías Desmonteix** | Desarrollo matemático — pseudoinversa, amplificación del ruido y número de condición (`02_desarrollo_matematico.md`, §§1–2); resultados de comparación de métodos | *(pendiente: revisar las fórmulas de `02_desarrollo_matematico.md` §§1–2 contra el PDF del profesor y las figuras 02–04, y anotar aquí qué verificó)* |
| **Paula Romero** | Desarrollo matemático — SVD truncada y Tikhonov (`02_desarrollo_matematico.md`, §§3–5); resultados de selección del parámetro (curva L y principio de discrepancia) | *(pendiente: revisar las fórmulas de `02_desarrollo_matematico.md` §§3–5 —incluida la monotonía del residuo de §5.2— y la sección corregida del principio de discrepancia en `04_resultados_preliminares.md` §4, y anotar aquí qué verificó)* |

## Reparto para lo que falta (ver también `cronograma.md`)
- Guion del video: borrador conjunto, tiempos a ajustar entre los tres tras
  un ensayo.
- Principio de discrepancia: error de programación corregido el 4 nov
  (Carlos, con asistencia de IA). Falta volver a correr `experimento.py`
  para confirmar los números y comparar contra la curva L.
- Revisión final de redacción y de que el repositorio esté completo antes de
  la entrega: los tres, antes de la fecha límite.
