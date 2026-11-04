# Actividades pendientes y cronograma

La fecha de entrega final "será comunicada por el docente en Blackboard"
(aún no confirmada al 4 de noviembre de 2026). El plan abajo cuenta desde el
2 de noviembre y debe ajustarse en cuanto se conozca la fecha real.

## Pendientes identificados en este avance

1. ~~**Arreglar el principio de discrepancia.**~~ **Resuelto (4 nov).** El
   fallo no se debía a la malla sino a un error en `regla_discrepancia`
   (devolvía el menor $\lambda$ que cumple en vez del mayor). Corregido, el
   criterio queda a 3–13 % del óptimo oráculo (ver
   `docs/04_resultados_preliminares.md` §4). Falta volver a correr
   `experimento.py` para confirmar los números. Para la entrega final queda
   compararlo contra la esquina de la curva L.
2. **Revisión cruzada de fórmulas.** Matías y Paula deben verificar,
   sección por sección, que `docs/02_desarrollo_matematico.md` corresponda
   exactamente al PDF del profesor (ver `docs/contribuciones.md`).
3. **Guion del video — ensayo cronometrado.** Grabar un ensayo y ajustar los
   tiempos de `guion_video.md` si algún bloque se alarga.
4. **Conclusiones finales.** Redactar la sección de conclusiones del
   proyecto completo (no solo del avance).
5. **Posible ampliación.** Si alcanza el tiempo: extender el experimento a
   una imagen pequeña (sección 10 del PDF) o agregar la penalización por
   diferencias ($\lambda\|Lx\|_2^2$).
6. **Entrega formal.** Crear la etiqueta `avance-parcial` en el commit de
   este estado, y que el representante entregue nombres, tema, título y
   enlace al repo en la actividad de Blackboard.
7. **Retroalimentación de la presentación inicial.** Describir qué nos
   observaron y qué cambiamos por ello (`docs/01_problema_objetivo.md` §7).
8. **Completar el experimento de la sección 9 del PDF.** Graficar los datos
   perturbados, graficar $|u_i^Tb|$ con los datos observados, agregar el
   residuo y la norma de la solución a las tablas y registrar la tolerancia
   de rango en el JSON (`docs/03_diseno_experimento.md` §8).
9. **Mejoras de método anotadas en el desarrollo.** Recorrer todos los $k$ de
   1 a $r$, resolver el sistema aumentado por QR, analizar la sensibilidad a
   $\tau$, localizar la esquina de la curva L y repetir con varias semillas
   (`docs/02_desarrollo_matematico.md` §§3.5, 4.8 y 5.2–5.3;
   `docs/03_diseno_experimento.md` §3.3).
10. **Verificaciones automáticas.** Convertir las comprobaciones en `assert`
    y agregar las pendientes ($A\mathbf 1=\mathbf 1$, escala del ruido,
    identidad sesgo más ruido, alcanzabilidad del umbral;
    `docs/03_diseno_experimento.md` §6).

## Cronograma propuesto (desde el 2 de noviembre)

| Fecha | Actividad | Responsable(s) |
|---|---|---|
| 4 nov (atrasada; prevista el 2–3 nov) | Revisión cruzada del desarrollo matemático (pendiente 2) | Matías, Paula |
| 4 nov | Redactar la retroalimentación de la presentación inicial (pendiente 7) | Los tres |
| 4 nov | Corregir el error del principio de discrepancia (pendiente 1) — hecho; falta volver a correr `experimento.py` | Carlos |
| 4 nov | Crear etiqueta `avance-parcial`, confirmar que el repo corre desde cero en una máquina limpia | Carlos |
| 4 nov | Entrega del avance parcial en Blackboard (nombres, tema, título, enlace) | Representante del grupo |
| Después del avance | Ensayo cronometrado del guion de video (pendiente 3) | Los tres |
| Después del avance | Figuras y tablas faltantes de la sección 9 y nueva ejecución de `experimento.py` (pendiente 8) | Carlos (propuesto) |
| Después del avance | Barrido completo de $k$ y sistema aumentado por QR (pendiente 9) | Matías (propuesto) |
| Después del avance | Sensibilidad a $\tau$, esquina de la curva L y varias semillas (pendiente 9) | Paula (propuesto) |
| Después del avance | Verificaciones con `assert` (pendiente 10) | Carlos (propuesto) |
| Después del avance | Redactar conclusiones finales (pendiente 4) | Los tres |
| Antes de la fecha final (a confirmar) | Ampliación opcional (pendiente 5), si el tiempo lo permite | Los tres |
| Antes de la fecha final (a confirmar) | Grabación y edición del video final | Los tres |
| Antes de la fecha final (a confirmar) | Revisión final del repositorio y entrega | Los tres |

*Actualizar la columna de fechas de las últimas tres filas en cuanto el
docente comunique la fecha de entrega final en Blackboard.*
