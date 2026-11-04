# Desarrollo matemático

**Proyecto:** Reconstrucción de señales con ruido: pseudoinversa, valores singulares y regularización de Tikhonov
**Curso:** Álgebra Lineal Aplicada (UP 2026-2), Prof. Oswaldo Velásquez Castañón
**Grupo:** Carlos Ledo, Matías Desmonteix, Paula Romero

Este documento sigue la formulación y la notación del material de referencia
del profesor [1]. Las demostraciones se presentan completas, incluidos los
pasos intermedios que [1] deja implícitos, y cada resultado se conecta con la
decisión del experimento que justifica. La redacción se preparó con
asistencia de IA; la revisión de cada sección por parte del grupo se registra
en `docs/contribuciones.md`.
El problema, los objetivos y el alcance están en `01_problema_objetivo.md`; el
modelo de señal, la matriz y los parámetros del experimento, en
`03_diseno_experimento.md`.

## Notación

- Los vectores son columnas. $\|\cdot\|$ es la norma euclídea para vectores y
  la norma espectral (inducida por la euclídea) para matrices.
- $A\in\mathbb R^{m\times n}$ es la matriz de observación; en el experimento
  $m=n=128$.
- SVD completa: $A=U\Sigma V^T$, con $U\in\mathbb R^{m\times m}$ y
  $V\in\mathbb R^{n\times n}$ ortogonales, columnas $u_i$ y $v_i$, y valores
  singulares $\sigma_1\ge\sigma_2\ge\cdots\ge\sigma_{\min(m,n)}\ge0$.
- $r$ es el rango (numérico, ver §1.5). $U_r$, $V_r$ son las primeras $r$
  columnas de $U$, $V$, y $\Sigma_r=\mathrm{diag}(\sigma_1,\dots,\sigma_r)$.
- $\operatorname{Im}A$ y $\ker A$ son la imagen y el núcleo de $A$; $P_S$ es
  la proyección ortogonal sobre el subespacio $S$.
- $x_{true}$ es la señal verdadera, $e$ el ruido, $\delta=\|e\|$ y
  $\eta=\|e\|/\|Ax_{true}\|$ el ruido relativo.
- $\lambda>0$ es directamente el coeficiente de la penalización de Tikhonov,
  no $\lambda^2$ [1, sec. 5].
- $\varepsilon=2^{-52}\approx2.22\times10^{-16}$ es la precisión de máquina en
  doble precisión.

## 1. Mínimos cuadrados y la pseudoinversa numérica

### 1.1 Existencia del mínimo y ecuaciones normales

El ajuste directo minimiza
$$f(x)=\|Ax-b\|^2=x^TA^TAx-2b^TAx+b^Tb.$$
Su gradiente y su hessiana son
$$\nabla f(x)=2\left(A^TAx-A^Tb\right),\qquad \nabla^2 f(x)=2A^TA.$$
Para todo $z$, $z^TA^TAz=\|Az\|^2\ge0$, así que $A^TA$ es semidefinida
positiva y $f$ es convexa. En una función convexa diferenciable todo punto
estacionario es un mínimo global, de modo que los minimizadores son
exactamente las soluciones de las **ecuaciones normales**
$$A^TAx=A^Tb.$$

*Estas ecuaciones siempre tienen solución.* Primero, $\ker(A^TA)=\ker A$: si
$A^TAz=0$, entonces $\|Az\|^2=z^TA^TAz=0$, luego $Az=0$; la otra inclusión es
inmediata. Como $A^TA$ es simétrica, su imagen es el complemento ortogonal de
su núcleo:
$$\operatorname{Im}(A^TA)=\ker(A^TA)^\perp=(\ker A)^\perp=\operatorname{Im}A^T.$$
Por tanto $A^Tb\in\operatorname{Im}(A^TA)$ y existe al menos una solución.

*Unicidad.* La solución es única si y solo si $\ker(A^TA)=\{0\}$, es decir,
si y solo si $\ker A=\{0\}$ (columnas linealmente independientes)
[1, sec. 2]. En otro caso, el conjunto de minimizadores es $x^*+\ker A$ para
cualquier minimizador $x^*$.

*Interpretación geométrica.* Las ecuaciones normales dicen que
$A^T(b-Ax^*)=0$: el residuo óptimo es ortogonal a $\operatorname{Im}A$, así
que $Ax^*$ es la proyección ortogonal de $b$ sobre $\operatorname{Im}A$.

### 1.2 La SVD y los subespacios fundamentales

De $A=U\Sigma V^T$ se obtiene $Av_i=\sigma_iu_i$ y $A^Tu_i=\sigma_iv_i$. Si
$r$ es el número de valores singulares positivos:

- $\{u_1,\dots,u_r\}$ es una base ortonormal de $\operatorname{Im}A$. En
  efecto, $Ax=\sum_{i=1}^r\sigma_i(v_i^Tx)\,u_i$ pertenece a su espacio
  generado, y cada $u_i=A(v_i/\sigma_i)$ está en la imagen.
- $\{v_{r+1},\dots,v_n\}$ es una base ortonormal de $\ker A$. En efecto,
  $Av_j=0$ para $j>r$, y si $Ax=0$ entonces $\sigma_i(v_i^Tx)=0$, luego
  $v_i^Tx=0$ para $i\le r$.
- En consecuencia, $\{v_1,\dots,v_r\}$ es una base ortonormal de
  $(\ker A)^\perp=\operatorname{Im}A^T$.

La SVD compacta conserva solo la parte con $\sigma_i>0$:
$$A=U_r\Sigma_rV_r^T=\sum_{i=1}^r\sigma_iu_iv_i^T,\qquad\sigma_1\ge\cdots\ge\sigma_r>0.$$

### 1.3 Solución de mínima norma

**Proposición 2.1** [1, sec. 2]. La solución de mínimos cuadrados de menor
norma es
$$x^\dagger=A^\dagger b=\sum_{i=1}^r\frac{u_i^Tb}{\sigma_i}v_i,\qquad A^\dagger=V_r\Sigma_r^{-1}U_r^T.$$

*Demostración.* Todo $x\in\mathbb R^n$ se escribe de forma única como
$x=V_rz+x_0$, con $z=V_r^Tx\in\mathbb R^r$ y $x_0=P_{\ker A}\,x$. Entonces
$Ax=U_r\Sigma_rz$. Descomponemos también los datos:
$$b=U_rU_r^Tb+b_\perp,\qquad b_\perp=(I-U_rU_r^T)\,b\perp\operatorname{Im}A.$$
Así, $Ax-b=U_r(\Sigma_rz-U_r^Tb)-b_\perp$, donde el primer término está en
$\operatorname{Im}A$ y el segundo es ortogonal a ella. Por el teorema de
Pitágoras, y porque $U_r$ tiene columnas ortonormales,
$$\|Ax-b\|^2=\sum_{i=1}^r\left(\sigma_iz_i-u_i^Tb\right)^2+\|b_\perp\|^2.$$
Cada sumando depende de una sola coordenada $z_i$ y se anula en
$z_i=u_i^Tb/\sigma_i$; el residuo mínimo es $\|b_\perp\|$, que ninguna $x$
puede reducir. Este valor no depende de $x_0$, así que todos los vectores
$V_rz+x_0$ con esas $z_i$ son minimizadores. Entre ellos,
$\|x\|^2=\|z\|^2+\|x_0\|^2$ porque $V_rz\perp x_0$, y la norma es mínima con
$x_0=0$. $\blacksquare$

*Qué significa para nuestro problema.* Como $A^\dagger A=V_rV_r^T$ es la
proyección sobre $(\ker A)^\perp$, al aplicar $A^\dagger$ a
$b=Ax_{true}+e$ queda
$$x^\dagger-x_{true}=-P_{\ker A}\,x_{true}+A^\dagger e.$$
El primer término es la información que la observación destruyó y que ningún
método recupera; el segundo es el ruido amplificado, que se estudia en §2. La
mínima norma es una elección dentro del conjunto de ajustes, no una garantía
de recuperar la señal verdadera [1, sec. 2].

### 1.4 Por qué no resolvemos las ecuaciones normales

De la SVD, $A^TA=V_r\Sigma_r^2V_r^T$: los valores singulares de $A^TA$ son
$\sigma_i^2$. Para rango completo esto da
$\kappa_2(A^TA)=(\sigma_1/\sigma_n)^2=\kappa_2(A)^2$ [1, sec. 2]. El error de
redondeo al resolver un sistema crece con el número de condición de la matriz
que se usa. En nuestro experimento $\sigma_1/\sigma_r\approx2.4\times10^{13}$
(§1.5), cuyo cuadrado, $5.7\times10^{26}$, supera con creces
$1/\varepsilon\approx4.5\times10^{15}$. Dicho de otro modo: al formar $A^TA$ en
punto flotante, toda dirección con $\sigma_i<\sqrt{\varepsilon}\,\sigma_1
\approx1.5\times10^{-8}$ queda por debajo de la precisión de máquina y se
pierde. Por eso trabajamos con la SVD de $A$ y, como contraste, con el sistema
aumentado de §4; nunca formamos $A^TA$ ni la inversa $(A^TA+\lambda I)^{-1}$.

### 1.5 Rango numérico, tolerancia y pseudoinversa numérica

En aritmética exacta, $r$ es el número de valores singulares positivos. En
punto flotante, la SVD calculada es la SVD exacta de una matriz $A+E$ con
$\|E\|$ del orden de $\varepsilon\|A\|$, y por la desigualdad de Weyl para
valores singulares, $|\sigma_i(A+E)-\sigma_i(A)|\le\|E\|$. Por tanto, un
$\sigma_i$ calculado del orden de $\varepsilon\sigma_1$ no se distingue de
cero. Usamos la regla por defecto de NumPy [4]:
$$\tau_{\text{rango}}=\max(m,n)\,\varepsilon\,\sigma_1,\qquad r=\#\{\,i:\sigma_i>\tau_{\text{rango}}\,\}.$$

Para nuestra matriz, $\sigma_1\approx1.001$ y
$\tau_{\text{rango}}\approx2.85\times10^{-14}$. Como
$\sigma_{89}\approx4.21\times10^{-14}>\tau_{\text{rango}}>\sigma_{90}\approx2.10\times10^{-14}$,
el rango numérico es $r=89$ y se descartan 39 de las 128 direcciones.

Esto tiene tres consecuencias:

1. **Lo que llamamos pseudoinversa es una pseudoinversa numérica.** Calcular
   $x^\dagger$ con $r=89$ equivale a una SVD truncada con $k=89$, con el corte
   fijado por la precisión de máquina y no por el ruido. Siguiendo
   [1, sec. 4], la llamamos *pseudoinversa numérica* e informamos su umbral,
   en lugar de presentarla como una solución "sin regularizar". Como ese corte
   no tiene relación con el nivel de ruido, sigue amplificándolo enormemente
   (§2).
2. **Número de condición.** [1, sec. 3] define $\kappa_2(A)=\sigma_1/\sigma_n$.
   Para nuestra $A$, $\sigma_1/\sigma_{128}\approx9\times10^{18}$, pero
   $\sigma_{128}\approx1\times10^{-19}$ es puro redondeo y ese cociente no
   tiene significado. Reportamos en su lugar el condicionamiento efectivo de
   la parte retenida, $\kappa_r=\sigma_1/\sigma_r\approx2.38\times10^{13}$,
   que es la cota de amplificación del ruido de la pseudoinversa numérica
   (§2).
3. **Piso de error común a todos los métodos.** Toda reconstrucción de este
   trabajo pertenece a $\operatorname{span}\{v_1,\dots,v_{89}\}$. Si $P$ es la
   proyección sobre $\operatorname{span}\{v_{90},\dots,v_{128}\}$, entonces
   cualquier reconstrucción $\hat x$ cumple
   $\|\hat x-x_{true}\|\ge\|Px_{true}\|$. Para nuestra señal,
   $\|Px_{true}\|/\|x_{true}\|\approx0.054$: ningún método del estudio puede
   bajar del 5.4 % de error relativo.

## 2. Por qué se amplifica el ruido

### 2.1 Coeficientes espectrales de los datos

Como $u_i^TA=\sigma_iv_i^T$, cada coeficiente de los datos en la base
$\{u_i\}$ mezcla señal y ruido de forma distinta:
$$u_i^Tb=\sigma_i\,(v_i^Tx_{true})+u_i^Te.$$
La señal llega atenuada por $\sigma_i$; el ruido, no. Al dividir entre
$\sigma_i$, como hace la pseudoinversa numérica, la señal se recupera
exactamente pero el ruido queda multiplicado por $1/\sigma_i$:
$$\frac{u_i^Tb}{\sigma_i}=v_i^Tx_{true}+\frac{u_i^Te}{\sigma_i},\qquad A^\dagger e=\sum_{i=1}^r\frac{u_i^Te}{\sigma_i}v_i\quad[1,\ \text{sec. }3].$$
Si $\sigma_i$ es pequeño, una perturbación pequeña en la dirección $u_i$
produce una oscilación grande en la dirección $v_i$ de la solución.

### 2.2 Cota absoluta

Como los $v_i$ son ortonormales,
$$\|A^\dagger e\|^2=\sum_{i=1}^r\frac{(u_i^Te)^2}{\sigma_i^2}\le\frac{1}{\sigma_r^2}\sum_{i=1}^r(u_i^Te)^2=\frac{\|U_r^Te\|^2}{\sigma_r^2}\le\frac{\|e\|^2}{\sigma_r^2},$$
donde la última desigualdad usa que proyectar sobre
$\operatorname{Im}A$ no aumenta la norma. Por tanto
$$\|A^\dagger e\|\le\sigma_r^{-1}\|e\|.$$
[1, sec. 3] enuncia esta cota para rango completo; la demostración vale igual
para cualquier $r$, con $\sigma_r$ el menor valor singular retenido. La
igualdad se alcanza cuando $e$ es múltiplo de $u_r$: el peor caso es un ruido
alineado con la dirección retenida más débil.

### 2.3 Cota relativa y número de condición

Como $\|A\|=\sigma_1$, se cumple $\|Ax_{true}\|\le\sigma_1\|x_{true}\|$.
Dividiendo la cota anterior entre esta desigualdad,
$$\frac{\|A^\dagger e\|}{\|x_{true}\|}\le\frac{\|e\|/\sigma_r}{\|Ax_{true}\|/\sigma_1}=\kappa_r\,\eta,\qquad\kappa_r=\frac{\sigma_1}{\sigma_r}.$$
La contribución del ruido al error relativo puede llegar a $\kappa_r$ veces el
ruido relativo. Con $\kappa_r\approx2.38\times10^{13}$ (§1.5), un ruido de
$\eta=10^{-4}$ admite un error relativo de hasta $2.4\times10^{9}$. Los
errores observados para la pseudoinversa numérica (§3 de
`04_resultados_preliminares.md`) equivalen a un factor de amplificación de
$8.3\times10^{11}$, por debajo de la cota: el ruido no está alineado con
$u_r$, sino repartido entre todas las direcciones.

El error es lineal en $e$. Como el experimento usa el mismo vector $z$ para
los tres niveles de ruido (ver `03_diseno_experimento.md`), el error de la pseudoinversa numérica crece
exactamente en un factor 10 de un nivel al siguiente. Esto explica la
regularidad de la tabla de resultados, y también indica que los tres niveles
no son realizaciones independientes del ruido.

### 2.4 Ejemplo diagonal (verificación)

Sea [1, sec. 3.1]
$$A=\begin{pmatrix}1&0\\0&10^{-3}\end{pmatrix},\qquad x_{true}=\begin{pmatrix}1\\1\end{pmatrix},\qquad b=\begin{pmatrix}1\\0.011\end{pmatrix}.$$
Su SVD es inmediata: $U=V=I$, $\sigma_1=1$ y $\sigma_2=10^{-3}$. Los datos
limpios son $Ax_{true}=(1,\,0.001)^T$, así que el ruido es $e=(0,\,0.01)^T$.

- **Solución exacta:** $x^\dagger=(1/1,\ 0.011/10^{-3})^T=(1,\,11)^T$, con
  residuo $\|Ax^\dagger-b\|=0$.
- **Ruido amplificado:** $A^\dagger e=(0,\ 0.01/10^{-3})^T=(0,\,10)^T$, y el
  error es $\|x^\dagger-x_{true}\|=10$.
- **Dirección que más amplifica:** la segunda, $v_2=(0,1)^T$, porque tiene el
  menor valor singular; multiplica el ruido por $1/\sigma_2=1000$.
- **Comprobación de la cota relativa:** $\kappa=1000$ y
  $\eta=0.01/\|(1,\,0.001)\|\approx0.01$, luego $\kappa\eta\approx10$. El error
  relativo observado es $10/\|(1,1)\|\approx7.07\le10$.

El ejemplo responde dos preguntas de defensa de [1, sec. 11]: una solución
puede tener residuo cero y ser mala porque ajusta también el ruido, y la
dirección que más amplifica es la de menor valor singular. La implementación
reproduce exactamente estos valores (`04_resultados_preliminares.md`, §1).

### 2.5 Condición discreta de Picard

Sin ruido, $v_i^Tx_{true}=u_i^T(Ax_{true})/\sigma_i$. Para que esos cocientes
se mantengan acotados, los coeficientes $|u_i^T(Ax_{true})|$ deben decrecer al
menos tan rápido como $\sigma_i$: es la condición discreta de Picard [3].

Con ruido, el término $u_i^Te$ no decrece con $i$. En nuestro experimento,
$e=\delta\,z/\|z\|$ tiene dirección uniformemente distribuida, y como $U$ es
ortogonal, también la tiene $U^Te$; por eso
$\mathbb E\big[(u_i^Te)^2\big]=\delta^2/m$ para todo $i$. Mientras
$|u_i^T(Ax_{true})|\gg\delta/\sqrt m$, el coeficiente $u_i^Tb$ es confiable;
cuando cae por debajo de ese nivel, $u_i^Tb$ es casi solo ruido, y dividirlo
entre $\sigma_i$ amplifica ese ruido. El índice donde los coeficientes de la
señal cruzan el nivel $\delta/\sqrt m$ anticipa hasta qué dirección conviene
confiar, es decir, el $k$ óptimo de §3. Por eso [1, sec. 9] pide mostrar el
espectro singular junto con los coeficientes $|u_i^Tb|$.

## 3. SVD truncada

### 3.1 Definición

Para $1\le k\le r$ [1, sec. 4],
$$x_k=\sum_{i=1}^k\frac{u_i^Tb}{\sigma_i}v_i=V_k\Sigma_k^{-1}U_k^Tb.$$
Equivalentemente, $x_k=A_k^\dagger b$ con
$A_k=\sum_{i=1}^k\sigma_iu_iv_i^T$. Por el teorema de Eckart y Young, $A_k$ es
la mejor aproximación de rango $k$ de $A$ en norma espectral, con
$\|A-A_k\|=\sigma_{k+1}$. La SVD truncada reemplaza $A$ por esa aproximación y
resuelve el nuevo problema por mínima norma.

### 3.2 Descomposición exacta del error: ruido propagado y sesgo

Por §2.1, $x_k=\sum_{i\le k}(v_i^Tx_{true})\,v_i+\sum_{i\le k}(u_i^Te/\sigma_i)\,v_i$,
y como $\{v_1,\dots,v_n\}$ es una base ortonormal de $\mathbb R^n$,
$x_{true}=\sum_{i=1}^n(v_i^Tx_{true})\,v_i$. Restando,
$$x_k-x_{true}=\sum_{i=1}^k\frac{u_i^Te}{\sigma_i}v_i-\sum_{i=k+1}^n(v_i^Tx_{true})\,v_i.$$
Las dos sumas usan conjuntos disjuntos de vectores ortonormales, así que por
el teorema de Pitágoras
$$\|x_k-x_{true}\|^2=\underbrace{\sum_{i=1}^k\frac{(u_i^Te)^2}{\sigma_i^2}}_{\text{ruido propagado}}+\underbrace{\sum_{i=k+1}^n(v_i^Tx_{true})^2}_{\text{sesgo por truncamiento}}.$$

- El ruido propagado no decrece al aumentar $k$; el sesgo no crece. Por eso
  el error no es monótono en $k$ y existe un $k$ que equilibra ambos términos.
- Pasar de $k$ a $k+1$ agrega $(u_{k+1}^Te)^2/\sigma_{k+1}^2$ y quita
  $(v_{k+1}^Tx_{true})^2$. Como $\sigma_{k+1}(v_{k+1}^Tx_{true})=u_{k+1}^T(Ax_{true})$,
  **agregar la dirección $k+1$ reduce el error si y solo si
  $|u_{k+1}^T(Ax_{true})|>|u_{k+1}^Te|$**, es decir, si en esa dirección la
  señal pesa más que el ruido. Con §2.5, el $k$ óptimo debe estar cerca del
  último índice en que $|u_i^T(Ax_{true})|$ supera $\delta/\sqrt m$.
- El sesgo incluye las direcciones $i>r$: es el piso de error de §1.5.

### 3.3 Residuo y norma

De la misma descomposición,
$$\|Ax_k-b\|^2=\sum_{i=k+1}^r(u_i^Tb)^2+\|b_\perp\|^2,\qquad\|x_k\|^2=\sum_{i=1}^k\left(\frac{u_i^Tb}{\sigma_i}\right)^2.$$
Al aumentar $k$, el residuo no crece y la norma no decrece: mejor ajuste a
cambio de una solución más grande. Ese es el compromiso que dibuja la curva L
(§5). El residuo por sí solo no sirve para elegir $k$, porque es mínimo en
$k=r$, que corresponde a la pseudoinversa numérica.

### 3.4 Ejemplo diagonal

Con $k=1$, $x_1=(1,\,0)^T$. El ruido propagado es
$(u_1^Te)^2/\sigma_1^2=0$ y el sesgo es $(v_2^Tx_{true})^2=1$, así que el
error es $1$, todo por sesgo: se descartó la segunda componente, cuyo valor
verdadero es 1. El residuo es $|u_2^Tb|=0.011$. Se acepta un peor ajuste a
cambio de una reconstrucción mucho más cercana a la verdad (error 1 frente a
10), como en [1, sec. 4].

### 3.5 Pseudoinversa numérica frente a truncamiento

Las dos son $x_k$: la pseudoinversa numérica toma $k=r=89$, fijado por la
tolerancia de máquina; la SVD truncada elige $k$ según el ruido. La diferencia
está en el criterio de corte, no en la fórmula. Esto responde otra pregunta de
defensa de [1, sec. 11].

Decisiones para el experimento:

- **Selección por índice.** Elegimos $k$ directamente, no un umbral sobre
  $\sigma_i$. Entre los valores retenidos no hay valores singulares repetidos
  (la menor brecha relativa entre consecutivos es 0.4 %, entre $\sigma_1$ y
  $\sigma_2$), así que $x_k$ está bien definido para cada $k$, como exige
  [1, sec. 4].
- **Barrido de $k$.** La lista actual de `experimento.py` usa nueve valores y
  dos de ellos (100 y 128) superan $r$, por lo que repiten la pseudoinversa
  numérica. Para ubicar el óptimo, el experimento debe recorrer todos los $k$
  de 1 a $r$ (pendiente de implementar).

## 4. Regularización de Tikhonov

### 4.1 El problema

Tikhonov busca un compromiso entre ajuste y magnitud [1, sec. 5]:
$$x_\lambda=\arg\min_x\left\{\|Ax-b\|^2+\lambda\|x\|^2\right\},\qquad\lambda>0.$$
Seguimos la convención de [1]: $\lambda$ es directamente el coeficiente de la
penalización, no $\lambda^2$.

### 4.2 Existencia y unicidad

**Teorema 5.1** [1, sec. 5]. Para todo $\lambda>0$ el problema tiene una única
solución,
$$x_\lambda=(A^TA+\lambda I)^{-1}A^Tb=\sum_{i=1}^r\frac{\sigma_i}{\sigma_i^2+\lambda}\,(u_i^Tb)\,v_i.$$

*Demostración de la unicidad.* Sea
$F(x)=\|Ax-b\|^2+\lambda\|x\|^2$. Su gradiente es
$\nabla F(x)=2A^T(Ax-b)+2\lambda x$ y su hessiana,
$2(A^TA+\lambda I)$. Para $z\ne0$,
$$z^T(A^TA+\lambda I)z=\|Az\|^2+\lambda\|z\|^2>0,$$
así que $A^TA+\lambda I$ es definida positiva. Por tanto es invertible, $F$
es estrictamente convexa y su único punto estacionario, la solución de
$(A^TA+\lambda I)x=A^Tb$, es su único minimizador global. $\blacksquare$

La clave está en el término $\lambda\|z\|^2$: es positivo incluso cuando
$z\in\ker A$, donde $\|Az\|^2=0$. Sin penalización ($\lambda=0$), cualquier
dirección del núcleo podía sumarse a una solución sin cambiar el residuo
(§1.1); con $\lambda>0$ esa ambigüedad desaparece. Esto responde la pregunta
de defensa sobre por qué $\lambda>0$ asegura unicidad [1, sec. 11].

### 4.3 Fórmula espectral

Con la SVD completa,
$A^TA+\lambda I=V(\Sigma^T\Sigma+\lambda I)V^T$ y $A^Tb=V\Sigma^TU^Tb$.
Escribiendo $x=Vy$ y multiplicando por $V^T$, el sistema se separa en
ecuaciones escalares:
$$(\sigma_i^2+\lambda)\,y_i=\sigma_i\,u_i^Tb\quad(i\le r),\qquad\lambda\,y_i=0\quad(i>r).$$
Luego $y_i=\sigma_iu_i^Tb/(\sigma_i^2+\lambda)$ para $i\le r$ y $y_i=0$ para
$i>r$: las componentes del núcleo se anulan porque no ayudan a ajustar los
datos pero sí aumentan la penalización. Esto da la fórmula del Teorema 5.1.

*Efecto de usar el rango numérico.* La implementación aplica la fórmula con
$r=89$ y omite las direcciones con $\sigma_i\le\tau_{\text{rango}}$. Cada una
de ellas aportaría un coeficiente de valor absoluto
$\sigma_i|u_i^Tb|/(\sigma_i^2+\lambda)\le(\sigma_i/\lambda)\,|u_i^Tb|$, así que
omitirlas cambia la solución en a lo más $(\sigma_{r+1}/\lambda)\,\|b_\perp\|$.
Para los $\lambda$ de la malla ($\lambda\ge10^{-10}\sigma_1^2$) y los
$\lambda$ óptimos del experimento, esa diferencia es despreciable.

### 4.4 Factores de filtro

La solución puede escribirse como la de la pseudoinversa numérica con cada
coeficiente multiplicado por un factor:
$$x_\lambda=\sum_{i=1}^rf_i(\lambda)\,\frac{u_i^Tb}{\sigma_i}\,v_i,\qquad f_i(\lambda)=\frac{\sigma_i^2}{\sigma_i^2+\lambda}\in(0,1).$$

- Si $\sigma_i\gg\sqrt\lambda$, entonces $f_i\approx1$: la dirección se
  conserva casi intacta.
- Si $\sigma_i\ll\sqrt\lambda$, entonces $f_i\approx\sigma_i^2/\lambda$: la
  dirección se atenúa fuertemente.
- Si $\sigma_i=\sqrt\lambda$, entonces $f_i=1/2$.

Así, $\sqrt\lambda$ funciona como un umbral suave sobre los valores
singulares. La SVD truncada es el caso de filtro abrupto, con $f_i=1$ para
$i\le k$ y $f_i=0$ para $i>k$ [1, sec. 5]. Esto predice que un $\lambda$
cercano a $\sigma_k^2$ se comporta de forma parecida a truncar en $k$.

### 4.5 Ejemplo diagonal

Con $\lambda=10^{-4}$, $\sqrt\lambda=0.01$ queda entre $\sigma_1=1$ y
$\sigma_2=10^{-3}$:
$$f_1=\frac{1}{1+10^{-4}}\approx0.99990,\qquad f_2=\frac{10^{-6}}{10^{-6}+10^{-4}}=\frac{1}{101}\approx0.0099.$$
Entonces
$$x_\lambda=\left(f_1\cdot\frac{1}{1},\ f_2\cdot\frac{0.011}{10^{-3}}\right)^T=\left(\frac{1}{1.0001},\ \frac{0.000011}{0.000101}\right)^T\approx(0.99990,\ 0.10891)^T,$$
como en [1, sec. 5]. El error es $\|(-0.0001,\,-0.89109)\|\approx0.8911$ y
el residuo, $\|((1-f_1)\cdot1,\ (1-f_2)\cdot0.011)\|\approx0.0109$.
Tikhonov no recupera el segundo valor verdadero, pero evita la enorme
amplificación del ruido: el error baja de 10 a 0.89.

### 4.6 Descomposición del error

Sustituyendo §2.1 en la fórmula de filtros,
$$x_\lambda-x_{true}=\sum_{i=1}^r\left[f_i\,\frac{u_i^Te}{\sigma_i}-(1-f_i)\,v_i^Tx_{true}\right]v_i-Px_{true},$$
con $P$ la proyección de §1.5. Por ortogonalidad,
$$\|x_\lambda-x_{true}\|^2=\sum_{i=1}^r\left[f_i\,\frac{u_i^Te}{\sigma_i}-(1-f_i)\,v_i^Tx_{true}\right]^2+\|Px_{true}\|^2.$$
Al aumentar $\lambda$, todos los $f_i$ disminuyen: el término de ruido
$f_iu_i^Te/\sigma_i$ se achica y el de sesgo $(1-f_i)v_i^Tx_{true}$ crece. A
diferencia de la SVD truncada, en cada dirección ambos efectos se combinan en
un solo término. Esto responde la pregunta de defensa sobre qué se gana y qué
se pierde al aumentar $\lambda$: se gana estabilidad y se pierde señal.

### 4.7 Estabilidad y sesgo

Si los datos cambian de $b$ a $b+e$, por linealidad
$$x_\lambda(b+e)-x_\lambda(b)=\sum_{i=1}^r\frac{\sigma_i}{\sigma_i^2+\lambda}(u_i^Te)\,v_i,$$
y por tanto
$$\|x_\lambda(b+e)-x_\lambda(b)\|\le\max_{s\ge0}\frac{s}{s^2+\lambda}\,\|e\|=\frac{\|e\|}{2\sqrt\lambda}\quad[1,\ \text{sec. }6].$$
El máximo se obtiene de $(s-\sqrt\lambda)^2\ge0$, que equivale a
$s^2+\lambda\ge2s\sqrt\lambda$, es decir, $s/(s^2+\lambda)\le1/(2\sqrt\lambda)$,
con igualdad en $s=\sqrt\lambda$. Mientras la pseudoinversa numérica
amplifica el ruido hasta por $1/\sigma_r\approx2.4\times10^{13}$, Tikhonov lo
limita a $1/(2\sqrt\lambda)$, que se controla eligiendo $\lambda$.

El precio es el sesgo: sin ruido, la componente $v_i^Tx_{true}$ se multiplica
por $f_i(\lambda)$ en lugar de por 1, y las componentes del núcleo siguen sin
recuperarse. Los límites confirman el compromiso: cuando $\lambda\to0^+$,
$f_i\to1$ y se recupera la pseudoinversa numérica; cuando
$\lambda\to\infty$, $f_i\to0$ y la solución tiende al vector cero. El valor
adecuado depende de la escala, el ruido y la señal, no de una constante
universal [1, sec. 6].

### 4.8 Sistema aumentado (verificación)

Para todo $x$,
$$\left\|\begin{pmatrix}A\\\sqrt\lambda\,I\end{pmatrix}x-\begin{pmatrix}b\\0\end{pmatrix}\right\|^2=\|Ax-b\|^2+\lambda\|x\|^2,$$
así que $x_\lambda$ es también la solución de mínimos cuadrados de este
sistema aumentado [1, sec. 9]. Su matriz $M$ cumple $M^TM=A^TA+\lambda I$, de
modo que sus valores singulares son $\sqrt{\sigma_i^2+\lambda}$: tiene rango
completo para $\lambda>0$ y
$$\kappa(M)=\sqrt{\frac{\sigma_1^2+\lambda}{\sigma_n^2+\lambda}}\le\sqrt{1+\frac{\sigma_1^2}{\lambda}}.$$
Incluso en el extremo de la malla, $\lambda=10^{-10}\sigma_1^2$, esto da
$\kappa(M)\le10^5$, frente a $\kappa_r\approx2.4\times10^{13}$ de $A$. Resolver
el sistema aumentado con un método estable es, por tanto, una forma
independiente y confiable de verificar la fórmula espectral.

*Sobre el método usado.* El código resuelve el sistema aumentado con
`np.linalg.lstsq`, que en LAPACK usa la rutina `gelsd`, basada en la SVD de
$M$ y no en QR. Es un solucionador estable, como admite [1, sec. 9], pero no
debe describirse como "resuelto por QR". Para que la verificación no dependa
de la SVD, conviene usar una factorización QR de $M$ (`np.linalg.qr`) y
resolver el sistema triangular resultante (pendiente de implementar).

## 5. Selección del parámetro sin conocer la verdad

### 5.1 Selección oráculo

En un experimento sintético podemos calcular
$$\lambda_{\text{or}}=\arg\min_{\lambda\in\text{malla}}\|x_\lambda-x_{true}\|,$$
y análogamente $k_{\text{or}}$ para la SVD truncada. Esta selección usa
$x_{true}$, que en un problema real no se conoce: es una referencia para
medir qué tan buena es una regla, no una regla aplicable [1, sec. 7]. Esto
responde la pregunta de defensa sobre qué información usa la selección
oráculo.

### 5.2 Principio de discrepancia

*Idea.* La señal verdadera tiene residuo $\|Ax_{true}-b\|=\|e\|=\delta$. Una
reconstrucción con residuo mucho menor que $\delta$ está ajustando el ruido.
Si se conoce $\delta$, se busca un $\lambda$ cuyo residuo sea del orden de
$\tau\delta$, con $\tau>1$ [1, sec. 7].

*La función residuo.* Como $b-Ax_\lambda=\sum_{i\le r}(1-f_i)(u_i^Tb)\,u_i+b_\perp$,
$$\rho(\lambda)^2=\|Ax_\lambda-b\|^2=\sum_{i=1}^r\left(\frac{\lambda}{\sigma_i^2+\lambda}\right)^2(u_i^Tb)^2+\|b_\perp\|^2.$$
Cada factor $\lambda/(\sigma_i^2+\lambda)$ es continuo y creciente en
$\lambda$, así que $\rho$ es continua y no decreciente, y es estrictamente
creciente si algún $u_i^Tb\ne0$. Sus límites son
$$\lim_{\lambda\to0^+}\rho(\lambda)=\|b_\perp\|,\qquad\lim_{\lambda\to\infty}\rho(\lambda)=\Big(\sum_{i=1}^r(u_i^Tb)^2+\|b_\perp\|^2\Big)^{1/2}=\|b\|.$$

*Existencia y unicidad.* Por el teorema del valor intermedio y la monotonía
estricta, si
$$\|b_\perp\|<\tau\delta<\|b\|,$$
existe un único $\lambda_\delta>0$ con $\rho(\lambda_\delta)=\tau\delta$. Esta
es la condición de que el umbral sea alcanzable que exige [1, sec. 7].

*Comprobación en nuestro experimento.* $b_\perp$ es casi solo ruido en las 39
direcciones descartadas. Con la distribución del ruido de §2.5, esperamos
$\|b_\perp\|\approx\sqrt{(m-r)/m}\;\delta=\sqrt{39/128}\;\delta\approx0.55\,\delta$,
y eso se observa en los tres niveles de ruido. Como
$0.55\,\delta<1.5\,\delta<\|b\|\approx4.68$, la condición se cumple para
$\tau=1.5$ y para cualquier $\tau>1$.

*La regla en una malla finita.* Por la monotonía,
$\{\lambda_j:\rho(\lambda_j)\le\tau\delta\}=\{\lambda_j:\lambda_j\le\lambda_\delta\}$.
La regla "el mayor $\lambda_j$ de la malla cuyo residuo no excede $\tau\delta$"
elige entonces $\lambda_g=\max\{\lambda_j\le\lambda_\delta\}$, que cumple
$$\frac{\lambda_\delta}{q}<\lambda_g\le\lambda_\delta,\qquad q=10^{10/39}\approx1.80,$$
donde $q$ es la razón entre puntos consecutivos de la malla logarítmica. La
regla en malla nunca regulariza más que la discrepancia exacta y queda a menos
de un paso de malla de ella, sin necesidad de resolver
$\rho(\lambda)=\tau\delta$. El conjunto no es vacío porque en el extremo
inferior de la malla $\rho\approx\|b_\perp\|<\tau\delta$. Además, la monotonía
justifica el algoritmo: recorrer la malla de mayor a menor y detenerse en el
primer $\lambda_j$ que cumple, que es el mayor (la versión anterior del
código no se detenía; ver `04_resultados_preliminares.md`, §4).

*Elección de $\tau$.* $\tau$ debe ser mayor que 1, porque incluso la señal
verdadera deja residuo $\delta$; pedir menos obliga a ajustar el ruido. Un
$\tau$ más grande regulariza más y aumenta el sesgo. Elegimos $\tau=1.5$ como
un margen moderado; es una decisión del grupo, y el efecto de cambiarla debe
reportarse con un análisis de sensibilidad (pendiente).

*Supuestos.* La regla necesita conocer $\delta$. En nuestro experimento lo
conocemos exactamente porque generamos el ruido; en un problema real habría
que estimarlo, y un error en esa estimación cambia el $\lambda$ elegido.
También supone que los errores del modelo no dominan al ruido [1, sec. 7].

### 5.3 Curva L

La curva L es el conjunto de puntos
$\big(\log\rho(\lambda),\ \log\|x_\lambda\|\big)$ para $\lambda>0$, con
$$\|x_\lambda\|^2=\sum_{i=1}^r\left(f_i(\lambda)\,\frac{u_i^Tb}{\sigma_i}\right)^2.$$
Como $f_i$ decrece con $\lambda$, la norma decrece mientras el residuo crece:
la curva va desde los $\lambda$ pequeños (residuo pequeño, norma enorme) hasta
los grandes (norma pequeña, residuo grande). La esquina, el punto de máxima
curvatura, marca el compromiso entre ambos. Localizarla es una heurística, no
una garantía [1, sec. 7], [3]. En este avance la curva está graficada y la
localización de la esquina queda pendiente.

### 5.4 Selección frente a evaluación

En todo el trabajo separamos dos momentos: primero se elige el parámetro con
una regla que no usa $x_{true}$ (discrepancia o curva L); después se evalúa
el resultado contra $x_{true}$ y se compara con el oráculo [1, sec. 7]. La
validación cruzada generalizada es otra regla posible y queda como
ampliación [3].

## Referencias

La numeración [n] remite a la lista única de fuentes, en formato IEEE, del
`README.md`.
