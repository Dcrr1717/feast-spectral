# TRASPASO — Paper NR (CompAMa) — 27 jul 2026

Pega el bloque **PROMPT** de abajo en la otra sesión de Claude Code, abierta en
`C:\Users\dcrr1\repos\feast-spectral`. Todo lo demás de este archivo es
contexto que el prompt ya referencia.

---

## ESTADO ACTUAL (verificado)

- Manuscrito: `C:\Users\dcrr1\repos\feast-spectral\paper\main.tex`
- Compila limpio: **14 páginas**, 0 overfull vbox, 0 em dashes, 0 refs
  indefinidas, 0 `Figure~\ref`/`Table~\ref` (todo en forma corta Fig./Tab.).
- Respaldo de la versión previa de 15 páginas: `paper\main.tex.bak-15pag`
- Último commit limpio en GitHub: `8119a25`. Los cambios de hoy posteriores a
  ese commit **no están commiteados todavía**.

### OJO — dato crítico sobre el conteo de páginas
Hubo un momento en que el PDF decía 13 páginas, pero **era falso**: se lograba
porque una página se desbordaba (`Overfull \vbox 28.85pt`) y el número de
página quedaba **impreso encima del pie de la Fig. 4**. Al corregir el
desborde el conteo honesto pasó a 14. **No aceptar 13 páginas sin verificar
que `grep -c "Overfull ..vbox" main.log` dé 0.**

---

## LO QUE YA SE HIZO HOY (no rehacer)

Reunión con el coautor Prof. Raúl Manzanilla (Yachay Tech) → 10 correcciones
aplicadas y verificadas:

1. Método renombrado **FEASTci** (nuestra implementación) vs **FEAST**
   (el algoritmo de Polizzi). Título, encabezados, resúmenes, keywords, Sec. 2.4.
2. FEAST explicado como *algoritmo*, no como paquete.
3. Filiales sin dirección postal (solo institución + país).
4. Marcas de autor numéricas (1, 2) en vez de símbolos.
5. Tabla de notación con columna de **Unidades**.
6. Tabla de **procedencia de métodos** (estándar / implementado por nosotros /
   modificado por nosotros).
7. β, estados espectroscópicos, Ha y a₀ explicados para no especialistas.
8. `Ref.~\cite{}` reemplazado por atribuciones descriptivas.
9. Pies de figura reducidos a títulos de una línea; la explicación se movió
   al texto corrido.
10. `wave function` separado; ℝ con `\mathbbm` (decisión del coautor, ver abajo);
    sumatoria con `\displaystyle`; sin em dashes; figuras regeneradas.

Decisiones ya tomadas por el usuario (NO revertir):
- **ℝ se escribe con `\mathbbm{R}`** (preferencia del profesor) aunque las
  instrucciones oficiales del journal piden `\mathbb{R}`. El usuario dijo:
  "hagámosle caso al profe, si me mandan a corregir yo corrijo luego".
  El macro `\realset` en el preámbulo permite cambiarlo en una línea.
- Las dos figuras de convergencia van **separadas** (ℓ=0 y ℓ=1), no unidas.
- No achicar figuras por debajo de 0.85\linewidth.
- El profesor relajó el límite a **12–13 páginas**.

---

## HALLAZGOS CONFIRMADOS PENDIENTES DE APLICAR

Salieron de una ronda adversarial de 4 revisores + refutación. Estos 7
**sobrevivieron** el intento de refutación.

**LOG COMPLETO** (con el texto LaTeX exacto de reemplazo de cada hallazgo,
en el campo `result.confirmed`), guardado en tres ubicaciones permanentes:
- `C:\Users\dcrr1\repos\feast-spectral\paper\LOG_diagnostico_hallazgos.json`  ← usar esta
- `C:\Users\dcrr1\OneDrive\Desktop\ReunionManzanilla\LOG_diagnostico_hallazgos.json`
- `C:\ReunionManzanilla\LOG_diagnostico_hallazgos.json`

### BLOCKER 1 — Sobreclamo de simetría compleja (el más grave)
El paper dice que la variante de contorno circular explota que H(θ) es
complejo-simétrica para reemplazar el par biortogonal de Kestyn. **El código
no hace eso**: `src/feast.py:128` usa `Hr = W.conj().T @ H @ W` (forma
sesquilineal, válida para cualquier matriz), no la forma bilineal `W^T H W`.
La simetría compleja no se usa en ningún lado. Es el claim de novedad principal
y está sobreclamado. → Reescribir la prosa para describir lo que el código
realmente hace; el aporte honesto es **la colocación del contorno fuera del eje
real sobre la resonancia anticipada**.

### BLOCKER 2 — "Nueve cifras" sin la restricción ℓ≥1 en los resúmenes
Ambos resúmenes (EN y ES) dicen que los parámetros críticos coinciden a nueve
cifras, sin decir "para ℓ≥1". La Tabla 4 del propio paper muestra 1s con
diferencia 6.0e-4 (cuatro cifras). Sec. 4.2 y conclusiones sí lo restringen.

### BLOCKER 3 — La certificación 1e-12 es falsa y no se ejecuta
El paper dice que max|S−I| "se mantiene por debajo de 1e-12 en todas las
corridas de producción" y que se verifica "at run time". Los datos del propio
repo (`results/E1_hydrogen_convergence.csv`, columna `overlap_err`) dan
1.303e-12 en N=150 y **1.510e-12 en N=180**, que es la base de producción de
la Tabla 3. Además `BasisSet.overlap_error()` solo se llama desde
`scripts/run_validation.py`; ni `run_resonances.py` ni `run_benchmarks.py` ni
`run_sweep.py` lo llaman, y `__init__` no lo verifica: **no existe el gate en
tiempo de ejecución que el resumen promete**. → Arreglar código (añadir el
gate real en `basis.py.__init__`, tolerancia 1e-9, el peor caso medido es
9.7e-12) **y** el texto (citar los valores medidos por tabla).

### MAJOR 4 — "Nueve cifras" tampoco es exacto para ℓ≥1
Solo 2p llega a nueve cifras; 3p y 3d llegan a ocho. Aparece 4 veces.
→ Reemplazar por una tolerancia con unidades: "dentro de 2e-9 a₀⁻¹, es decir
ocho cifras significativas".

### MAJOR 5 — Bloque de prueba declarado real y descrito como complejo
`Y \in \realset^{N×m_0}` seguido de "complex for the circular-contour variant"
se contradice. El código usa bloque genuinamente complejo en `feast_region`.
→ Añadir macro `\complexset` y declarar `Y \in \complexset^{N×m_0}` con la
aclaración de que entradas reales bastan para la variante de intervalo.

### MAJOR 6 — Piso de precisión contradice la Fig. 2
El texto dice que las curvas se aplanan en ~1e-10 Ha, pero la Fig. 2 (ℓ=1)
se aplana entre 1e-15 y 1e-13. Además "para ℓ=1 los errores llegan a 1e-14"
no se cumple con N=150, β=1/6 (el CSV da 2.644e-13).
→ Escribir el piso como dependiente del canal.

### MAJOR 7 — No existe el "explicit residual gate"
En `feast_region` el residuo solo es criterio de parada, no filtro de
aceptación (la selección es puramente geométrica, `feast.py:131`). La variante
de intervalo hace lo mismo, así que ni siquiera distingue.
→ Quitar "explicit residual gate" o implementarlo de verdad.

### REFUTADO (no tocar)
"La tabla de procedencia asigna mal el rol de la variante circular" — se
verificó que sí produce números reportados.

---

## HALLAZGOS SIN VERIFICAR (la refutación murió por límite de gasto)
Existen pero **nadie los atacó**, así que pueden ser falsos positivos.
Verificar antes de aplicar:
- [presentación] Las 4 figuras mezclan dos familias tipográficas
- [presentación] Fig. 4: la caja de leyenda sobresale del eje
- [presentación] Fig. 3: rótulos del inset demasiado pequeños
- [presentación] Tablas 3, 4 y 5 sin unidades en los encabezados de columna
- [fidelidad] Radio de contorno FEASTci declarado vs código
- [fidelidad] "Exacto a 1e-10 Ha o mejor" en Tabla 3
- [fidelidad] Acuerdo declarado con los benchmarks
- [fidelidad] Filas ℓ=0 de la Tabla 4
- [claridad] "Entradas de Tab. 3 exactas a…"

---

## PROMPT (pegar en la otra sesión)

```
Trabajo en C:\Users\dcrr1\repos\feast-spectral sobre el manuscrito
paper\main.tex (revista Bulletin of Computational Applied Mathematics).
Lee primero, completo, el archivo
C:\Users\dcrr1\OneDrive\Desktop\ReunionManzanilla\TRASPASO_OTRA_CUENTA.md
que tiene el estado, las decisiones ya tomadas y los hallazgos pendientes.

Los detalles completos de cada hallazgo, con el texto LaTeX exacto de
reemplazo propuesto, están en este JSON (campo result.confirmed):
C:\Users\dcrr1\repos\feast-spectral\paper\LOG_diagnostico_hallazgos.json
Léelo completo antes de editar nada.

TAREA, en este orden:

1. Aplica los 3 BLOCKERS y los 4 MAJORS confirmados. Para cada uno verifica tú
   mismo la evidencia contra el código en src/ y los datos en results/ antes de
   editar: los revisores pueden equivocarse y tú eres responsable del resultado.
   El Blocker 3 requiere tocar src/basis.py además del texto.

2. Verifica los 9 hallazgos NO verificados que lista el traspaso. Descarta los
   que no se sostengan; aplica los reales.

3. Deja el paper en 12 o 13 páginas. CRÍTICO: verifica siempre con
   grep -c "Overfull ..vbox" main.log  que dé 0. Un 13 obtenido desbordando la
   página es un fallo, no un éxito: el número de página termina impreso encima
   del texto. Si hay que elegir, prefiere 14 páginas correctas a 13 rotas, y
   dime que no se pudo.

4. Restricciones que NO puedes violar:
   - Nunca cambies un valor numérico de resultados.
   - Prohibido el em dash "---": usa comas, dos puntos, punto y coma o paréntesis.
   - Referencias en prosa siempre en forma corta: Fig.~\ref, Tab.~\ref,
     Sec.~\ref, Eq.~(\ref).
   - Los pies de figura son títulos de una línea; la explicación va en el texto.
   - Las dos figuras de convergencia van separadas, nunca unidas, y ninguna
     figura por debajo de 0.85\linewidth.
   - ℝ se queda con \mathbbm{R} (decisión del coautor), aunque el journal pida
     \mathbb{R}.
   - No agradezcas a CIDED ni a ninguna institución: el usuario lo pidió
     explícitamente porque no aportaron.

5. SKILLS: úsalas a tu libre discreción, y también workflows/subagentes si te
   ayudan. Las que encajan aquí:
   - scientific-writing  -> para reescribir los pasajes de los hallazgos.
   - math_daniel (o math-rigor) -> úsala SÍ O SÍ para el Blocker 1. El claim
     en disputa es matemático ("como H(theta) es complejo-simétrica, sus
     subespacios invariantes izquierdo y derecho coinciden salvo conjugación,
     luego basta un subespacio derecho"). Hay que decidir con rigor si la
     premisa es cierta, si la conclusión se sigue, y si el código la usa o no.
   - humanizer_academic (o avoid-ai-writing) -> el manuscrito recibió varias
     pasadas de edición asistida hoy; conviene un barrido de "AI-isms" antes
     de enviar.
   - academic-paper-reviewer -> opcional, como pasada final de árbitro.
   - latex-guardian -> útil por sus chequeos de compilación y estructura, PERO
     exige confirmación del usuario antes de cada edición .tex; no la actives
     si vas a trabajar en modo autónomo sin él presente.

6. Al terminar: compila dos veces, reporta el conteo de páginas, confirma las
   verificaciones mecánicas (overfull vbox, em dashes, refs cortas, refs
   indefinidas), y haz commit + push describiendo lo aplicado y lo rechazado
   con su razón. Sincroniza también la carpeta
   C:\Users\dcrr1\OneDrive\Desktop\ReunionManzanilla\Paper_NR_FEAST_fuente\
   (main.tex, main.pdf, figures/) y regenera ahí el zip
   Paper_NR_FEAST_fuente_27-7-2026.zip.

Contexto del usuario: César Reinoso, matemático de ESPOCH. Háblale en español,
directo, sin adular. El paper se envía a CompAMa y el coautor Manzanilla ya
revisó; estos son los arreglos de la ronda de diagnóstico posterior.
```

---

## ARCHIVOS CLAVE
- Manuscrito: `repos\feast-spectral\paper\main.tex`
- Respaldo 15 pág: `repos\feast-spectral\paper\main.tex.bak-15pag`
- Código: `repos\feast-spectral\src\{feast,basis,reference_fd}.py`
- Datos: `repos\feast-spectral\results\`
- Generador de figuras: `repos\feast-spectral\scripts\replot_paper_figures.py`
- Instrucciones oficiales del journal: en el scratchpad,
  `CompAMa_original\Instructions_journal0.tex` (también hay copia en el zip
  `CompAMa_plantilla_oficial.zip` del Escritorio)
- Carpeta de entrega: `Escritorio\ReunionManzanilla\Paper_NR_FEAST_fuente\`
