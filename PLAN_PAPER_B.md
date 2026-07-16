# PLAN PAPER B — Resonancias relativistas: escalamiento complejo del Dirac radial

**Título tentativo (reposicionado tras scout):** *Benchmark Dirac–Debye
resonances from a certified exact-structure spectral method: spin-symmetry
breaking of autoionization widths at strong coupling.*

**Tesis en una frase:** casar los dos motores validados del proyecto
(escalamiento complejo con θ-trayectorias + base γ exacta con pencil
desplazado) para producir tablas de referencia sistemáticas de
resonancias Dirac–Debye sobre la grilla (Z, λ) con precisión certificada
de 10+ dígitos, y demostrar la **ruptura de la simetría de espín de las
anchuras a acoplamiento fuerte**: los dobletes (p₁/₂, p₃/₂) que Shi et
al. encontraron casi degenerados en su régimen se separan dramáticamente
cuando Zα → 1 — la contraparte en Γ del desdoblamiento λc de FINDINGS
§11.

**Revista objetivo:** PRA (ruta subscription, $0). **Autoría:** decidir
con la regla del Paper A (César primero; coautores solo con rol real).

## Literatura verificada (scout 2026-07-16) — POSICIONAMIENTO

**Claims MUERTOS (no reclamar):**
- "Primeras resonancias Dirac–Yukawa": **Shi–Liu–Niu–Guo, Few-Body
  Syst. 55:135 (2014)** — CSM sobre Dirac radial + Yukawa, tablas de
  E_r Y Γ, límite no-rel verificado vs J-matrix, y comparación
  explícita de dobletes j = l±1/2 ("buena simetría de espín" en su
  régimen). Citarlo, reproducirlo y superarlo (patrón E-A1).
- "Método nuevo": CSM-Dirac tiene 20+ años — Šeba LMP 16:51 (1988,
  teoría), Balslev–Helffer AAM 13:186 (1992), Ivanov–Ho PRA 69:023407
  (2004, Stark relativista), Alhaidari PRA 75:042707 (2007, base
  Laguerre tridiagonal relativista CSM). La novedad es la precisión
  certificada + la física, no el escalamiento.

**Claims VIVOS (defendibles):**
1. Tablas de REFERENCIA sistemáticas sobre grilla (Z, λ) en
   parametrización atómica/plasma (Dirac–Debye electrostático) — la
   literatura Dirac–Debye (Bielińska-Wąż/Karwowski/Saha PRE 2004 y
   secuelas) cubre ligados y continuo, NUNCA anchuras de resonancia.
2. **El régimen Zα → 1**: Shi et al. trabajaron donde los dobletes son
   casi degenerados; el régimen de ruptura fuerte + trayectorias de
   resonancias al cruzar λc + comportamiento de umbral en +c² está
   intocado.
3. Precisión benchmark con doble método independiente y ensamblaje
   certificado (ellos: precisión de paper de física; nosotros: 10+
   dígitos adjudicables).

**Anclas de validación:** Shi et al. 2014 (¡obtener el texto — paywall
Springer; pedir a autores si hace falta!); Alhaidari–Bahlouli–
Abdelmonem JPA 41:032001 (2008) para el límite no-rel; Ivanov–Ho 2004
(Stark relativista) como cross-check en otra perturbación;
Ackad–Horbatsch PRA 78:062711 (2008) para resonancias supercríticas
(potencial distinto, mismo instrumento a prueba).

---

## Fase F0 — Factibilidad y teoría del escalado (1 sesión) — CON COMPUERTA

El operador escalado, r → re^{iθ}: TODOS los términos de O± escalan
uniformemente (d/dr → e^{−iθ}d/dr, κ/r → e^{−iθ}κ/r), así que

    h_θ = [[V_θ + c², c e^{−iθ}O₊],[c e^{−iθ}O₋, V_θ − c²]],
    V_θ(r) = −Z e^{−λ r e^{iθ}} e^{−iθ}/r   (Yukawa escalado).

El término de masa NO se escala. Matrices complejo-SIMÉTRICAS (no
hermitianas): eig genérico, producto bilineal (transpuesta, no
conjugada).

- **Hipótesis H-B1 (verificar antes de nada):** la ecuación indicial en
  r→0 es invariante bajo el escalado (todos los términos singulares
  comparten el factor e^{−iθ}) ⇒ **γ no cambia y la base γ sigue siendo
  exacta en la cúspide bajo escalamiento complejo.** Si se confirma, es
  un selling point del paper (exact structure survives complex
  scaling); verificarla simbólica y numéricamente.
- E-B0a: continuo rotado — espectro de h_θ libre discreto: dos ramas
  desde ±c² con ángulo −2θ en la variable p² (en binding: rama
  electrónica desde 0 como el no-rel; rama positrónica desde −2c²).
  Compuerta: los conos salen limpios y estables en nb.
- E-B0b: estados ligados θ-independientes (teorema de Balslev–Combes
  discreto): binding de 1s/2s/2p a Z=80 invariante en θ∈[0.05,0.4] a
  ~1e-12 (con el pencil desplazado complejo: ensamblar h_θ − c²S).
- E-B0c: certificaciones complejas — C1 (overlap real, no cambia),
  C4 versión bilineal ⟨f|O₊g⟩ = ⟨O₋f|g⟩ con transpuesta. Definir
  tolerancias para aritmética compleja.

**COMPUERTA G0:** si la rama positrónica contamina la ventana
electrónica (autovalores espurios complejos mid-gap) o los ligados no
son θ-estables → detenerse y diagnosticar con la anatomía T1 antes de
seguir. Criterio de aborto de fase: >2 sesiones sin ligados θ-estables.

## Fase F1 — Primera resonancia y anclaje no-relativista (1 sesión)

- E-B1a: la 2p de Yukawa a Z=1, λ=0.25: límite no-relativista conocido
  de NUESTRO propio motor validado (E_r=0.00306445, Γ=4.369e-3,
  triple-verificado). La versión Dirac (κ=+1 y −2) debe reproducirlo a
  O((Zα)²) ~ 5e-5 relativo. θ-trayectorias: buscar estacionariedad
  |dE/dθ| mínima, protocolo idéntico al no-rel.
- E-B1b: barrido β y nb: independencia de parámetros no lineales
  (mismo estándar que E7 del motor no-rel).
- **Nota κ>0:** el espurio de Tupitsyn–Shabaev ahora es un autovalor
  complejo cerca del 1s del canal compañero — identificarlo y
  documentar su firma compleja (bonus de anatomía, conecta con Paper T).

**COMPUERTA G1:** coincidencia no-rel a mejor que 1e-4 relativo en E_r
y Γ. Si no: revisar convención de escalado del potencial antes de
producir nada.

## Fase F2 — Validación independiente (1-2 sesiones) — LA FASE CRÍTICA

La validación externa es escasa (confirmar con scout de literatura; a
la fecha del plan no conocemos tablas Dirac–Yukawa de resonancias).
Mitigación con DOS métodos independientes internos:

- E-B2a: **shooting de onda saliente en el plano complejo**: integrar la
  ODE radial de Dirac con condición de frontera saliente (momento
  complejo de la rama electrónica) y cazar el polo de la matriz S por
  Newton en E complejo. Método totalmente independiente de la base
  (integrador ODE + análisis complejo). 3-5 resonancias a 10+ dígitos
  como anclas internas.
- E-B2b: límite λ→λc⁺: la resonancia debe nacer del estado ligado que
  muere en λc(Z,κ) — CONSISTENCIA con el mapa E-A3 (FINDINGS §11):
  E_r(λ→λc) → 0 y Γ → 0 con la ley de umbral correcta (p=2 en ℓ=0,
  p=1 en ℓ=1 — ya la conocemos).
- E-B2c: si el scout encuentra CUALQUIER número publicado (CAP-Dirac,
  estabilización relativista, R-matrix), reproducirlo.

**COMPUERTA G2:** acuerdo base-vs-shooting ≥ 10 dígitos en ≥3 casos.
Ese es el estándar "certificado" del proyecto; sin esto no hay paper.

## Fase F3 — Producción de física (2 sesiones)

- E-B3a: **tabla de referencia** E_r − iΓ/2 para Dirac–Yukawa:
  Z ∈ {1, 20, 40, 60, 80, 92}, κ ∈ {−1, +1, −2}, 2-3 λ por caso
  (protocolo de dígitos honestos: deriva en nb/θ/β como barra de error).
- E-B3b: **el resultado insignia — desdoblamiento de estructura fina de
  las anchuras**: Γ(2p₁/₂) vs Γ(2p₃/₂) al mismo λ y Z alto. Predicción
  cualitativa del mapa λc: el 2p₁/₂ (más resistente al apantallamiento)
  debería tener Γ menor cerca del umbral común. Cuantificarlo vs Zα.
- E-B3c: trayectorias de resonancias en el plano complejo al barrer λ
  (figura estrella: familias κ separándose con Zα).
- E-B3d (opcional, si el tiempo da): FEAST complejo para ventanas
  espectrales — solo si nb crece tanto que eig denso moleste (con el
  techo de cuadratura eliminado podemos ir a nb~1000; eig complejo
  2000×2000 sigue siendo trivial, así que probablemente FEAST NO es
  necesario aquí — no forzarlo, honestidad ante todo).

## Fase F4 — Verificación adversarial + manuscrito (1-2 sesiones)

- Workflow de verificación (patrón E-A3): un refutador numérico con
  protocolo independiente + un verificador del shooting.
- Manuscrito REVTeX/PRA reutilizando la infraestructura del Paper A;
  figuras con el estilo de casa; FINDINGS como fuente.

## Riesgos y mitigaciones

| Riesgo | Señal | Mitigación |
|---|---|---|
| Estados salientes mal representados por base real | Γ no converge en nb | β complejo (escalado exterior suave) o θ mayor; el motor no-rel ya pasó por esto |
| Crowding del continuo rotado sobre resonancias anchas | trayectorias θ sin punto estacionario | ventanas espectrales + filtrado por residuo; subir θ; mostrar θ-estabilidad POR CADA estado tabulado (exigencia del scout) |
| Cancelación de masa en espectro complejo | ruido 1e-11 en Γ pequeñas | pencil desplazado complejo (h_θ − c²S) desde el día 1 |
| Polución no-hermitiana nueva | autovalores complejos sin θ-traza | anatomía T1 extendida (¡material para Paper C/T v2!); Zaytsev et al. PRA 100:052504 documenta balanceo anti-espurios en CSM relativista |
| s₁/₂ se vuelve estado virtual cerca del umbral (sin barrera centrífuga) | la "resonancia" κ=−1 no aparece | claims de splitting SOLO en pares (p₁/₂,p₃/₂) y (d₃/₂,d₅/₂) — resonancias de forma genuinas |
| Terminología de comunidad nuclear | referee del área Guo | enganchar explícitamente el análisis de simetría de espín de Shi et al. y la literatura pseudospin (Lu–Zhao–Zhou PRL 109:072501) |
| Tablas de Shi et al. inaccesibles (paywall) | sin ancla directa | acceso institucional ESPOCH / mail a autores; hedge: límite no-rel Alhaidari 2008 + ligados Dirac–Debye PRE 69:016404 |
| Z > 137 mal definido a núcleo puntual | γ imaginario | tablas hasta Z ≤ 137; NO entrar al supercrítico (colisiona con Ackad–Horbatsch y cambia el carácter del paper) |

## Presupuesto y dependencias

- **Total estimado: 6-8 sesiones.** Depende de: nada externo. El código
  parte de `src/dirac.py` (rama compleja nueva: `dirac_cs.py`) +
  `find_resonance` del motor no-rel.
- Entregables por fase quedan en `exploration/dirac_cs/` con FINDINGS
  propio, mismo protocolo de certificación + pre-registro de
  predicciones que nos ha salvado ya cuatro veces.
