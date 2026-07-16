# PLAN PAPER B — Resonancias relativistas: escalamiento complejo del Dirac radial

**Título tentativo:** *Complex-scaled exact-structure spectral method for
relativistic resonances: reference Dirac–Yukawa widths and the
fine-structure splitting of autoionization by screening.*

**Tesis en una frase:** casar los dos motores validados del proyecto
(escalamiento complejo con θ-trayectorias del motor no-relativista +
base γ exacta con pencil desplazado del motor Dirac) para producir las
primeras tablas de referencia de resonancias Dirac–Yukawa, incluyendo
la contraparte en anchuras Γ del desdoblamiento de estructura fina que
descubrimos en el mapa λc (FINDINGS §11).

**Revista objetivo:** PRA (ruta subscription, $0). **Autoría:** decidir
con la regla del Paper A (César primero; coautores solo con rol real).

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
| Crowding del continuo rotado sobre resonancias anchas | trayectorias θ sin punto estacionario | ventanas espectrales + filtrado por residuo; subir θ |
| Cancelación de masa en espectro complejo | ruido 1e-11 en Γ pequeñas | pencil desplazado complejo (h_θ − c²S) desde el día 1 |
| Polución no-hermitiana nueva | autovalores complejos sin θ-traza | anatomía T1 extendida (¡material para Paper C/T v2!) |
| Scoop 2026 | scout encuentra competidor directo | reposicionar a "benchmarks + splitting" (los dígitos y el splitting siguen siendo nuestros) |

## Presupuesto y dependencias

- **Total estimado: 6-8 sesiones.** Depende de: nada externo. El código
  parte de `src/dirac.py` (rama compleja nueva: `dirac_cs.py`) +
  `find_resonance` del motor no-rel.
- Entregables por fase quedan en `exploration/dirac_cs/` con FINDINGS
  propio, mismo protocolo de certificación + pre-registro de
  predicciones que nos ha salvado ya cuatro veces.
