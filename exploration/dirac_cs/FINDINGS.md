# FINDINGS — Paper B: Dirac complejo-escalado (dirac_cs)

## F0 — Factibilidad (2026-07-16) · COMPUERTAS G0a/G0b/G0c: VERDES

**Ensamblaje:** el pencil desplazado complejo se construye con las
MISMAS matrices reales certificadas del solver Dirac:
H_θ − c²S = [[e^{−iθ}V_d, ce^{−iθ}K],[ce^{−iθ}Kᵀ, e^{−iθ}V_d − 2c²S]]
(complejo-SIMÉTRICO; Coulomb escala exacto con e^{−iθ}; Yukawa añade
e^{−λre^{iθ}} dentro de la cuadratura convergida). θ=0 reproduce el
solver real a 2.7e-11 (P-F0.3 ✓).

**H-B1 confirmada + teorema gratis (P-F0.2, G0a ✓):** la ecuación
indicial es invariante (todo término singular comparte e^{−iθ}) y —más
fuerte— **la fibración SVD de L1' sobrevive al escalado**: la fibra es
h_θ(σ) = [[c², cσe^{−iθ}],[cσe^{−iθ}, −c²]] y el espectro libre
discreto cae EXACTO en ±√(c⁴+c²σ_j²e^{−2iθ}) (medido: 3.2e-10 absoluto
= ~2e-14 relativo a c²). Los dos conos rotados desde ±c² salen de la
misma álgebra. Toda la maquinaria del paper del teorema se extiende al
caso complejo-escalado — remark para Papers B y T.

**Balslev–Combes discreto (P-F0.4, G0b ✓):** 1s de Z=80 (base γ,
nb=100): parte real = Sommerfeld exacto en todos los dígitos impresos
(−3532.192093491), |Im| ≤ 6e-11, deriva ≤ 1.4e-10 en θ ∈ [0.05, 0.4].

**G0c ✓ (adjudicado en F1):** el único punto sospechoso
(−33.09−0.32j a nb=100) es el ligado n=10 (exacto: −33.151) deformado
por truncación de base real en estado difuso: a nb=160 converge a
−33.1506+0.0000j. NO es espurio mid-gap. Lección de instrumento: los
Rydberg altos adquieren ancho falso bajo rotación si nb no los
resuelve — el detector de espurios debe comparar |Im| contra el error
de truncación por estado.

## F1 — Ancla no-relativista (2026-07-16) · COMPUERTA G1: VERDE

Resonancia 2p de Yukawa, Z=1, λ=0.25, contra NUESTRO benchmark no-rel
triple-verificado (E_r = 0.00306445, Γ = 4.369e-3):

| canal | E_r | Γ | dE_r/E_r | dΓ/Γ | ‖dE/dθ‖ |
|---|---|---|---|---|---|
| 2p₁/₂ (κ=+1) | 0.00306432 | 4.368987e-3 | −4.3e-5 | −3.1e-6 | 1.2e-11 |
| 2p₃/₂ (κ=−2) | 0.00306443 | 4.369406e-3 | −5.6e-6 | +9.3e-5 | 5.1e-12 |

Mesetas θ-estacionarias limpias (1e-11). Las desviaciones vs no-rel
tienen la magnitud O(α²) y el SIGNO correcto de la física relativista:
el 2p₁/₂ se corre más (−4.3e-5) que el 2p₃/₂ (−5.6e-6) — no son error,
son el corrimiento relativista asomando. **Desdoblamiento a Z=1:
dE_r = −1.1e-7, dΓ = −4.2e-7** — la semilla de la ruptura de simetría
de espín que debe crecer dramáticamente hacia Zα → 1 (el resultado
insignia planificado, F3).

**Siguiente (F2, la fase crítica):** shooting complejo independiente a
10+ dígitos + consistencia λ→λc con el mapa E-A3 + reproducir las
tablas de Shi–Liu–Niu–Guo 2014 cuando se consiga el texto.

## F2 — Validación independiente (2026-07-16) · COMPUERTA G2: VERDE

**Método independiente (f2_shooting.py):** shooting de Gamow sobre la
ODE radial de Dirac — arranque indicial + serie de primer orden en
r₀=1e-8 (error O((2cr₀)²)~1e-11), integración en s=ln r hasta r=1
(mata la rigidez 1/r) y en r hasta R; condición de resonancia =
anulación del coeficiente ENTRANTE, extraído como wronskiano con el
**spinor libre saliente exacto de Riccati–Hankel** (G derivado de la
propia ODE — sin adivinar signos); secante compleja.

**Dos bugs cazados por las compuertas en el primer intento:**
(1) el matching con onda plana pura ignora el centrífugo κ/(pR)~0.1 a
R finito → desviaciones 0.4-1.8; la cura es Riccati–Hankel exacto.
(2) el selector de la base por ventana+mínimo agarró un ARTEFACTO
ESTÁTICO (E_r idéntico en todo λ — punto del continuo rotado): el
selector correcto es por **θ-estabilidad** (autovalor que NO se mueve
entre dos ángulos de rotación).

**P-F2.1 ✓ (G2: ≥3 casos a ≥10 dígitos relativos):**
| caso | base | shooting | dev rel |
|---|---|---|---|
| Z=40 κ=+1 | 2.941245536171−1.355551909824j | ...536305−...909813j | 4.2e-11 |
| Z=40 κ=−2 | 3.413196281651−1.806500986794j | ...281388−...986677j | 7.4e-11 |
| Z=80 κ=−2 | 12.673150466581−6.188018766904j | ...466566−...766901j | **1.1e-12** |
| Z=80 κ=+1 | 0.848891759808−0.079848876487j | ...757494−...877035j | 2.8e-9 |
| Z=1 (ambos) | — | — | 1.6e-10 ABSOLUTO (|E|~0.004: 8 dígitos rel) |
R-estabilidad del shooting: 1e-10..8e-14 (P-F2.2 ✓).

**P-F2.3 ✓ nacimiento en el umbral (Z=40, κ=−2):** E_r crece desde
0.027 (s=0.221) linealmente; intercepto s_c = 0.2207 vs 0.2160 de E-A3
(consistente con el ~2% de incertidumbre 2p documentado en FINDINGS
§11); **exponente de Wigner medido 1.530 vs 1.5 teórico (ley p-wave)**.

**LA FÍSICA INSIGNIA APARECIÓ EN F2 (anticipo de F3):** al MISMO
apantallamiento s=0.24:
- Z=40 (Zα=0.29): Γ(2p₁/₂)=2.71 vs Γ(2p₃/₂)=3.61 → split de anchuras 25%
- Z=80 (Zα=0.58): E_r = 0.849 vs 12.67; **Γ = 0.16 vs 12.38 → −99%**
La ruptura de simetría de espín de la autoionización es de DOS ÓRDENES
DE MAGNITUD a acoplamiento fuerte — el compañero resistente (2p₁/₂,
mayor λc) apenas cruzó su umbral mientras el 2p₃/₂ ya es ancho. El
desdoblamiento λc de E-A3 hecho dinámica. (Números preliminares de un
solo (θ,β); F3 les pone el protocolo de dígitos honestos.)

**Siguiente (F3):** producción — tabla (Z, κ, λ) con protocolo de
dígitos honestos, el mapa del splitting vs Zα, trayectorias en el plano
complejo; y conseguir Shi et al. 2014 (pendiente operativo de César).

## E-B2c — Reproducción y adjudicación de Shi et al. 2014 (2026-07-16)

**El paper (conseguido por César):** Shi–Liu–Niu–Guo, Few-Body Syst.
55:135 (2014). Convenciones IDÉNTICAS a las nuestras (ℏ=M=1, c físico) —
cero fricción de unidades. Sus casos: dobletes de ℓ ALTO (κ=+5/−6 [ℓ=5],
+10/−11 [ℓ=10]) con acoplamientos fuertes A=110–200 (A/c=0.8–1.46),
base de oscilador armónico, θ=1.0. Por eso ven "buena simetría de
espín": el ℓ alto la protege. Nuestro régimen atómico 2p a Zα→1 queda
intacto como novedad. Reportan 4 dígitos.

**Reproducción (eb2c_shi_reproduction.py, doble método):**
- **8/12 estados reproducidos a su redondeo** (dev 1.4e-5..5.4e-5),
  con base-vs-shooting a 1e-9..5e-10: sus 4 dígitos → benchmarks
  nuestros de 9-10 dígitos. Ej.: su 2.2834−i0.0549 → nuestro
  2.283366856−i0.054871772 (shooting concuerda a 1.5e-9).
- **ADJUDICACIÓN — los estados angostos a A=170 están mal en su
  tabla:** su 0.4648−i0.0066 es nuestro 0.479288180−i0.000051244
  (doble método a 1.5e-11); su 0.6613−i0.0107 es 0.669863069
  −i0.000331590 (9.9e-11). Posiciones corridas ~1.5e-2 y anchuras
  infladas ~100×. El patrón se auto-delata en su propio paper: su
  límite no-rel también discrepa del J-matrix en los angostos (p.ej.
  anchura 0.0168 vs 1.09e-6 del J-matrix). Diagnóstico: base HO
  sub-convergida para estados angostos cerca del umbral — exactamente
  el régimen donde la estructura exacta gana (L3b en acción).
- **Añadimos las anchuras que ellos no resuelven:** sus "−0.0000" en
  A=200 son Γ = 3.04e-7 (κ=10) y 5.04e-7 (κ=−11).
- Limitación honesta del shooting: para resonancias MUY anchas
  (Γ ≳ 4) la extracción del coeficiente entrante muere en doble
  precisión (e^{2|Im p|R}); esos estados quedan validados por
  θ-estabilidad + acuerdo a 4 dígitos con Shi.

**Para el paper B:** sección de validación externa lista — reproducción
+ adjudicación + extensión de dígitos, el sello de la casa. Contacto
con Guo: estrategia de dos etapas (mail_guo_draft.md) — pospuesto por
decisión de César hasta terminar la producción.

## F3 — Producción: tabla y mapa del splitting (2026-07-16)

**Tabla (f3_production.py → f3_table.csv):** doblete 2p en
Z ∈ {1..130}, s ∈ {0.24, 0.30}, protocolo de dígitos honestos (grilla
2θ×2nb×2β; spreads 1e-11..1e-10 en todos los estados limpios) + ancla
de shooting donde Γ ≲ 4. Corrección en vivo: los estados MUY anchos
(Γ/E_r ≳ 2) requieren θ = (0.65, 0.80) — el par (0.45, 0.60) no los
expone y la fila Z=1/s=0.30 original era un artefacto de umbral
(reemplazada; el CSV consolidado tiene los valores buenos). Nota de
instrumento: el ancla de shooting en la tabla usó R=18/λ (cola ~1e-8);
para anclas finales usar R=30/λ como en F2.

**EL MAPA INSIGNIA (E-B3b, figB_splitting):** razón
Γ(2p₃/₂)/Γ(2p₁/₂) al MISMO s:
| Zα | s=0.24 (cerca del umbral) | s=0.30 (continuo profundo) |
|---|---|---|
| 0.007 | 1.000 | 1.000 |
| 0.146 | 1.064 | 1.012 |
| 0.292 | 1.333 | 1.051 |
| 0.438 | 2.404 | 1.138 |
| 0.584 | **77.5** | (2p disueltos) |
| ≥0.67 | **∞ (2p₁/₂ LIGADO, 2p₃/₂ resonante)** | — |

**Dos descubrimientos físicos:**
1. **La ruptura de simetría de espín de las anchuras se AMPLIFICA cerca
   del apantallamiento crítico**: a s=0.24 escala hasta 77× y diverge
   (compañero ligado desde Z=92 — exactamente donde E-A3 pone
   s_c(2p₁/₂) > 0.24 ✓ consistencia entre proyectos); a s=0.30 el
   mismo doblete apenas rompe (≤1.14). Es un efecto de UMBRAL: la
   anchura es exponencialmente sensible a la posición relativa al
   s_c propio de cada compañero. El desdoblamiento λc hecho dinámica.
2. **Inversión del signo de dE_r entre regímenes**: cerca del umbral
   E_r(3/2) > E_r(1/2); en el continuo profundo se invierte
   (Z=40, s=0.30: 9.23 vs 8.98) — eco relativista de las inversiones
   que Shi et al. reportan entre dobletes distintos.

**Trayectorias (E-B3c, f3_traj.csv):** Z=80, s = 0.245→0.30: el 2p₁/₂
recorre (5.3, Γ=2.6) → (39.1, Γ=87.9) mientras el 2p₃/₂ se disuelve
tras s=0.26 — figura de plano complejo lista.

**Estado del plan B: F0 ✓ F1 ✓ F2 ✓ E-B2c ✓ F3 ✓. Falta F4:
verificación adversarial + manuscrito.**
