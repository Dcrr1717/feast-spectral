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
