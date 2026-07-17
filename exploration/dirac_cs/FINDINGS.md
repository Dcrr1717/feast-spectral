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
