# FINDINGS — Paper E: ECSC + Stark (ecsc)

Contexto: paper #1 de la secuencia del GAP_SCAN_2026H2. Extiende
Wang-Chen-Guo JPB 52:025001 (2019, demo en puntos dispersos) al primer
mapa sistemático (F, λ) bajo ECSC + manifolds + comparación Debye-vs-
ECSC. Sesión autónoma 2026-07-24.

## E0 — Compuertas del motor (2026-07-24) · PASADAS

**Identidad clave (cero código nuevo de operadores):**
V_ECSC = −(Z/r)e^{−λr}cos(λr) = −(Z/2r)[e^{−μ₊r} + e^{−μ₋r}] con
μ± = λ(1±i) — media suma de dos Yukawas con apantallamiento complejo.
Bajo r→re^{iθ} el decaimiento exige cosθ−sinθ>0 → **tope analítico
θ<π/4** (assert en ecsc_core.py; nuestros ángulos 0.18–0.42 caben).

- P-E0.0 hidrógeno λ=0: dev 6.8e-11 ✓ (raise-on-fail).
- P-E0.1 regresión Yukawa vs motor Paper D: max|ΔH| = 2.2e-16 ✓.
- P-E0.2 expansión pequeño-λ EXACTA verificada: e^{−x}cos x =
  1 − x + x³/3 − x⁴/6 (¡SIN término x²!) → E_1s = −1/2 + λ − λ³ + O(λ⁴)
  (⟨r²⟩=3); medido a λ=0.01: −9.88e-7 vs −1e-6 ✓. Contraste Yukawa:
  E_1s = −1/2 + λ − (3/4)λ² — el ECSC liga MÁS a λ chico.
- P-E0.3 curva 1s: λc(1s) ≈ 0.720273 (bisección variacional nb=160;
  con difusividad de umbral como en el paper NR). Deriva entre grillas
  ≤7.5e-11 en todo el rango ligado.

## E1 — Estructura de polos field-free a través de λc (2026-07-24)

**ℓ=1 (el resultado limpio): λc(2p) = 0.148205** (vs Yukawa 0.220217 —
la cola oscilatoria desliga antes). El polo 2p se despega como
**resonancia genuina con θ-estabilidad 1e-13** (¡máquina!) cerca del
umbral: λ=0.153: +0.0023587−0.0002779j (q 1e-13); se ensancha rápido
(λ=0.218: q 3e-9; λ≳0.26: q>5e-3, estados muy anchos requieren más
base — régimen citado con cautela).

**ℓ=0 (honestidad): NO hay resonancia s limpia.** El polo 1s cruza λc
moviéndose suave (dE/dλ≈0.044 estable, q(base)~1e-5) PERO la
trayectoria θ deriva monótonamente sin meseta (drift 9e-4 a λ=0.75,
7e-3 a λ=1.0) — es una **estructura ancha de umbral** (Γ≳2E_r) apenas
aislada: la primera capa repulsiva de la cola (altura ~e^{−π}λ/π ≈
0.014λ) es demasiado baja para atrapar el polo. El contraste físico
con Debye: en Yukawa ℓ=0 ni siquiera existe esa estructura (polo
virtual puro, sin barrera). El candidato de P-E0.4 (E_r≈4e-4 estático
en λ) ERA un artefacto — cazado por el criterio de movimiento
(lección Paper B aplicada).

**Lección de método:** para polos s cerca del umbral con colas
oscilantes haría falta exterior complex scaling o continuación de
Jost — declarado como límite honesto, no perseguido aquí.

## E2 — Mapa (F, λ) + doblete n=2 · corriendo

Tracker del Paper D v2 con ancla ECSC (curva F=0 mientras exista,
extrapolación compleja más allá de λc); mapa λ∈[0,0.9]×F∈{0.02..0.10}
+ doblete n=2 a F=0.02. Resultados al aterrizar.

## E3 — Catálisis a ligadura igualada · pendiente de correr

Anclas E_b ∈ {−0.40,−0.30,−0.20,−0.12}: bisectar λ_Y y λ_E con E_1s
igual a 1e-9, comparar Γ(F=0.04, 0.08) — aísla el efecto de forma de
barrera/cola del efecto trivial "el apantallamiento sube el nivel".

## Pendientes de literatura

- Wang-Chen-Guo 2019 full text (paywalled; fijar frontera de novedad
  exacta y comparar sus puntos si tabulan).
- λc(1s) ECSC: cross-check contra el JPA 57 (2024) de parámetros
  críticos de la familia generalizada (nuestro caso = θ=π/4 de su
  parametrización, β=λ√2) — buscar arXiv.
- Clásico PRA 6:1391 (1972) energías 1s ECSC.
