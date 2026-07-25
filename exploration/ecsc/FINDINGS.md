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

## E1b — Reconciliación con Jiao-Ho EPJD 75:313 (2021) (2026-07-24)

Su paper (encontrado durante la redacción): λc=0.720524 (GPS alta
precisión; el nuestro 0.720273 = sesgo de difusividad documentado),
λr+=0.90035 donde nacen "well-defined shape resonances", y el rango
intermedio 0.7205<λ<0.9003 declarado **"unexplained"** (conjeturan
estados virtuales).

**Nuestra respuesta a su pregunta abierta:** la continuación del polo
1s muestra que cruza el umbral y se FUNDE con el continuo rotado.
Diagnóstico con huella nueva: el cociente (Γ/2)/E_r del polo rastreado
converge a tan(2θ) del ángulo de rotación (λ=0.75: 1.54 ≈ tan(1.0);
λ=1.0: 1.55) — la firma inequívoca del rayo del continuo. Cerca de λc
(0.715–0.73) el polo es genuino (q~5e-5, se mueve con λ); más allá de
λ≈0.75 el tracking sigue estados del rayo. CONSISTENTE con su
conjetura de estado virtual: no hay resonancia s aislada en el rango
intermedio resoluble por escalado uniforme.

**Cazas del polo agudo sobre λr+ (2 intentos, ambos negativos con
diagnóstico):** (1) escaneo θ-estable → cuerda λ-independiente =
artefacto; (2) discriminador de dos escalas (β=0.15 vs 0.06) → la
cuerda densa del rayo empareja entre bases por densidad, no por
física (todos los "hits" λ-independientes y sobre el rayo tan(2θ)).
El polo de cola de Jiao (si E_r pequeño con Γ diminuto) requiere
exterior scaling o su mapeo GPS — fuera de alcance del setup uniforme;
lo citamos sin reclamar reproducción. HERRAMIENTA NUEVA para el kit:
el test del cociente tan(2θ) como fingerprint de artefacto de rayo.

## E2 — Mapa (F, λ) + doblete n=2 (2026-07-24) · COMPLETO

**95/95 puntos (5F × 19λ, λ∈[0,0.9]) SIN pérdidas + 26/26 doblete.**
CSV: e2_map.csv, e2_manifold.csv. Peor q: 1e-4 (F=0.10, λ=0.9).

**RESULTADO ESTRELLA (más fuerte que el análogo del Paper D): el
estado FUNDAMENTAL cruza E_r=0 y sobrevive como resonancia de energía
positiva en los CINCO campos**, con cruce en λ*(F): ≈0.686 (F=0.02),
0.681 (0.04), 0.673 (0.06), ~0.66 (0.08, 0.10 — leve corrimiento con
F), continuando suave hasta λ=0.9 (F=0.02: E_r=+0.0493, Γ=0.0966).
Bajo Debye el fundamental JAMÁS cruza en nuestro rango (Paper D:
E_r<0 en todo el mapa). La conexión con E1: a F=0 el polo se funde
con el continuo al cruzar; el campo lo mantiene aislado (la resonancia
Stark es de otro carácter — el campo abre el canal de túnel que la
define) — física nueva limpia.

**Doblete n=2 ECSC (F=0.02):** azul cruza E_r=0 en λ≈0.0936 (vs 0.112
Debye — el ECSC desliga antes, consistente con λc(2p) 0.148 vs 0.220);
roja se acerca al umbral sin cruzar hasta λ=0.3 (E_r=−0.008,
Γ=0.21). Ratio Γ_red/Γ_blue: 1.97 → 2.63 (λ=0.3).

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
