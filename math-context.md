# math-context.md — Teorema 1 (no-polución Dirac radial) · v1.0 (2026-07-11)

Tabla de símbolos, hipótesis globales y convenciones para toda derivación
del track teórico. Toda prueba del proyecto debe referirse a este archivo.
Numeración de resultados: L1', L2, L3 (lemas), T1/T2 (teoremas), según
`TEOREMA_NOPOLUCION_roadmap.md` y `TEOREMA1_enunciado.md` (raíz del
proyecto padre).

## 1. Símbolos

| Símbolo | Significado | Tipo/Dominio |
|---|---|---|
| c | velocidad de la luz (unidades atómicas Hartree) | c > 0; numéricamente 137.035999084 (CODATA 2022) salvo indicación |
| κ | número cuántico relativista | κ ∈ ℤ∖{0} |
| ν | acoplamiento de Coulomb Z/c = Zα | ν ∈ [0, 1) |
| γ | exponente de cúspide √(κ²−ν²) | γ > 0 |
| H | espacio de Hilbert radial L²((0,∞), dr)² (componentes reales) | — |
| O± | operadores radiales ∓d/dr + κ/r | dominio denso ⊂ L²(0,∞) |
| h_κ | operador de Dirac radial [[V+c², cO₊],[cO₋, V−c²]] | V = −Z/r (o apantallado) |
| H₀ | operador libre (V ≡ 0) | — |
| Λ± | proyectores espectrales de H₀ sobre energía ≥ c² / ≤ −c² | ortogonales, Λ₊+Λ₋ = I |
| p | módulo de momento del modo | p > 0 |
| t | parámetro adimensional p/(2c) | t > 0 |
| E(p) | √(c⁴ + c²p²) (rama positiva) | — |
| r(t) | razón exacta de componentes: 2t/(1+√(1+4t²)) | 0 < r(t) < t para t > 0 (se prueba) |
| {f_j} | familia escalar ortonormal (Laguerre generalizada, prefactor r^μ) | j = 1..N |
| K | matriz de acoplamiento discreta K_ij = ⟨f_i, O₊ f_j⟩ | K ∈ ℝ^{N×N}, NO simétrica |
| V_N | subespacio de prueba 2N-dimensional ⊂ H | según familia (naive/RKB/DKB±) |
| Π_N | proyector ortogonal sobre V_N | — |
| δ_N | defecto de conformidad ‖(I−Π_N)Λ₊Π_N‖₂ | — |
| ℓ(u) | fuga de sector de un vector de prueba u (Def. 3 abajo) | ∈ [0,1] |

## 2. Hipótesis globales

- H1. Todos los espacios son reales; producto interno ⟨u,v⟩ = ∫₀^∞ (u₁v₁ + u₂v₂) dr; ‖·‖ inducida. En dimensión finita: euclídeo estándar, ‖·‖₂ para matrices = norma espectral.
- H2. Identidad IBP: ⟨f, O₊ g⟩ = ⟨O₋ f, g⟩ para f, g en el dominio común con términos de frontera nulos en 0 e ∞ (garantizado por prefactor r^μ, μ > 1/2, y decaimiento exponencial; certificación C4 del solver la verifica numéricamente en cada build).
- H3. Modelo discreto de referencia: H₀ᴺ = [[c²I, cK],[cKᵀ, −c²I]] ∈ ℝ^{2N×2N}, simétrico. Λ±ᴺ = proyectores espectrales sobre autovalores positivos/negativos de H₀ᴺ. Se asume K de rango completo salvo indicación (los σ_j = 0 se tratan como caso frontera aparte).
- H4. Los enunciados "por modo" usan la SVD real K = UΣVᵀ, Σ = diag(σ_1..σ_N), σ_j ≥ 0; u_j, v_j columnas de U, V.
- H5. Convención de fibra: e₁ = (a, 0), e₂ = (0, b) con a, b ∈ ℝᴺ unitarios; "coordenadas de fibra (x, y)" significa x·e₁ + y·e₂.

## 3. Definiciones

- Def. 1 (familias de prueba por modo, coordenadas de fibra): balanceada u_bal(t) = (1, t); dividida u_spl = (1, 0); positrónica balanceada u_pos(t) = (−t, 1); positrónica anti-balanceada u_anti(t) = (+t, 1).
- Def. 2 (matriz de fibra): h(p) = [[c², cp],[cp, −c²]] ∈ ℝ^{2×2}.
- Def. 3 (fuga de sector): para u de tipo electrónico, ℓ(u) := ‖Λ₋u‖/‖u‖; para u de tipo positrónico, ℓ(u) := ‖Λ₊u‖/‖u‖ (el proyector del sector "equivocado" en cada caso).

## 4. Convenciones

- Expansiones asintóticas: t → 0⁺; O(t^k) uniforme en compactos de t.
- "Exponente de fuga" = exponente del término dominante de ℓ(u(t)) cuando t → 0⁺.
- Verificación numérica: tolerancia relativa 1e-9 salvo indicación; muestras aleatorias uniformes en los dominios declarados.
- Prohibido introducir símbolos nuevos sin definirlos aquí o en el punto de uso con referencia explícita.

## 5. Estado de resultados

- L1'-fibra: **PROBADO** (2026-07-11, FINDINGS §9, verificación en `exploration/dirac/l1prime_fiber_lemma.py`). Exponentes exactos por modo 3/1/3/1 con constantes 1/1/1/2; identidad clave t − r = 4t³/(1+w)². Sustituye los exponentes AJUSTADOS de FINDINGS §6 (2.15/0.95/3.60/1.59), contaminados por mezcla de modos (columnas Laguerre ≠ modos SVD).
- L2 (defecto global exacto): **PROBADO** (2026-07-11, FINDINGS §10, verificación en `exploration/dirac/l2_global_defect.py`). δ_N(familia única sin truncar) = max_j ℓ_j√(1−ℓ_j²), idéntico para Λ₊ y Λ₋; espacio DKB par completo: δ = 0 exacto. Ley de escala medida: σ_max ≈ 0.64·βN (q = 0.992, builds certificados).
- C-umbral (defecto×Z → aparición de polución): **CONJETURA, sin prueba** — enunciado preciso en FINDINGS §10.
- L3, T1, T2: ver roadmap; pendientes.
