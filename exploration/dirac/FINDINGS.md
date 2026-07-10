# Diagnóstico DKB radial — veredicto cruzado (2026-07-09)

**Método:** derivación propia + derivador independiente (workflow adversarial,
con verificación numérica propia del agente) + tres discriminadores
empíricos (D1 directo-vs-ensamblado, D2 deriva de cuadratura, D3
condicionamiento de S). Certificaciones de ensamblaje activas en cada build.

## 1. Bugs de código — encontrados y muertos (3)

| Bug | Detector | Estado |
|---|---|---|
| `O₋O₊f` con centrífugo `κ(κ+1)` (es `κ(κ−1)`; el código calculaba `O₊O₋`) | Derivación propia, **confirmado por agente independiente** con prueba numérica (reproduce 1.2e-7 / 1.6e-9 / colapso Z=80; corregido da 2.2e-11) | FIXED |
| `d2q` usando `dq[j]` a medio construir (refactor) | **Certificación C2/C3 automática** en su primer uso (0.958 rel.) | FIXED |
| Tolerancia de certificación 1e-9 vs ruido real 2e-9 | Calibración (bugs estructurales O(0.1–1); ruido <1e-8) | FIXED (1e-7) |

Trampa didáctica del bug 1: `κ(κ+1) = ℓ(ℓ+1)` para AMBOS signos de κ — el
término erróneo "parece" el centrífugo correcto. El orden del producto
importa: `[O₋,O₊] = −2κ/r²`.

## 2. La convención de signo s en Ψ⁻ = (s·(1/2c)O₊φ, φ) — NO es un bug, es física

- Flip del espinor COMPLETO: congruencia, espectro invariante (agente ✓).
- Flip de UNA componente (s=+1 ↔ s=−1): **cambia el subespacio generado**.
  Ambos spans dan pares de Galerkin legítimos (D1 lo confirma), pero:
  - **s=−1** (Shabaev, balance positrónico físico P≈−(1/2c)O₊Q): span
    conforme; limpio a Z=1 en todos los canales (1e-11–1e-12), min(S)=1.00.
  - **s=+1** (anti-balanceado): poluciona — estado espurio mid-gap visible
    desde nb≈60 INCLUSO a Z=1 (el agente a nb=40 no lo veía; nb=60 sí:
    −2.5e4). Marca de polución genuina: aparece al crecer N.

## 3. Los "fallos" restantes a Z=80 NO son bugs — son los fenómenos del Teorema 1

Con ensamblaje certificado fiel (D1: todos los bloques coinciden con la
aplicación directa del operador a 1e-11–1e-15), a Z=80 (ν=0.58, subcrítico):

**(a) κ=−1 (ℓ=0): violación del dominio de forma — DKB+Coulomb puntual.**
La componente grande de Ψ⁻ es O₊φ, que NO se anula en r=0 para ℓ=0 →
⟨O₊φ|Z/r|O₊φ⟩ es log-divergente. Empíricamente (D2b): el elemento crece
−59732 → −68547 → −77364 al duplicar nq — **incrementos de −8815, −8816:
lineal exacto en log(nq)**. El "autovalor" deriva (−1.56 → −1.77 → −1.95).
El DKB con núcleo puntual y ℓ=0 contiene funciones FUERA del dominio de
forma. Consistente con que Shabaev formula DKB para núcleo finito.
→ Para el teorema: DKB-puntual requiere prefactor r^γ adaptado a la cúspide
(o núcleo finito). Acopla (L1) con (L3).

**(b) κ=+1: colapso variacional de partición no conforme — polución pura.**
h~O(r) en el origen (dominio de forma OK), y el estado espurio a
binding≈−22466 (mid-gap, E≈−3687) es **estable frente a cuadratura**
(idéntico a 5 cifras con nq=2000/4000/8000). No es divergencia ni ruido:
es un autovalor convergido del subespacio que no corresponde a ningún
estado físico. El operador de Dirac no está acotado por abajo → la cota
"variacional" en la brecha SOLO vale si la partición conforma (min-max de
DES/Morozov–Müller). Su violación aquí es el detector de no-conformidad.

**(c) κ=−2, +2 (ℓ≥1): esencialmente limpios** (−4e-8, −1.5e-4 moderado) —
la severidad decae con ℓ, como predice el análisis del origen.

**(d) D3:** condicionamiento de S benigno en todos los casos (≤1.3e4) —
descartado como mecanismo.

## 4. Estado del instrumento

- Z=1: precisión de máquina en 4 canales × 2 métodos. RKB variacional
  SIEMPRE (todas las Z probadas) — consistente con que su span está en el
  dominio de forma y conforma mejor.
- Convergencia (B2): Z=1 instantánea; Z=80/120 lenta y algebraica —
  limitada por la cúspide r^γ (γ=0.81/0.48), el fenómeno (L3).
- El laboratorio queda CERTIFICADO: integrales fieles + fenómenos de
  polución reproducibles bajo control. Es exactamente el instrumento que
  el Teorema 1 necesita — produce el fenómeno sin artefactos.

## 5. Mejoras pendientes (del agente derivador, adoptar en próxima sesión)

1. Pencil desplazado H−c²S ensamblado analíticamente (mata el piso 1e-11
   por cancelación de masa en reposo: 2c²·eps_mach).
2. Certificación I1: ME(f,h) = ME(g,f) (IBP; caza errores de signo O±) — AÑADIDA.
3. Chequeo I7: ensamblar H₋₊ por la ruta IBP opuesta y comparar con Hpmᵀ
   ANTES de simetrizar (el reuso de Hpm.T enmascara esa clase de error).
4. Cuadratura graduada en el origen para el canal ℓ=0.
5. Prefactor r^γ (base adaptada a la cúspide) para Coulomb puntual a Z alto.

## 6. Medición δ_N / Lema (L1) — sesión 2026-07-10

**Setup:** defecto de conformidad respecto a la partición espectral libre
Λ± (canal κ=−1, β=0.9, referencia Laguerre μ_ref=0.5 con nb_ref=300,
C1_ref=1e-11; residuos de completitud monitoreados por columna).

**Global (δ_N = ‖(I−Π_N)Λ₊Π_N‖):** ordenamiento robusto en las tres
referencias probadas: DKB s=−1 ≪ naive ≡ RKB < DKB s=+1 — reproduce
exactamente el orden de polución observado. PERO δ_N crece con N en TODOS
los spans (a β fijo): el enunciado crudo "δ_N→0" del roadmap es FALSO tal
cual. El lema debe ser POR MODO (L1').

**(L1') Por columna — amplitudes de fuga de sector, medidas (nb=80):**

| familia | fuga (amplitud) | exponente ajustado |
|---|---|---|
| DKB s=−1 electrónica (f, +(1/2c)O₋f) vs Λ₋ | (p/2c)² | **2.15** |
| dividida (f, 0) vs Λ₋ | (p/2c)¹ | **0.95** |
| DKB s=−1 positrónica (−h, f), masa en Λ₊ | (p/2c)⁴ | **3.60** |
| DKB s=+1 positrónica (+h, f), masa en Λ₊ | (p/2c)² | **1.59** |

Ratio de conformidad s=−1 vs s=+1: 2300× (modos bajos) a 33× (j=80).
Confirma numéricamente el balance positrónico P=−(1/2c)O₊Q (signo menos).
Escala de momento por columna: p_j/2c ~ β·√j/c (0.0033 → 0.034, j=1→80).

**Contenido para el teorema:** (L1') las familias balanceadas conforman a
segundo orden en p/2c, las divididas/anti-balanceadas a primer/orden cero;
δ_N global ~ (β√N/c)^s. La no-polución debe formularse con el defecto por
modo relativo a la brecha (o compresiones à la Lewin–Séré), no con δ_N→0.
Queda por conectar: tasa de defecto × acoplamiento Z → umbral de aparición
de polución (explicaría "limpio a Z=1, poluido a Z=80" cuantitativamente).

## 7. E-A1: validación externa y arbitraje de anclas Dirac-Yukawa (2026-07-10)

**Cadena de referencia:** Krauthauser–Hill (Can. J. Phys. 80:181, 2002;
benchmarks de 40 dígitos, unidades c=1) → Baye–Filippin–Godefroid (PRE
89:043305, 2014; malla de Lagrange) → Chen–Shao–Wu (arXiv:2604.11063,
2026; ALLSM, Tabla 4 en Hartree).

**Convenciones descubiertas (no declaradas en ningún paper de la cadena):**
- c = 137.0359895 (CODATA **1986**) — deducida del patrón de residuos
  (2·Δc/c·|E| = 1.07e-5 exacto) y confirmada independientemente por el
  verificador a 10 dígitos desde ambos casos.
- **ERRATA en CSW Tabla 4:** "Case 2: V0=0.4c, λ=0.07c" no corresponde a
  su propia energía; los dígitos 4↔7 están transpuestos. El caso real
  (Baye Tabla VI) es **V0=0.7c, λ=0.04c**. Consistente en prosa, tabla y
  figura de CSW.

**Resultados (cusp, nb=80, invariante adimensional f=E/c²):**
| Caso | nuestro f − KH(40d) | Baye LM | ALLSM 2026 | veredicto |
|---|---|---|---|---|
| 1 (0.1c, 0.01c) | ~1e-14 | (= ref, 5.4e-13 de KH) | off **4.4e-11** | KH ✓ |
| 2 (0.7c, 0.04c) | **−1.9e-15** | off 3.7e-11 | off **7.6e-9** | KH ✓ |

Ambos arbitrajes favorecen al benchmark clásico de 40 dígitos; nuestro
método de estructura exacta se sienta sobre él a precisión de máquina en
el invariante adimensional. Bonus técnico: a γ=0.714 (caso 2) la
cuadratura Legendre se atasca en 1e-4 mientras la regla cusp exacta da
deriva 3e-11 — demostración en vivo del claim de cuadratura exacta.

**Piezas para el paper:** sección de validación externa completa
(reproducción de 40 dígitos + doble arbitraje + errata detectada en la
competencia + convención oculta identificada). El pencil desplazado
(pendiente) bajaría el margen del arbitraje del Caso 1 de 4× a ~1000×.
