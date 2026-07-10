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

## 8. Pencil desplazado H−c²S + tabla E-A2 (2026-07-10, verificado adversarialmente)

**Instrumento nuevo (src/dirac.py):** `solve_naive(shifted=True)` y
`solve_rkb(shifted=True)` ensamblan el pencil H−c²S analíticamente (la
masa sobrevive solo como −2c² en el bloque inferior) y refinan cada
autovalor de brecha con el cociente de Rayleigh sumado exactamente
(math.fsum). El bloque O(c²) se contrae contra las componentes PEQUEÑAS
del estado electrónico → cada término de la forma cuadrática es
O(|E_binding|).

**Veredicto del refutador (PARTIAL → adoptado):** álgebra de bloques
confirmada (derivación independiente + numérica, ambas rutas de código);
piso real ≈ **1e-14 relativo (50–100·eps)**, NO ~eps como decía el
docstring original (ya corregido): lo fijan el pre-redondeo de los
productos x_i·H_ij·x_j antes del fsum y el redondeo de ensamblaje de las
entradas O(c²). Ganancia real: Z=1 gana 3–4 órdenes (1e-11 →
2e-16..5e-15); en estados profundos (Z=92, |E|~4861) el piso
no-desplazado ya era ~1.7e-15 rel y ambos comparten el piso de
ensamblaje ~1e-14. "Nunca peor que el no-desplazado" vale con holgura
3e-11 (piso backward de LAPACK n·eps·2c², n=120), no con 1e-12. Para
RKB, cond(Sss)~1e6–1e7 infla el piso del eigensolver a ~1e-9 absoluto —
no citar 5e-11 ahí.

**Tabla E-A2 (ea2_point_table.py/.csv):** Z ∈ {92,100,110,120,130,136},
κ = ±1,±2, 3 estados por canal, base γ + cuadratura cusp + pencil,
c = 137.035999084 DECLARADO, referencia Sommerfeld en decimal a 40
dígitos, β elegido por deriva mínima (nunca por cercanía al exacto).
**72 niveles, TODOS entre 1e-16 y 1.3e-13 relativo con nb=45.**
Z=136 (ν=0.9924, γ=0.123, pleno régimen de extensión distinguida):
1s = −16474.1124540939 Ha a 1e-13. CSW puntual llega a Z=92.
Consistencia interna: κ=+1 y κ=−1 dan el nivel n=2 idéntico a 2e-14 rel
(canales independientes). Certificaciones C1/C4 ≤ 6.2e-13 en 144 builds.
**Verificación independiente (CONFIRMED):** protocolo distinto (nb=55,
β ∈ {0.65,0.95}·Z) reproduce todo a ≤5e-13.

**Espurio κ>0 — teorema conocido, medición nueva:** el canal naive κ>0
carga EXACTAMENTE un estado espurio degenerado con el 1s del canal −κ.
Medido: degeneración a 4e-15–8e-13 relativo en los 12 casos (ambos
verificadores). NO es novedad el fenómeno: espurios conocidos desde
Drake–Goldman PRA 23:2093 (1981), clasificación Pestka Phys. Scr. 68:254
(2003), **teorema de degeneración exacta: Tupitsyn–Shabaev, Opt.
Spectrosc. 105:183 (2008)** (hipótesis: misma base escalar para ambas
componentes y ±κ — exactamente nuestro caso, pues μ=γ depende de κ²);
cura DKB en Shabaev PRL 93:130405 (2004). Nuestra contribución:
confirmación al piso de redondeo + conexión con la anatomía L1'.
Citar así, no presentar la degeneración como novedad.

**RKB con base γ es mal instrumento en κ<0 puntual:** rel 5.4e-3 en el
1s de Z=92 κ=−1 (por eso el cross-check da nan ahí) — la prescripción
O₋f no reproduce la estructura de cúspide r^γ; en κ>0 cruza a 1e-9.
Dato de anatomía: el balance cinético restringido y la base adaptada a
la cúspide son prescripciones en CONFLICTO en el núcleo puntual.

**E-A1 con pencil (margen actualizado):** Caso 1: f−KH = **+1.4e-16**
(eps de máquina contra el benchmark de 40 dígitos); en Yukawa domina la
truncación de base (~2e-12 Caso 1, ~2e-10 Caso 2 — la cuadratura cusp no
es exacta para e^{−λr} y roots_genlaguerre limita nq ≲ 500), así que
ambos márgenes de arbitraje quedan en **~25×** — decisivos, honestamente
cuantificados. La proyección "1000×" del §7 era para Coulomb puro; en
apantallados no aplica.

**Límites del instrumento anotados:** scipy roots_genlaguerre desborda
para nq ≳ 500 (alternativa futura: Golub–Welsch propio); el refinamiento
RQ está verificado hasta brechas de 1e-4 Ha entre niveles y sin probar
por debajo.
