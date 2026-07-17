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

## 9. Lema L1'-fibra PROBADO: exponentes exactos 3/1/3/1 (2026-07-11)

**Corrección mayor a §6.** Los exponentes ajustados de §6
(2.15/0.95/3.60/1.59) eran artefactos de mezcla de modos: las columnas
Laguerre NO son modos del operador libre, y el proxy p_j ~ β√j era
heurístico. El enunciado correcto es exacto y algebraico
(síbolos/hipótesis en `math-context.md` v1.0; verificación completa en
`l1prime_fiber_lemma.py`; protocolo math-rigor).

**Lema L1'-fibra.** Sea H₀ᴺ = [[c²I, cK],[cKᵀ, −c²I]] con
K_ij = ⟨f_i, O₊f_j⟩ y SVD K = UΣVᵀ. Entonces:
(A) ℝ²ᴺ = ⊕_j span{(u_j,0),(0,v_j)}, cada plano es H₀ᴺ-invariante con
matriz de fibra h(σ_j) = [[c², cσ_j],[cσ_j, −c²]], y los proyectores
Λ±ᴺ se descomponen por bloques (cálculo funcional de matrices
block-diagonales).
(B) Con t = σ/2c, w = √(1+4t²), r = 2t/(1+w) (razón exacta de
componentes; r < t pues t−r = t(w−1)/(w+1) > 0), las fugas de sector
por modo son EXACTAS:
- dividida (1,0):        ℓ = r/√(1+r²)                    = t − (3/2)t³ + O(t⁵)
- balanceada (1,t):      ℓ = 4t³/[(1+w)²√(1+r²)√(1+t²)]   = t³ + O(t⁵)
- positrónica (−t,1):    ℓ = idéntica a la balanceada      = t³ + O(t⁵)
- anti-balanceada (t,1): ℓ = (t+r)/[√(1+r²)√(1+t²)]       = 2t − 3t³ + O(t⁵)
(La identidad clave t − r = 4t³/(1+w)² es algebraica, no asintótica.)

**Verificación (5 capas):** simbólica SymPy (todas las identidades → 0);
numérica 2×2 aleatoria, 3 escalas de c + fronteras t∈{1e-8, 0.5, 5}
(1.1e-16); reducción discreta con K aleatoria n∈{2,3,5}×20 (7.2e-16);
**sistema real** (Laguerre κ=−1, β=0.9, nb=240): |medido − fórmula| =
5.0e-16, exponentes ajustados 0.9999/2.9999/2.9999/0.9999, constantes
0.9996/0.9992/0.9992/1.9991. La reducción es EXACTA en el discreto: no
hay límite continuo que tomar.

**Consecuencias para el teorema:**
1. Las familias balanceadas conforman a TERCER orden por modo (no
   segundo): mejor de lo que creíamos. La razón de conformidad
   anti/balanceada por modo es 2/t² + O(1).
2. El contraste correcto es 3-vs-1 (balanceada vs dividida/anti), con
   constantes exactas 1 y 2 — enunciado limpio para (L1') del roadmap.
3. δ_N global por familia = max sobre modos de las fugas → tasas
   (β√N/c)³ vs (β√N/c)¹ a β fijo; el crecimiento con N observado en §6
   queda explicado por t_max = σ_max/2c creciente.
4. El "exponente 3.60" de la positrónica s=−1 era el t³ verdadero
   contaminado; el "1.59" de s=+1 era el 2t. El ordenamiento de §6
   (s=−1 ≪ naive≡RKB < s=+1) sigue siendo cierto y ahora tiene
   constantes exactas.

## 10. Lema L2 PROBADO: defecto global exacto y ley de escala (2026-07-11)

**Lema L2** (verificación: `l2_global_defect.py`, protocolo math-rigor).
Para una familia de prueba de UN vector por fibra SVD (prescripción sin
truncar), (I−Π_N)Λ±Π_N es block-diagonal por fibras y

    δ_N = max_j ℓ_j·√(1−ℓ_j²)     (EXACTO, mismo valor para Λ₊ y Λ₋)

con ℓ_j las fugas de L1'-fibra (§9). Prueba por fibra en 3 líneas,
coordinate-free: w = (I−bbᵀ)Λ₊b = g(v₊−gb) con g=⟨v₊,b⟩ ⇒ ‖w‖² =
g²(1−g²); la simetría Λ₊/Λ₋ es ⟨v₊,b⟩²+⟨v₋,b⟩²=1.

**Corolario L2c (conformidad exacta del DKB):** el espacio DKB par
completo (ambas familias, sin truncar) contiene ambos vectores de cada
fibra ⇒ Π restringe a la identidad en cada fibra ⇒ **δ = 0 exacto**.
Medido en el sistema real: 3.8e-15. El defecto no nulo del DKB s=−1 en
§6 es enteramente un efecto de TRUNCACIÓN (nb < nb_ref, mezcla de
modos), no de la prescripción.

**Verificación:** simbólica (cadena escalar 0 + parametrizada 0);
K aleatoria n∈{2,3,5}×20, ambos proyectores (1.3e-15); sistema real
κ=−1, β=0.9, nb=240: spl 0.3549, bal 0.0718, pos 0.0718, anti 0.4943 —
|medido − fórmula| ≤ 3.3e-16 en los cuatro.

**Ley de escala (V4, builds certificados):** σ_max(N) ≈ 0.64·βN con
exponente ajustado q = 0.992 (N = 60..360). El proxy p_j ~ β√j de §6
estaba mal también en la escala: t_max = σ_max/2c crece LINEAL en N.
Consecuencia: δ_N^bal ~ (0.32βN/c)³ y δ_N^spl ~ 0.32βN/c a β fijo —
el crecimiento de §6 queda explicado con constantes. Nota de
instrumento: N=480 con cuadratura Legendre por defecto FALLA la
certificación C1 (0.31) — el sobre de validez termina antes; la
certificación lo detectó sola.

**Conjetura C-umbral (SIN PRUEBA — próximo objetivo del track):**
para el Galerkin de h = H₀ + V sobre una familia truncada con defectos
por modo {d_j}, la aparición de un impostor mid-gap ocurre cuando el
acoplamiento efectivo del sector peligroso (modos con t_j = O(1))
a través de V excede su distancia al borde de la brecha; en unidades
ν = Z/c, el umbral empírico observado (§2-3: limpio a Z=1, impostor
DKB s=+1 desde nb≈60; κ=+1 naive impostor anclado al 1s del compañero)
debe salir de max_j d_j·|⟨V⟩_j|/dist_j ~ 1. Formalizar con
compresiones à la Lewin–Séré es exactamente el paso T1 del roadmap.

## 11. E-A3: mapa relativista de apantallamiento crítico λc(Z,κ) (2026-07-11)

**Instrumento (ea3_critical_screening.py):** V = −Z e^{−λr}/r; base γ +
cusp + pencil desplazado (imprescindible: cerca del umbral binding→0);
para κ>0 se usa ev[1] (ev[0] = espurio de §8, inerte al apantallamiento).
Observable adimensional s_c = λc/Z. Extrapolación E_b = −A(λc−λ)^p con
p LIBRE en la ventana cercana al umbral (auto-bracketing).

**Ley de umbral (corregida tras verificación adversarial):** ℓ=0 (κ=−1)
→ p = 2.00 exacto (ley estándar de pozo corto, confirmada a todo Z);
ℓ=1 (2p) → **p = 1 estándar** (ley lineal de corto alcance ℓ≥1). El
"exponente anómalo 1.25–1.47" reportado antes era ARTEFACTO de ajuste de
potencia libre en ventana ancha (pesaba puntos profundos); el
verificador lo demostró: p→1 al acercarse al umbral. Con p=1 fijo (por
ℓ): κ=−1 s_c(Z=1) = 1.19050 (dev 1.1e-4 de 1.190612); 2p s_c(Z=1) =
0.2156 (dev 4.6e-3 de 0.220216, del Valle–Nader — subextrapolación
cerca del umbral, ver nota). El 1s es preciso; los valores ABSOLUTOS del
2p cargan ~2% de incertidumbre de extrapolación cerca del umbral (base
difusa vs artefacto de sobre-ligado a β<0.02Z). El verificador con
método más fino (regresión log-log, más cerca del umbral) obtuvo 2p(Z=1)
= 0.2203 (dev 1e-4). El desdoblamiento (relativo) es robusto.

**Resultado principal — realce relativista del apantallamiento crítico
(s_c = λc/Z, p fijo por ℓ):**
| Zα | 1s (κ=−1) | 2p₁/₂ (κ=+1) | 2p₃/₂ (κ=−2) |
|---|---|---|---|
| 0.007 | 1.1905 | 0.2156 | 0.2156 |
| 0.657 | 1.2775 | 0.2395 | 0.2177 |
| 0.949 | 1.4783 (+24%) | 0.3076 (+43%) | 0.2204 (+2.2%) |

**Resultado insignia secundario — desdoblamiento de estructura fina del
apantallamiento crítico (CONFIRMADO por verificación independiente,
concuerda a ~4%):** 2p₁/₂ y 2p₃/₂ arrancan degenerados (no-rel) y se
SEPARAN con Zα: a Zα=0.949 el 2p₁/₂ requiere ~40% más apantallamiento
que el 2p₃/₂ para desligarse (0.308 vs 0.220; separación absoluta
≈0.087, verificador ≈0.083). Es un efecto puramente relativista, sin
análogo no-relativista, y nadie lo ha reportado (CSW da un solo punto
crítico, esencialmente no-rel). El 2p₁/₂ (j=1/2, más ligado
relativistamente) es más resistente al apantallamiento — consistente con
su corrimiento Dirac más profundo. El desdoblamiento NO es artefacto de
base (Eb idéntico a 5 cifras en β∈{0.1..1.0}·Z) ni de la elección ev[1]
(en Z=1 ev[1] da el 2p físico, no el 1s espurio).

**Nota (verificación E-A3, veredicto CONFIRMED):** dos correcciones de
claims secundarios adoptadas arriba — (i) exponente ℓ=1 es p=1 estándar,
no anómalo; (ii) el espurio κ>0 NO es "inerte al apantallamiento" (se
mueve cientos–miles de Ha bajo screening), solo degenerado con el 1s de
−κ en λ=0; permanece muchísimo más profundo que el umbral 2p, así que
usar ev[1] es operativamente correcto. Ambos resultados insignia (subida
1s, desdoblamiento 2p) sobreviven robustos.

## 12. TEOREMA T1-acotado PROBADO: dicotomía de polución para V acotado (2026-07-12)

**La conjetura C-umbral de §10 queda PROBADA para potenciales acotados**
(protocolo math_daniel; método declarado: Courant–Fischer/Weyl +
compresión no expansiva + block-diagonalización de L2; verificación
completa en `t1_bounded_window.py`; Defs. 4-5 y H6 en math-context v1.1).

**Teorema T1-acotado.** Bajo H1, H3–H6 (V_d simétrica, ‖V_d‖ ≤ v,
W = diag(V_d,V_d)), con G_v = (−c²+v, c²−v) y las familias de Def. 1
sin truncar:

(a) *Clasificación exacta de Rayleigh por fibra* (identidades simbólicas,
residuo 0): R_bal − c² = 2c²t²/(1+t²) ≥ 0; R_pos + c² = −2c²t²/(1+t²)
≤ 0; R_spl = c²; **R_anti = c²(5t²−1)/(1+t²) está ESTRICTAMENTE dentro
de la brecha para todo 0 < t < 1/√2** (c²−R_anti = c²(2−4t²)/(1+t²),
R_anti+c² = 6c²t²/(1+t²)) **y cruza el centro exacto en t = 1/√5.**

(b) *El padre no tiene espectro en G_v*: |eig_k(H₀ᴺ+W) − (±E(σ_k))| ≤ v
(Weyl) y |±E(σ)| ≥ c² ⇒ spec(H₀ᴺ+W) ∩ G_v = ∅.

(c) *Ventana protegida (no-polución)*: el Galerkin de la familia
balanceada tiene TODOS sus autovalores ≥ c²−v; la positrónica ≤ −c²+v;
la dividida en [c²−v, c²+v]. Ninguna produce autovalores en G_v — para
CUALQUIER potencial acotado, sin condición de pequeñez.

(d) *Polución garantizada y localizada*: los autovalores ordenados del
Galerkin anti-balanceado están a distancia ≤ v de los R_anti(t_j)
ordenados; cada fibra con |R_anti(t_j)| < c²−2v FUERZA un autovalor
dentro de G_v, donde por (b) no hay espectro del padre: espurio genuino,
con posición predicha por la curva R_anti.

**Prueba (5 pasos):** (1) Rayleigh por expansión directa de ⟨u,h(2ct)u⟩
por familia. (2) Compresión no expansiva: ‖QᵀWQ‖ ≤ ‖W‖ = ‖V_d‖ ≤ v
(Cauchy–Schwarz con columnas ortonormales). (3) QᵀH₀ᴺQ = diag(R(t_j))
por la fibra-invarianza de L2 paso 1. (4) Weyl (|eig_k(D+P) − eig_k(D)|
≤ ‖P‖, consecuencia de Courant–Fischer) aplicado dos veces: a D=QᵀH₀ᴺQ
para (c)-(d) y a D=H₀ᴺ para (b). (5) |R|<c²−2v ⇒ |eig| ≤ |R|+v < c²−v
⇒ eig ∈ G_v; con (b), espurio. □

**Verificación (3 capas):** simbólica (todas las identidades de (a) → 0,
R_anti(1/√5)=0 exacto); aleatoria n∈{2,3,5}×20 con V_d de norma
controlada: todos los márgenes ≥ 0 (peor caso −8e-11 = redondeo dentro
de la tolerancia declarada 1e-9); **sistema real** (Laguerre κ=−1,
nb=240, pozo acotado V₀=500 Ha): ‖V_d‖ = 492.75 ≤ 500 ✓ (bound de
compresión verificado), márgenes +20/+102/+5.6/+493, y el teorema
predice y encuentra exactamente **5 impostores garantizados** (fibras
con |R_anti| < c²−2v; el más profundo en R = +4031 Ha, fibra
t = 0.504).

**Alcance honesto:** T1-acotado NO cubre Coulomb (V no acotado, H6
falla). Los fenómenos observados a Coulomb puntual (espurio profundo
naive κ>0 §8; polución DKB s=+1 §2) quedan para **T1-Coulomb
(ABIERTO)**: requiere Weyl relativo a formas o compresiones à la
Lewin–Séré. El mecanismo direccional, sin embargo, ya está identificado
y probado: es la dirección anti-balanceada barriendo la brecha
(R_anti(t) recorre (−c², c²) cuando t recorre (0, 1/√2)).

**Para el paper del teorema:** T1-acotado + L1' + L2 forman un bloque
autocontenido y completamente verificado; T1-Coulomb es el capítulo
final que falta. Para el Paper A (PRA), la Sec. V puede citar (c)-(d)
como "the balanced prescription is pollution-free in the
potential-free window for every bounded potential, while the
anti-balanced direction sweeps the gap — proven in the companion
paper".

## 13. T1-Coulomb, asalto 1: la ruta discreta queda CERRADA con diagnóstico (2026-07-13)

**Herramientas probadas (t1c_coulomb_step.py, protocolo math_daniel):**

**(P1) Ventana de inclusión por resolvente.** a_N := ‖W(H₀ᴺ)⁻¹‖. Para
E real con |E| < c²(1−a_N), E ∉ spec(H₀ᴺ+W). Prueba: factorización
H₀+W−E = (I + W(H₀−E)⁻¹)(H₀−E) + Neumann, con ‖H₀(H₀−E)⁻¹‖ =
c²/(c²−|E|) por cálculo espectral. Verificada (0 autovalores en ventana,
60 casos aleatorios).

**(P2) Protección condicional de la familia balanceada.** α_N :=
‖|H₀|^{−1/2}W|H₀|^{−1/2}‖ y g(t) := R_bal/E = (1+3t²)/((1+t²)√(1+4t²)).
Entonces min eig(Galerkin_bal) ≥ min_j E_j(g(t_j) − α_N); en particular
α_N < g(t_max) ⇒ ningún autovalor balanceado en o bajo el centro de la
brecha. Clave estructural: |H₀ᴺ| = ⊕E(σ_j)·I es ESCALAR por fibra
(h² = E²I), así que la forma ⟨·,|H₀|·⟩ diagonaliza exactamente donde
vive la familia. Verificada (margen ≥ 0 en 60 casos + sistema real).

**MEDICIÓN a Coulomb (base γ + cusp, κ=−1, integrales exactas,
Z ∈ {40,92,118,130}, N ∈ {40,80,160}): las hipótesis DIVERGEN.**
- a_N ~ N^{0.46} (≈√N): Z=40: 0.41→0.56→0.77; Z=130: 3.3→4.6→6.5.
  La ventana P1 es vacía ya en N=40 para Z ≥ 92 (a_N > 1).
- α_N ~ log N (incrementos constantes por duplicación: ~0.06 a Z=40,
  ~0.18 a Z=130): 1.60→1.78→1.96 a Z=130 — divergencia LOGARÍTMICA.
- g(t_max) ~ 3c/(0.64βN) decae como 1/N ≫ más rápido que log N crece:
  la condición α_N < g(t_max) falla irreparablemente (solo Z=40, N=40
  la cumple). El Galerkin balanceado de una sola familia a Coulomb tiene
  autovalores profundamente negativos (reales, no espurios del solver
  físico: el solver usa el espacio COMPLETO 2N donde esas direcciones se
  resuelven contra las positrónicas).

**DIAGNÓSTICO (el hallazgo del asalto):** la obstrucción es estructural
y precisa — **(ΠH₀Π)⁻¹ ≠ ΠH₀⁻¹Π**. La desigualdad de Kato/Hardy que
controla Z/r en el continuo usa el crecimiento ilimitado de |p|; el
operador libre DISCRETO tiene el momento tapado en σ_max ≈ 0.64βN
(ley de escala de L2) mientras la base γ sigue resolviendo la
singularidad 1/r al crecer N. Ninguna cota relativa respecto al libre
discreto puede ser uniforme en N. **Conclusión: T1-Coulomb NO puede
probarse con cotas relativas discretas; las compresiones al nivel del
operador continuo (Lewin–Séré) son NECESARIAS, no opcionales.** Esto
convierte la elección de framework del roadmap en un teorema de
imposibilidad de la alternativa.

**Ruta v2 (plan preciso):** (i) formular la protección con Λ± del
operador libre CONTINUO y el defecto por modo de L1' (las fugas t³/t
son exactamente las compresiones malas cuantificadas); (ii) para el
control de Z/r usar la desigualdad de Hardy relativista en el sentido de
FORMAS respecto al |D₀| continuo (constantes sharp por canal: el techo
Talman ν≈0.906 y el régimen distinguido ν<1 de DES/Esteban–Lewin–Séré);
(iii) la compatibilidad de dominio de forma de la base γ (r^γ exacto) es
lo que permite (ii) hasta ν→1. Este es el capítulo final del paper del
teorema — el track funcional de César.

## 14. Lema L3 PROBADO: representación exacta y dicotomía de tasa (2026-07-13)

**El lema de la cúspide — la teoría de aproximación detrás de todo el
paper (protocolo math_daniel; verificación en `l3_cusp_rate.py`).**

**(L3a) Representación exacta [PROBADO].** Para κ < 0, n_r = 0, ambas
componentes del estado base Dirac-Coulomb son proporcionales al MISMO
perfil r^γ e^{−λr} con λ = Z/|κ|, E = γc²/|κ| (sustitución directa:
las dos filas de la ecuación dan condiciones lineales en (a,b)
CONSISTENTES — verificado simbólicamente para κ=−1,−2). Con μ=γ y
β=Z/|κ| el estado vive en el espacio de prueba de nb=1: **el solver
reproduce el estado base exacto con DOS grados de libertad** — medido a
3e-16..7e-15 relativo para Z ∈ {40,92,130}, κ ∈ {−1,−2}. (Z=130 con
2 DOF: −12838.9194966482.) Demo definitiva de estructura exacta.

**(L3b) Dicotomía de tasa con exponente UNIVERSAL [PROBADO].** Para el
perfil de cúspide u = x^γ e^{−x/2} expandido en la familia ortonormal
q_j^{(μ)}, los coeficientes tienen forma cerrada EXACTA (identidad
∫₀^∞ e^{−x}x^s L_j^{(α)}dx = Γ(s+1)Γ(j+α−s)/(j!Γ(α−s))):

    c_j = N_j · Γ(μ+γ+1)·Γ(j+μ−γ) / (j!·Γ(μ−γ)),

y por Stirling |c_j| ~ [Γ(μ+γ+1)/|Γ(μ−γ)|]·j^{−(1+γ)}. Consecuencias:
1. **El exponente de decaimiento −(1+γ) NO depende de μ**: un prefactor
   desajustado — por cercano que esté — solo mejora la CONSTANTE
   (∝ 1/|Γ(μ−γ)|), nunca la tasa. "Estructura exacta o tasa algebraica
   fija; no hay término medio."
2. Anulación idéntica ⇔ μ−γ ∈ {0,−1,−2,...} (exactitud en μ=γ).
3. Cola L² tras N términos ~ N^{−(1/2+γ)}.
Verificación: fórmula exacta a 1e-9 (nivel de cuadratura) en
μ ∈ {0.9, 1, 2}; razón de constantes μ=2/μ=1 medida 10.1660 vs fórmula
10.1660; μ=γ anula todo (<2.3e-12); pendientes de decaimiento
convergiendo monótonamente al exponente universal (−1.735/−1.724 en
j∈[200,390) → −1.7411). **Este es el teorema que explica por qué
ALLSM (aproximación adaptativa de γ) no puede alcanzar convergencia
exponencial: cualquier μ≠γ da la misma tasa algebraica.**

**(L3c) Tasa de autovalores [MEDIDO; ley CONJETURADA].** Base entera
(μ=1) a Z=80, κ=−1: error relativo del autovalor ~ N^{−1.605} (ratios
por duplicación 3.01/3.05/3.06 — power law limpio). La conjetura
ingenua N^{−(2γ−1)} = N^{−0.62} queda REFUTADA por la medición; la ley
consistente es **N^{−2γ}** = N^{−1.624} (tasa dual
‖(I−Π)ψ‖·‖(I−Π)Hψ‖ ~ N^{−(1/2+γ)}·N^{−(γ−1/2)}). Probar L3c con la
estimación dual de autovalores para operadores de primer orden es el
remate pendiente del lema.

**Tablero del teorema tras L3: L1' ✓, L2 ✓, T1-acotado ✓, L3a ✓,
L3b ✓; abiertos: L3c (ley N^{−2γ}, medida) y T1-Coulomb (ruta v2
delineada en §13).** El paper del teorema tiene su espina dorsal
completa.

## 15. Campaña experimental del teorema: E-T1/E-T2/E-T3 (2026-07-13)

**E-T1 (et1_l3c_universality.py) — universalidad de L3c CONFIRMADA.**
Exponente del error de autovalor vs 2γ_κ predicho: Z=40 (−1.901 vs
−1.913), Z=80 (−1.605 vs −1.624), Z=120 (−0.928 vs −0.966); y con
**μ fraccionario γ+0.3: −1.593** — mismo exponente que μ entero, solo
cambia la constante, exactamente como manda L3b. κ=−2 consistente en
ventana preasintótica (−3.43 vs −3.83). Mecanismo (corregido en vivo:
el producto dual ingenuo es degenerado para autofunción exacta pues
hψ=Eψ): **dE_N ≈ (0.43–0.55)·⟨res,|H₀|res⟩** — cuadrático en la norma
de forma del residuo de proyección, ratio O(1) estable (ventana
N ≤ nb_ref/4; el punto N=160 sube a 0.92 por truncación de la
referencia). L3c queda: ley medida en 5 configuraciones + mecanismo
verificado; falta solo el enunciado/prueba del estimador cuadrático.

**E-T2 (et2_continuum_form_bound.py) — la premisa de la ruta v2
VALIDADA; la adjudicación de constantes, honestamente ambigua.**
α̃(nb; nb_ref) = ‖M^{−1/2}WM^{−1/2}‖ con M = compresión del |H₀| de
referencia grande (proxy del continuo):
- Q1 ✓: α̃ CONVERGE en nb_ref (0.8633→0.8687, incrementos a la mitad
  cada paso) — vs la divergencia log del α discreto (§13). La
  obstrucción ERA el |H₀| discreto; el control continuo es viable.
- Q3: α̃ sube lentamente con nb (0.828→0.896, nb=20→80) hacia el sup
  continuo del canal.
- Q2: α̃/ν crece con ν (0.99 → 1.375): la cantidad medida (restringida
  al span de la base γ, que cambia con Z) queda ENTRE Talman sharp
  (1.104) y Kato genérico (π/2 = 1.571); los datos se inclinan a "Kato
  suavizado por el término de masa", NO a la hipótesis Talman-sharp
  pre-registrada. α̃ < 1 hasta ν ≈ 0.75 a nb práctico → ahí es
  formulable un teorema condicional de protección continua; más allá se
  necesitan las formas descompuestas de DES (consistente con la
  literatura).

**E-T3 (et3_splus_onset.py) — el mapa del impostor s=+1 CORRIGE §2.**
Con detector fino (estado en brecha a rel > 1e-3 de TODOS los niveles
exactos): **no hay umbral en nb — el impostor anti-balanceado existe
desde nb=20 a todo Z ∈ {1..80}**, nace pegado al borde INFERIOR de la
brecha (Z=1, nb=20: binding −36905 ↔ E ≈ −0.965c²) y BUCEA hacia el
centro al crecer nb (−36905 → −35081 → −32173 → −28314 para nb =
20/40/60/80). El "onset en nb≈60" de §2 era artefacto del detector
(inspección de medio-gap); registro corregido. La trayectoria es la
firma del mecanismo probado en T1-acotado: la dirección anti-balanceada
barre R_anti(t) = c²(5t²−1)/(1+t²) desde −c² hacia arriba conforme
crece t_max. Cuantitativamente t_eff(nb) > t_max(nb) a Coulomb (el
potencial no acotado acelera el barrido) — ESA brecha cuantitativa es
exactamente lo que T1-Coulomb v2 debe predecir. Regresión del refactor
solve_dkb(s): s=−1 en Z=1, nb=60 da 6.5e-11 ✓ (récord histórico).

**Mejoras de instrumento en esta campaña (src/dirac.py):**
1. solve_dkb expone s = ±1 (default físico −1; s=+1 solo experimentos).
2. Guardia ruidosa en la regla cusp: scipy roots_genlaguerre produce
   NaN para nq ≳ 350–400, y los pesos wk underflowean peso GENUINO de
   cola para x > 709 → **techo real de la regla cusp: nb ≲ 170**
   (C1 lo caza: 0.307 a nb=240). Corrige el "≳500" de §8. Superarlo
   requiere Golub–Welsch propio en log-espacio (mejora futura).
3. Certificación NaN-proof: `not (max ≤ tol)` en vez de `max > tol`
   (un NaN antes pasaba silencioso — hueco cerrado).

## 16. Golub–Welsch log-espacio + balance atómico DKB (2026-07-16)

**(a) EL TECHO nb≲170 CAYÓ (cusp_rule + laguerre_funcs_scaled).**
Nodos del Jacobi tridiagonal (eigvalsh_tridiagonal, sin NaN a ningún
nq); pesos por la identidad de Christoffel W = 1/Σφ_k² con recurrencia
de exponente corrido (nada under/overflowea jamás); y la evaluación de
la base con **√W fusionado en las funciones** (todos los valores O(1) en
todo nodo; las relaciones de derivadas son lineales con coeficientes
por-nodo independientes de j, así que el factor pasa limpio). Nodos con
peso verdadero > e^600 descartados exactamente (contribución
< e^{−4(nq−nb)} con el margen estándar; C1/C4 certifican cada build).
**Validación: nb=45 regresión intacta (4.6e-16); nb=400 → C1=1.2e-13,
1s a 3.1e-14; nb=800 → 1s a 1.2e-14.** El límite práctico ahora es el
eigh denso, no la cuadratura. Regresiones test_shifted_pencil y
l1prime sin cambios.

**(b) Balance atómico DKB puntual (atomic_balance_dkb.py) — el
pendiente de §3a CERRADO.** Factor regularizado w(r) = c/(2c²+Z/r)
(→ cr/Z en 0, → 1/2c en ∞) aplicado a ambas direcciones del par DKB:
- A1: dominio de forma RESTAURADO — ⟨wO₊f|Z/r|wO₊f⟩ estacionario a
  4e-9 frente a nq (el libre diverge log-lineal, drift −6.6e8 por
  duplicación — la huella exacta de §3a).
- A2: espectros limpios y variacionales (sin violación) a Z=80/120 en
  base entera; el error residual (+6.5e-4 / +1.6e-2) es la truncación
  de cúspide de la base entera — la ley N^{−2γ} de L3c — no del balance.
- A3: autovalor atómico estacionario (−3529.905714 idéntico en
  nq=2000/4000/8000) mientras el libre deriva (−3532.93 → −3533.15,
  no convergente; pasa "cerca" del exacto −3532.19 por accidente y POR
  DEBAJO — no variacional).
Moral: a núcleo puntual con base entera, la cura correcta del dual es
la regularización de Coulomb; con la base γ el problema ni aparece
(L3a). Ambas rutas quedan implementadas y medidas.

**Interpretación física:** el estado ligado relativista es más profundo
(corrimiento de Dirac), así que hace falta más apantallamiento para
empujarlo al continuo; el efecto crece monótonamente con Zα y es máximo
para el estado más ligado de cada j. Es la contraparte "de umbral" del
mapa espectral E(Z,κ) de E-A2.

## 17. E-C1: las constantes del lema puente — la ruta del teorema decidida (2026-07-17)

**Medición (ec1_bridge_constants.py; κ=−1, Coulomb, base γ, referencias
hasta 480 gracias al Golub–Welsch log-espacio):** tres candidatos a
operador de control en el espacio pequeño: A₁ = |ΠH₀Π| (fibra-diagonal,
la ruta P2), A₂ = (ΠH₀²Π)^{1/2} (intermediario de convexidad),
A₃ = Π|H₀|Π (compatible con Hardy).

| Z | nb/ref | α₁ | α₂ | α₃ | c_bridge |
|---|---|---|---|---|---|
| 92 | 40/160 | 1.0005 | 0.8276 | 0.8687 | 3.405 |
| 92 | 80/320 | 1.1298 | 0.8662 | 0.9003 | 4.069 |
| 92 | 120/480 | 1.2075 | 0.8847 | 0.9153 | 4.495 |

**Veredicto E-C1 (lo que César necesitaba para decidir):**
1. **α₂ y α₃ CONVERGEN** (incrementos a la mitad; límites ~0.90 y ~0.93
   a ν=0.67, ambos < 1) mientras α₁ diverge — las DOS rutas continuas
   son viables en constantes; la convexidad (A₂) es la mejor.
2. Órdenes de operador verificados: A₁ ≤ A₂ y A₃ ≤ A₂ (convexidad de
   x² + monotonía de √, confirmadas a redondeo).
3. **El hallazgo estructural: el puente es GRATIS en la diagonal de
   fibras** — ⟨b_j, A₃ b_j⟩/E_j(A₁) = 1.000 en todo el rango físico de
   t (y 0.854 ≤ 1 en t≈29, la dirección buena). El crecimiento de
   c_bridge (3.4→4.5) vive ÍNTEGRAMENTE en los términos CRUZADOS entre
   fibras.

**Ruta para E-C2/T1-Coulomb (concreta y medida):** descomponer
A₃ = D₃ + O₃ en la base de fibras; D₃ ≤ E_j·(1+ε) con ε→0 (medido);
el teorema se reduce a controlar ‖O₃‖ (o su acción sobre la familia
balanceada) — un objeto acotable con la estructura SUSY del continuo.
Alternativa: trabajar con A₂ directamente (mejor constante) explotando
que ΠH₀²Π es computable por canales (−d²/dr² + centrífugo + c⁴).

## 18. E-C2: lemas puente + la ventana partida — T1-Coulomb estructurado (2026-07-17)

**Lemas estructurales (ec2_bridge_lemma.py):**
- **C2a (reducción por componentes) PROBADO:** H_ref² = diag(c⁴+c²KKᵀ,
  c⁴+c²KᵀK) por la estructura SUSY ⇒ |H_ref|, A₂ y A₃ son diagonales
  por componentes (verificado: off-blocks ~1e-15 relativo); A₁ también
  (via (PHP)²). El problema puente de Dirac se reduce a DOS problemas
  escalares tipo Schrödinger.
- **C2b (la cadena tiene nombres) PROBADO por teoremas clásicos:**
  A₁ ≤ A₂ es Kadison–Schwarz ((PHP)² ≤ PH²P) + monotonía de Löwner de
  √; A₃ ≤ A₂ es la desigualdad de Hansen (concavidad de operador de √).
  Citables, verificados numéricamente en E-C1.

**REFUTACIÓN que afila (P-C2.2):** la ruta intrínseca ω_N =
‖E^{-1/2}W_bal E^{-1/2}‖ (sin referencia) CRECE (Z=92: 0.74→0.98;
>1 en ν=0.75) — no solo |ΠH₀Π| falla (§13): CUALQUIER escala de
energía del espacio pequeño falla. El continuo es obligatorio, ahora
con dos evidencias independientes.

**C2e (localización):** el costo del puente Δ = A₃−A₁ escalado por
√(E_iE_j) se concentra en fibras de t GIGANTE (t=19.7, 7.8 — energías
≫ c², lejísimo de la brecha); los modos físicos t ≲ 1 están limpios.

**EL RESULTADO — la ventana partida con ley medida:** restringiendo al
sector t_j ≤ 1:

| Z | ν | ω_low (nb=80/120) | ω_low/ν | η_low |
|---|---|---|---|---|
| 40 | 0.292 | 0.162 / 0.170 | 0.575 | 0.15-0.21 |
| 92 | 0.671 | 0.3859 / 0.3827 | **0.575** | 0.21-0.23 |
| 118 | 0.861 | 0.4915 / 0.4947 | **0.571** | 0.21-0.24 |

**ω_low ≈ 0.57·ν, ESTABLE en nb y CONSTANTE en Z ⇒ ω_low < 1 para
TODO ν < 1** — el sector físico está controlado fibra a fibra hasta el
borde subcrítico. η_low ~ 0.2 (manso).

**Estructura final de T1-Coulomb (para el track de César):**
1. Partir en t* = 1: sector bajo (físico) + sector alto (|E| ≥ √5c²).
2. Sector bajo: protección à la T1-acotado con la forma de Coulomb
   controlada por ω_low ≈ 0.57ν (constante MEDIDA, estable) y el puente
   diagonal barato (ratio 1.000, E-C1).
3. Sector alto: eliminación por complemento de Schur/Feshbach con
   denominador de energía ≥ (√5−1)c² ≈ 1.24c².
4. **E-C3 (lo que falta):** medir/acotar el acople cruzado
   ‖E_low^{-1/2} W_cross E_high^{-1/2}‖ entre sectores — el último
   número del teorema.

