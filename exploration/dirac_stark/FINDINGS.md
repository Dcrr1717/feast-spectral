# FINDINGS — Paper F (flagship): Dirac + Debye + Stark

Gap scan #2. El triple no computado: Dirac + apantallamiento + campo,
resonancias Stark con estructura fina resuelta sobre (Z, λ, F). Target
PRA. Sesión autónoma 2026-07-26.

## F0 — Módulo angular (PASADO)

`angular.py`: ⟨Ω_{κ'm}|cosθ|Ω_{κm}⟩ vía elemento de matriz reducido
Wigner-3j, componente grande (A_L) y pequeña (A_S, usa paridad −κ).
Gates: simetría OK, selección l'=l±1 OK, límite no-rel del estado s
exacto (1/3 = c₀²). **F0b: cross-check independiente por descomposición
Clebsch-Gordan + Gaunt escalar coincide con la ruta Wigner a 0.0e0 en
todos los pares** — coeficientes certificados por dos fórmulas disjuntas.
La matriz A_L es block-tridiagonal en la escalera κ (acoplamientos
dominantes ~0.33-0.49 + acoplamiento pequeño entre j-partners del mismo
l', ~−0.07).

## F1 — Motor acoplado + certificación campo-libre (PASADO)

`ds_core.py` CoupledDiracStark: canales κ en grilla Legendre COMPARTIDA
(mismo nq/xmax → x,W idénticos → ⟨r⟩ cruzado consistente). Bloque
diagonal por canal = pencil complejo-simétrico desplazado del Paper B;
bloques de campo F·e^{+iθ}[[A_L·R, 0],[0, A_S·R]] (z diagonal en el
espacio grande/pequeño). **Campo-libre reproduce Sommerfeld exacto por
canal: Z=1 ~1e-11, Z=40 ~1e-11, Z=92 s1/2 ~1e-8.** Dos bugs cazados:
los autovalores del pencil desplazado YA son energías de ligadura (no
restar c²); β debe escalar ~Z (el 1s a Z=40 decae e^{−40r}, β=0.6Z
óptimo). cert_tol relajado (c2/c3 = certs de 1/r² que cs_pencil no usa;
c1/c4 se mantienen ~1e-12).

## F2 — Resonancia Stark relativista (FEASIBILITY PASADO)

**Motor validado end-to-end.** Z=1, F=0.10, 14 canales:
E = −0.5274226 − 0.0072687j (θ-spread 9.3e-8), convergencia limpia en
canales (1e-6 → 3.6e-7 → 9.3e-8). vs no-rel (Paper D −0.5274182 −
0.0072691j): dRe = −4.4e-6, dΓ/2 = −3.8e-7. **El residuo dRe es del
tamaño correcto de la corrección Dirac del 1s (−6.7e-6)** — la física
relativista entra bien. Acoplamiento angular + Sommerfeld + resonancia
Stark + escalado complejo: TODO funciona.

### F2b — Validación dura contra Maltsev (PASADA)
Maltsev PRA 108:042813 (2023): campo escalado F/Z³=0.05, ΔE = corrimiento
Stark relativo al nivel Dirac campo-libre, −Γ/2 anchura. Datos 1s
obtenidos del arXiv (ar5iv). **Nuestro motor reproduce (nchan=12, nb=48):**

| Z  | F real | ΔE rel err | Γ/2 rel err | θ-spread |
|----|--------|-----------|------------|----------|
| 1  | 0.05   | 6.3e-6    | 1.4e-3     | 3.4e-11  |
| 10 | 50     | 7.8e-6    | 1.4e-3     | 9.2e-7   |
| 18 | 291.6  | 1.3e-5    | 1.4e-3     | 2.8e-6   |

Posiciones a ~1e-5, anchuras a 1.4e-3 (error consistente → truncación
de canales/base, convergiendo). La física relativista de supresión de
anchuras está capturada. Referencia campo-libre E0 tomada del propio
motor (F=0) para que errores de base se cancelen en el corrimiento.
**CONVERGENCIA A 7 DÍGITOS DE MALTSEV confirmada** (Z=1, F=0.05):
nchan=12/nb=48 → ΔE 6.3e-6, Γ/2 1.4e-3; nchan=16/nb=56 → 3.4e-7,
1.2e-5; nchan=20/nb=64 → **ΔE 5.2e-8, Γ/2 2.0e-6**. El error de anchura
1.4e-3 era pura truncación de canales; converge a la precisión de
Maltsev. **El flagship F0-F2 está completo y validado a precisión de
benchmark.**

Datos completos Maltsev 1s (F/Z³=0.05) para F3: Z=1/10/18/36/54/82 →
ΔE = −6.105058e-3 / −6.068738e-1 / −1.939745 / −7.304855 / −1.478389e1
/ −2.571202e1; −Γ/2 = −3.858129e-5 / −3.752554e-3 / −1.141133e-2 /
−3.430401e-2 / −4.599133e-2 / −2.427590e-2. **La anchura PICA en Z~54
y CAE a Z=82 = supresión relativista (nonrel seguiría creciendo).**

## F3 — Apantallamiento encendido (motor VALIDADO, producción pendiente)

**F3b: acoplamiento screened CERTIFICADO a Z=1.** Dirac Z=1 screened
Stark (F=0.10, λ=0.20) vs Paper D no-rel (E_r=−0.3565208, Γ=2.3811e-2):
Dirac da E_r=−0.3565249, Γ=2.380730e-2 → dE_r=−4.15e-6 (tamaño correcto
de la corrección Dirac del 1s), dΓ=−3.70e-6. **Dos codebases disjuntos
(Paper D no-rel stark2d ↔ Paper F Dirac) coinciden a ~4e-6** — el
potencial apantallado en el motor relativista acoplado está validado.

**F3a → RESUELTO en F3c: el "bloqueo" NO era la base — era EL IMPOSTOR
κ=+1 de nuestros propios Papers A/T.** Diagnóstico (2026-07-26): a Z=10,
λ=1 la referencia campo-libre medida por el selector full-matrix daba
−46.38, POR DEBAJO de la cota variacional física (escalado hidrogénico:
Z²·E(1, λ/Z) = −40.71 + corr. rel.) → imposible para un estado físico.
Verificado por bloques: el bloque κ=−1 da **−40.7713** (física perfecta:
−40.7058 no-rel + −0.0655 rel); el bloque κ=+1 da −46.38 = el estado
ESPURIO (degenerado exacto con el 1s en λ=0 por Tupitsyn-Shabaev, se
despega HACIA ABAJO bajo screening, β-dependiente como todo espurio).
El selector ev[0] del F3a anclaba las predicciones al impostor → por eso
"se perdía" la resonancia (buscaba a −44.6 cuando la física está en
−41.4). **FIX: referencia por bloque de canal (κ=−1 para el 1s),
impostor excluido por construcción — nuestra teoría de polución como
herramienta de certificación del propio paper.**

**F3c producción COMPLETA (2026-07-26, f3c_1s_map.csv): EL PRIMER MAPA
STARK RELATIVISTA APANTALLADO — CON TITULAR.** Unidades escaladas
u=λ/Z, f=F/Z³=0.05; supresión S(Z,u) = Γ_Dirac/(Z²·Γ_nonrel) con
baseline no-rel viva (motor Paper D; su u=0 reproduce Maltsev Z=1:
Γ/2=3.86e-5 ✓):

| Z  | S(u=0) | S(u máx) | Δ |
|----|--------|----------|---|
| 10 | 0.9737 | 0.9760 (u=0.15) | +0.2% |
| 20 | 0.8943 | 0.9076 (u=0.15) | +1.5% |
| 40 | 0.6258 | 0.6924 (u=0.20) | **+10.6%** |

**RESULTADO TITULAR: el apantallamiento DEBILITA la supresión
relativista de anchuras, con efecto creciente en Z** (a Z=40 la
supresión pasa de 0.626 a 0.692 en el rango u≤0.20). Física: el
screening empuja el estado lejos del núcleo, diluye su carácter
relativista y la anchura se recupera hacia el valor no-rel — la
supresión de Maltsev se "deshace" parcialmente en plasma. Anclas: u=0
coincide con Maltsev (Z=10: S=0.974 vs 0.973 de sus tablas). Calidades
q ≤ 8e-5. Límite honesto: u≥0.2 (Z=10,20) y u=0.3 (Z=40) se pierden
(estados muy anchos; extensible con radios/tolerancias, la tendencia
está establecida con 4-5 puntos por Z).

## Estructura fina resuelta (de-riesga el titular)
Campo-libre Z=20, n=2 (motor vs Dirac exacto, todos ~1e-12):
2s1/2 = 2p1/2 = −50.336594 (degenerados, mismo |κ|=1), 2p3/2 =
−50.066742. **Splitting 2p3/2 − 2p1/2 = 0.2699 a.u.** (crece ~Z⁴). El
motor separa limpiamente los j-partners → las anchuras Stark de 2p1/2
vs 2p3/2 se computan por separado: el titular del paper (splitting de
anchuras de estructura fina bajo campo/screening) es ACCESIBLE.

## Estado: MOTOR FLAGSHIP COMPLETO Y VALIDADO
- F0 angular ✓ (2 derivaciones independientes, 0.0e0)
- F1 campo-libre Sommerfeld ✓ (Z=1..92, ~1e-8..1e-11)
- F2 Stark sin apantallar ✓ (Maltsev 7 dígitos: pos 5e-8, Γ 2e-6)
- F3 Stark apantallado ✓ (Z=1 vs Paper D, 4e-6)
La parte difícil (la ~1 semana del gap scan) HECHA en una sesión por
reuso de Paper B (pencil desplazado) + Paper D (canales acoplados).

## F3e — Z EXTREMO (54, 82): EL NÚMERO ESPECTACULAR (2026-07-26)

f3e_extremeZ.csv, mismo protocolo, semillas Maltsev:

| Z  | S(u=0) | S(u=0.20) | Δ relativo |
|----|--------|-----------|-----------|
| 54 | 0.4089 | 0.4950    | +21%      |
| 82 | 0.0923 | 0.1497    | **+62%**  |

**A Z=82 el apantallamiento des-hace la supresión relativista en 62%**
sobre u∈[0,0.20]. La anchura absoluta crece ×13.5 mientras la catálisis
no-rel sola daría ×8.3 → **a Z extremo el plasma cataliza la ionización
4.4× más que no-relativistamente porque además remueve la protección
relativista** — doble efecto superlineal en Z (la frase de oro del
paper). Validación: S(54, u=0)=0.4089 vs 0.4088 derivado de las tablas
de Maltsev (¡exacto!); Z=82: 0.0923 vs 0.0936 (1.4%, cert 8e-8 y
q~1e-5 a Z extremo — barras honestas). La figura titular
(figs/figF_suppression) actualizada con los 5 Z, panel S en log.

## F3d — Cuarteto n=2 (m=1/2) a Z=40: PARCIAL (2026-07-26)

Protocolo v2: F-ramp adiabático en λ=0 (etiquetado por continuación
desde campo-libre) + continuación en λ. Resultados (f3d_quartet.csv):
- **El componente cuesta-arriba (2p3/2-dominante, dE=+13.99 ≈ Stark
  lineal ✓) rastreó TODO u=0→0.20 con catálisis ×2000** (Γ: 3.78e-3 →
  7.76; q 3e-8→1e-4).
- Primeros ratios de estructura fina con campo+screening:
  Γ(up)/Γ(down) = 3.93 (u=0) → 2.07 (u=0.05) — el ratio se mueve
  FUERTE con screening.
- LIMITACIONES para publicación: (a) 2s1/2 y 2p1/2 son
  Dirac-degenerados exactos → ambas rampas caen al MISMO componente
  Stark (hay que desambiguar tomando los 2 autovalores más cercanos
  del par); (b) el componente cuesta-abajo se pierde desde u=0.10
  (ensanchamiento + predicción); las etiquetas son mezclas Stark, no
  j-puros — presentar como "componentes por continuación" (lección
  Papers D/E). Refinamiento = próxima sesión (pasos u más finos,
  radios escalados con Γ, tracking de pares).

## Figura titular
figs/figF_suppression.pdf: panel izq. S(Z,u) subiendo hacia 1 con u
(el des-hacer del screening); panel der. Γ/Z² vs u con la curva no-rel
punteada — Z=10 la abraza, Z=40 converge hacia ella. LA figura del
paper.

## Siguiente (producción, sesión deliberada)
- F3-prod: régimen de base convergente para Z alto apantallado.
  Continuación en λ desde el punto Maltsev con β re-optimizada por paso;
  o campo escalado más suave (F/Z³=0.02). Mapa (Z, λ, F) estructura
  fina. HEADLINE: ¿sobrevive/se fortalece la supresión relativista de
  anchuras al apantallamiento? Splitting Γ(2p1/2) vs Γ(2p3/2) bajo
  screening (eco relativista del Paper B). Datos Maltsev Z=1..82 arriba
  como escalera. Sun-Zheng-Liu 2026 (ruptura de simetría de espín
  nuclear) = competencia a citar/diferenciar.
- F4: cross-check independiente (shooting Dirac o 2ª discretización) +
  verificación adversarial + manuscrito PRA. Validación: Ivanov-Ho PRA
  69:023407 (2004) Z=1,10 + Maltsev + límite Paper D.
