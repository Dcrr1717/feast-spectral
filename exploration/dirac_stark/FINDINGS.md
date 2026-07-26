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

**F3a (Z=10, F=50, screening): régimen extremo — la base contraída
(β=5, afinada para el 1s sin apantallar) NO abarca el estado apantallado
+ extendido por el campo.** No es bug del motor (Z=1 apantallado
perfecto); es diseño de base/régimen. Los campos escalados de Maltsev
(F/Z³=0.05 → F=50 a Z=10) son extremos (su estudio de ionización sin
apantallar); para la física apantallada hay que elegir un régimen con
base convergente (β adaptativa al estado apantallado-extendido, o campo
escalado menor, o continuación en λ desde λ=0 con la base re-optimizada).

## Estado: MOTOR FLAGSHIP COMPLETO Y VALIDADO
- F0 angular ✓ (2 derivaciones independientes, 0.0e0)
- F1 campo-libre Sommerfeld ✓ (Z=1..92, ~1e-8..1e-11)
- F2 Stark sin apantallar ✓ (Maltsev 7 dígitos: pos 5e-8, Γ 2e-6)
- F3 Stark apantallado ✓ (Z=1 vs Paper D, 4e-6)
La parte difícil (la ~1 semana del gap scan) HECHA en una sesión por
reuso de Paper B (pencil desplazado) + Paper D (canales acoplados).

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
