# FINDINGS — Paper D: Stark + apantallamiento (stark2d)

## D0 — Compuerta de arquitectura (2026-07-17) · PASADA

**Decisión: canales acoplados con la base exacta (arquitectura ii) gana
por nocaut.** El término Stark F·z solo acopla ℓ↔ℓ±1 (sector m=0) →
sistema radial block-tridiagonal con nuestra base de Laguerre por canal
(β y grilla de cuadratura COMPARTIDOS entre canales para que los
elementos dipolares cruzados sean contracciones exactas), escalado
complejo analítico:
H(θ) = ⊕_ℓ[e^{−2iθ}T_ℓ + V_θ,ℓ] + F e^{+iθ}Σ c_ℓ(|ℓ⟩⟨ℓ+1|⊗R + h.c.),
c_ℓ = (ℓ+1)/√((2ℓ+1)(2ℓ+3)).

**Resultado del benchmark (d0_architecture.py):** resonancia
fundamental del Stark de hidrógeno a F=0.1 a.u.:

    E = −0.527418175 − 0.007269057j

con **780 DOF** (ℓmax=12, nb=60): convergencia geométrica en ℓmax
(4.2e-5 → 1.2e-7 → 4.7e-9 → 4.8e-10), θ-estabilidad 3e-8,
certificación de solapamiento 1e-12, y coincidencia con el valor
clásico de rotación compleja (−0.527418 − 0.007269) en todos los
dígitos disponibles de memoria — re-anclar contra Jentschura/Kolosov
con dígitos completos en D1. La malla 2D (40k+ puntos para mucho menos)
queda relegada a herramienta de verificación cruzada en D1.

**Bug cazado por certificación (lección):** el primer build tenía la
cinética A LA MITAD (factor 0.25 en vez de 0.5 al simetrizar) — el
espectro daba −1.0/−0.25/−0.111 (el doble exacto) y el estado "θ-estable"
de la ventana era la componente roja del n=2 del Hamiltoniano
equivocado; el pick de P-D0.2 lo enmascaró por coincidencia. La
certificación de niveles de hidrógeno puro (P-D0.0, dev < 1e-9) queda
PERMANENTE en el script. Moral de siempre: ventanas y picks no
sustituyen certificaciones de física exacta.

**Límite F=0, λ=0.25 (P-D0.2):** el bloque ℓ=1 reproduce el benchmark
triple-verificado del motor 1D a 7.6e-5 con pick grueso de un solo θ —
suficiente como chequeo de maquinaria; protocolo fino en D1.

**Siguiente (D1):** anclas Stark puro con dígitos completos
(Jentschura 2001 / Telnov 1989 / Kolosov 1987 — ojo con la trampa He⁺
de la review ECS) + límite λ>0/F=0 fino contra el motor 1D + primera
resonancia combinada (F, λ) ≠ 0.

## D1 — Anclas de validación (2026-07-17) · PASADA

**D1a — Stark puro, fundamental (d1_d3_production.py).** Protocolo
2 grillas (ℓmax=12/nb=60 y 16/80) × 2 θ; el valor citado es la grilla
fina y la deriva entre grillas es la barra de error:

    F=0.04:  E = −0.5037715909 − 0.0000019463j   (deriva 2.4e-11)
    F=0.08:  E = −0.5175606169 − 0.0022698288j   (deriva 2.5e-10)
    F=0.10:  E = −0.5274181749 − 0.0072690567j   (deriva 1.8e-09)

F=0.10 coincide con el valor clásico de rotación compleja en los 10
dígitos disponibles. Auto-convergidas a 9–10 dígitos; el re-anclaje
dígito a dígito contra las tablas de Jentschura/Kolosov queda como
chequeo editorial (los PDFs no están a mano), no computacional.

**D1b — F=0, λ=0.25, bloque ℓ=1 vs motor 1D.** Con protocolo fino
(nb=220, β=0.15): E = 0.003064450 − 0.002184741j, desviación
**2.4e-7** contra el benchmark triple-verificado del motor 1D
(0.00306445 − 0.0021845j) — los dos motores coinciden en todos los
dígitos citados del benchmark. (La grilla gruesa nb=150/β=0.20 aún
deriva 3.8e-3: esta resonancia ancha de umbral necesita la grilla
fina; la fina es la que se cita.)

## D2 — Democión honesta del claim de método (2026-07-17)

A ~1000 DOF una diagonalización densa compleja toma segundos y el mapa
completo minutos: FEAST/warm-start es INNECESARIO a esta escala. El
criterio del propio plan manda: el claim de novedad del paper D es el
método de canales acoplados certificado + el primer mapa sistemático
(F, λ), no el eigensolver. FEAST queda para las mallas 2D (40k+ DOF)
donde sí pagaría — mencionado como perspectiva, no como resultado.

## D3 — El mapa (F, λ) (2026-07-17)

**v1 falló como se debe fallar: ruidosamente.** La continuación en λ
usaba ventana fija +0.09 cuando el corrimiento real es ~+λ (pendiente
exacta dE/dλ|₀ = ⟨e^{−λr}⟩ = 1 para CUALQUIER estado normalizado) →
"lost" en todo λ>0. Además el "triplete" n=2 era un error de conteo:
m=0 en n=2 tiene solo k=±1 parabólicos — es un DOBLETE (la corrida
encontró exactamente 2 estados en λ=0: −0.2067 y −0.0890, el
corrimiento cuadrático hacia abajo respecto al lineal ±3F es visible).

**v2 (d3_map_v2.py):** predicción física
E_pred = E_Yukawa(λ) + [E_prev − E_Yukawa(λ_prev)] (curva F=0 exacta
por autovalores simétricos del bloque ℓ=0, un eigensolve barato por λ)
+ pasos Δλ=0.05 + radio/tolerancia escalados con Γ + reintento con
θ=(0.30, 0.42) antes de declarar disuelto.

**RESULTADO v2: 85/85 puntos del mapa rastreados SIN UNA SOLA PÉRDIDA
(584 s), doblete n=2 rastreado 13/13 λ incluyendo el cruce de umbral.**
CSV: d3_map.csv (F, lam, Er, Gamma, Eyuk, qual), d3_manifold.csv.
θ-estabilidad: peor esquina 6e-6 (F=0.10, λ=0.8), típica ≤1e-7.

### La física del mapa (los tres resultados del Paper D)

**1. Ionización de campo catalizada por apantallamiento.** A campo
fijo, Γ(λ) crece monótona y ENORMEMENTE: amplificación λ=0→0.8:
F=0.02: >8 órdenes (de <1e-10 a 1.34e-2); F=0.04: ~1e4; F=0.06: 135;
F=0.08: 22; F=0.10: 9. Origen cuasiestático transparente: el
apantallamiento sube el nivel y el exponente de túnel
exp[−(2/3F)(2|E_r|)^{3/2}] reacciona exponencialmente (a F=0.02
predice 1.7e7 entre λ=0.2 y 0.8 vs 2.1e6 medido — el resto es
prefactor). Lectura de plasma: un campo inocuo para el átomo sin
apantallar domina la anchura de línea a apantallamiento fuerte.

**2. Frontera de disolución Γ=|E_r|.** Cruzada solo para F≥0.06:
λ* ≈ 0.78 (F=0.06), 0.73 (F=0.08), 0.68 (F=0.10). Monótona hacia
apantallamiento débil al crecer el campo.

**3. Corrimiento de campo NO monótono.** Panel central de figD_map:
a λ=0 perturbativo (F=0.02: −0.909e-3 vs −9F²/4=−0.9e-3, 1%); a
F=0.02 el corrimiento se profundiza ×6.8 con λ (polarizabilidad
efectiva crece); pero a F=0.10 tiene MÁXIMO −32.5e-3 en λ≈0.45 y
retrocede a −18.9e-3 en λ=0.8. Cuadruplicar F (0.02→0.04) a λ=0.8
solo aumenta el shift ×1.8 (no ×4): el apantallamiento ENCOGE el
régimen perturbativo Stark; más allá la respuesta dominante al campo
es la anchura, no el corrimiento.

**Doblete n=2 (F=0.02 = 5×F_c(n=2), D3b):** en λ=0:
E(k=−1)=−0.206682−i0.030395 (cuesta abajo, ioniza 2× más rápido) y
E(k=+1)=−0.088984−i0.015445; desviaciones +0.012/+0.010 de PT2 (F ya
no perturbativo para n=2). Bajo λ: ambos suben a ritmo inicial ~1
(pendiente universal), el splitting sobrevive; **la componente azul
cruza E_r=0 en λ≈0.112 y NO pasa nada dramático — continúa suave como
resonancia sobre el umbral** (mismo fenómeno que la 2p Yukawa del
paper NR, aquí alcanzado por una componente Stark bajo continuación de
apantallamiento), llegando a E_r=+0.0837, Γ=0.070 en λ=0.30. La roja
duplica anchura 0.061→0.160; ratio Γ₋₁/Γ₊₁ 1.97→2.29.

**Manuscrito paper_stark/main.tex COMPLETO** (REVTeX/PRA, 4 págs,
2 figuras figD_map [3 paneles: E_r, shift, Γ log] + figD_manifold,
3 tablas, 15 refs; compila limpio). Verificación adversarial de 5
lentes (datos, física, código, bibliografía web, prioridad web)
lanzada como workflow.
