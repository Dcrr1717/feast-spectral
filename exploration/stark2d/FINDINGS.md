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
