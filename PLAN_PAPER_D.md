# PLAN PAPER D — Stark + apantallamiento: el problema 2D donde FEAST es imprescindible

**Título tentativo:** *Contour-integral eigensolvers for non-separable
atomic resonances: the Stark problem in a Debye plasma.*

**Tesis en una frase:** el átomo hidrogenoide apantallado EN CAMPO
ELÉCTRICO no es separable → el problema es genuinamente 2D (r,θ) o
acoplado en ℓ → la diagonalización densa es inviable y FEAST
disperso complejo-escalado pasa de decorativo a habilitante; la física
nueva es el mapa (F, λ) de los manifolds Stark n=2,3 bajo
apantallamiento (competencia ionización-por-campo vs
deslocalización-por-apantallamiento).

**Revista objetivo:** CPC (track software: el diferenciador es el
paquete abierto certificado; la competencia 1D no publica código).
**Motor:** evolución del motor no-relativista del proyecto.
**Cuándo:** después del B (o en paralelo lento); es el más pesado.

## Literatura verificada (scout 2026-07-16) — POSICIONAMIENTO

**Claims MUERTOS (no reclamar):**
- "Primeras resonancias Stark apantalladas": Yu & Ho, Phys. Plasmas
  (2005, estado base, F≤0.12, pocos λ) y sobre todo **Lin & Ho, J.
  Phys. B 44:175001 (2011)** — CSM + malla Lagrange-Laguerre, campo +
  apantallamiento, base Y excitados, con la narrativa "competencia
  campo-vs-apantallamiento" YA en su abstract. Es el near-miss del
  concepto entero: hay que superarlo explícitamente o el referee
  (plausiblemente Ho) lo tumba.
- "Primer solver de contorno para resonancias CS": Mizusaki–Myo–Kato
  PTEP (2014) ya combinó Sakurai–Sugiura + CSM (nuclear). Y FEAST v4
  documenta subspace-restart como feature estándar → el warm-start se
  vende como CONTRIBUCIÓN DE WORKFLOW FÍSICO (tracking de resonancias
  en el plano (F,λ) con economía de factorizaciones, con evidencia de
  tiempos), no como innovación de álgebra lineal.

**Claims VIVOS (defendibles):**
1. **El mapa sistemático (F, λ)** — toda la literatura previa reporta
   puñados de puntos; el mapa denso con curvas críticas λc(F) no existe.
2. **Los manifolds n=2 y n=3 bajo apantallamiento** (cruces evitados,
   topología del manifold, diagrama de fases ionización-vs-
   deslocalización) — intocados.
3. Precisión benchmark moderna (Yu–Ho 2005 usó base Slater L≤11; Lin–Ho
   tiene 15 años) + paquete abierto certificado.

**Advertencia técnica de oro (Fernández & García, arXiv:1505.01327):**
demostraron que el esquema CRLM (malla Laguerre + rotación compleja) de
Fernández-Menchero & Summers dio anchuras MAL POR ÓRDENES DE MAGNITUD
por precisión insuficiente de elementos de matriz. Nuestro esquema
hereda exactamente ese modo de fallo — y nuestra disciplina de
certificación de ensamblaje es literalmente el antídoto: convertirla en
selling point del paper (sección de certificación en CPC).

**Anclas Stark puro (corregidas por scout):** usar Jentschura (2001) /
Telnov (1989) / Kolosov (1987). OJO: en la review de ECS
(arXiv:2103.05058) la tabla F=0.02–0.06 que parece hidrógeno base es
He⁺ (n=2) — trampa de ~5e-3 a.u. si se usa a ciegas.

**Urgencia:** competencia activa 2022–2026 en el espacio físico exacto
(Das–Das–Ghoshal Phys. Plasmas abril 2025; Mukherjee–Roy 2025; el grupo
GPS incumbente arXiv:2504.17317 publica varios papers/año y su salto a
2D les tomaría meses) → si D se hace, hacerlo rápido.

---

## Fase D0 — Decisión de arquitectura (1 sesión) — LA COMPUERTA GRANDE

Dos arquitecturas compiten; decidir con un benchmark de 1 día en vez de
suposiciones:

- **Opción (i) — malla 2D directa** (la de RUTA_Q1: FD/FEM en (ρ,z) o
  (r,cosθ), complejo-escalada, dispersa; verificado N=80k en 7 s con
  LU dispersa). Pros: simple, escala. Contras: la cúspide r=0 en 2D no
  tiene truco de estructura exacta → convergencia algebraica (¡la ley
  L3b nos persigue!); necesita mallas grandes.
- **Opción (ii) — canales acoplados con la base exacta:** V_Stark = F·z
  acopla ℓ ↔ ℓ±1 solamente ⇒ sistema radial acoplado por bloques
  (tridiagonal en ℓ) con NUESTRA base de Laguerre por canal (cúspide
  exacta, cuadratura exacta). Matriz dispersa por bloques
  N_ℓ × nb; FEAST/shift-invert por bloques. Pros: hereda toda la
  precisión de estructura exacta; matrices mucho menores. Contras:
  truncación en ℓ (el Stark fuerte mezcla muchos ℓ — medir cuántos).
- E-D0: mismo caso de prueba en ambas (Stark puro débil, resonancia del
  estado base a F=0.04-0.1 a.u., benchmarks clásicos de dilatación
  compleja) → decide precisión/DOF/tiempo. **Compuerta G-D0: elegir
  arquitectura con datos.**
- **[RESUELTA 2026-07-17, exploration/stark2d/FINDINGS.md]: gana (ii)
  por nocaut — F=0.1: E = −0.527418175 − 0.007269057j con 780 DOF,
  convergencia geométrica en ℓmax, θ-estab 3e-8, coincide con el valor
  clásico en todos los dígitos. La malla queda como verificación
  cruzada de D1. Bug de cinética a la mitad cazado por la certificación
  de niveles de hidrógeno (ahora permanente, P-D0.0).** Mi prior: (ii) gana en precisión por DOF y
  (i) gana en robustez a F grande; el paper puede usar (ii) como motor
  principal y (i) como verificación cruzada — DOS métodos
  independientes = el estándar de validación del proyecto.

## Fase D1 — Validación con Stark puro (1-2 sesiones)

- E-D1a: resonancias Stark del hidrógeno (λ=0): reproducir los
  benchmarks clásicos de alta precisión (dilatación compleja; valores
  del scout) en E_r y Γ para F ∈ {0.04, 0.08, 0.1} a.u., estado base
  y n=2. Objetivo: 8-10 dígitos.
- E-D1b: límite λ>0, F=0: reproducir NUESTRO propio motor 1D validado
  (λc's y resonancias de Yukawa del CompAMa-engine) desde el código 2D.
  Consistencia interna de las dos ramas del proyecto.
- **COMPUERTA G-D1:** ambos límites a ≥8 dígitos. Sin esto no hay mapa.

## Fase D2 — La novedad de método: continuación warm-start (1 sesión)

- El mapa (F, λ) requiere ~centenares de puntos; cada punto = contorno
  FEAST con factorizaciones LU dispersas. Novedad a cuantificar:
  **reutilizar factorizaciones/subespacios a lo largo del camino de
  continuación** (el subespacio de resonancia cambia poco entre puntos
  vecinos) — medir speedup vs cold-start (objetivo honesto: ≥5×).
- E-D2: benchmark de continuación sobre una línea de F con λ fija.
  Si el speedup es <2×, la "novedad de método" se degrada a "detalle
  de implementación" — reposicionar el paper hacia la física (el mapa)
  y el software. Honestidad ante todo.

## Fase D3 — Producción de física: el mapa (F, λ) (2-3 sesiones)

- E-D3a: estado base: E_r(F,λ), Γ(F,λ) — ¿el apantallamiento estabiliza
  o desestabiliza frente al campo? (competencia entre supresión de la
  barrera Stark y deslocalización Debye — la pregunta física central).
- E-D3b: manifolds n=2, n=3: los sub-niveles Stark bajo apantallamiento;
  cruces evitados; dónde muere cada rama (conexión con λc del motor 1D).
- E-D3c: figuras estrella: superficie Γ(F,λ) del estado base +
  abanico del manifold n=2.

## Fase D4 — Paquete CPC + manuscrito (2 sesiones)

- CPC exige el programa como entregable: empaquetar el motor (1D+2D)
  con la disciplina de certificación como feature documentada, tests,
  y el Program Summary. Zenodo DOI al aceptar.
- Manuscrito con el formato CPC (elsarticle — vendorizar la clase como
  hicimos con REVTeX si MiKTeX sigue roto).

## Riesgos

| Riesgo | Señal | Mitigación |
|---|---|---|
| Truncación en ℓ explota a F grande | E_r no converge al añadir canales | limitar el dominio F del mapa (honesto); la física interesante está a F moderado |
| Benchmarks Stark puros no se reproducen | >1e-6 de desviación | es bug nuestro: los clásicos están a 10+ dígitos; depurar antes de seguir |
| Scoop del grupo pseudospectral 1D | scout encuentra 2D suyo 2024-26 | nuestro diferenciador se mantiene: FEAST+código abierto+certificación; reposicionar énfasis |
| El mapa es físicamente aburrido (monótono) | sin estructura en (F,λ) | improbable (cruces evitados n=2 casi garantizados); si pasa, el paper es método+software CPC igual |

**Presupuesto: 8-12 sesiones** (el más caro). Depende de: nada duro;
conviene que B esté cerrado para reutilizar su rama compleja madura.

---

## Orden global del portafolio (para no perdernos)

1. **B** (6-8 sesiones, cero competencia, cierra la trilogía Dirac) ← SIGUIENTE
2. **C** en paralelo lento (ritmo funcional de César; Claude asiste)
3. **D** al final (el más pesado; hereda la rama compleja de B)

Cada plan tiene compuertas con criterios de aborto — la regla de la
casa: ningún experimento sin predicción pre-registrada, ninguna fase
nueva con la compuerta anterior roja.
