# Gap scan Q1 2023–2026 — próximos papers (2026-07-24)

Deep-research de 9 agentes (8 scouts direccionales con ≥6 búsquedas web
c/u + abogado del diablo con spot-checks adversariales propios; 183
llamadas de herramientas). Ningún veredicto refutado. Resumen ejecutivo
priorizado; evidencia completa en el output del workflow (sesión
2026-07-24).

## Secuencia recomendada (el DA la razonó por reloj de scoop)

1. **ECSC + Stark (mapas)** — semanas 1–2.5 → J. Phys. B / PRE
2. **Stark relativista apantallado (Dirac+Debye+F)** — semanas 3–7 → PRA ← *flagship*
3. **Floquet no-hermítico + apantallamiento** — semanas 8–11 → PRA

Lógica: el ECSC es la victoria rápida que reusa el pipeline del Paper D
casi literal Y extiende el paper 2019 del grupo de Guo (sinergia con la
relación Paper B); el relativista lleva el único reloj de scoop real
(extensión obvia de Maltsev 2023 una vez enunciada; arXiv:2603.05865
2026 muestra que el espacio "relativista+screening+Stark" se está
calentando); el Floquet tiene el riesgo de scoop más bajo (cero
movimiento 2015–2026) y va último.

## Los 8 veredictos

### #1 ECSC + campos (PARTIAL — el mejor ratio impacto/esfuerzo)
- Existe: Wang-Chen-Guo JPB 52:025001 (2019) = demo del método en
  puntos dispersos; Jiao field-free (EPJD 2021, CTP 2024); PRE
  107:045201 (2023) delimita explícitamente el régimen perturbativo
  (= agujero reconocido por la literatura).
- Falta: mapa (F,λ) denso, manifolds excitados n=2/3, comparación
  Debye-vs-ECSC de la catálisis a apantallamiento igualado, y la física
  cualitativamente nueva: **la cola oscilatoria del ECSC soporta
  resonancias de forma sin análogo Debye → cruces evitados con las
  resonancias Stark** (figura estrella: trayectorias en el plano
  complejo).
- Bonus técnico: ECSC = Re[Yukawa con λ complejo] → posiblemente CERO
  código nuevo de operadores.
- Objeción de árbitro: "incremental, Paper D con otro potencial" →
  encabezar con la rama de resonancias de forma, no con el método.
- Pre-tarea: leer el full text de Wang-Chen-Guo (fijar frontera de
  novedad exacta) y verificar la reducción ECSC-en-argumento-complejo.

### #2 Dirac + Debye + F (PARTIAL — el flagship)
- Existe: Ivanov-Ho PRA 69:023407 (2004, Dirac+F Z=1,10); **Maltsev-
  Tumakov-Popov-Shabaev PRA 108:042813 (2023): Dirac+F, Z=1–82,
  1s/2s/2p₁/₂/2p₃/₂, supresión relativista de anchuras hasta 1 orden —
  SIN screening, ~0 citas, sin follow-up**; Bielińska-Wąz 2004
  (Dirac+Debye ligados); Shi 2014 (Dirac+Yukawa sin campo). El triple
  Dirac+screening+campo: NADIE.
- Física: ¿sobrevive la supresión relativista al apantallamiento?
  competencia screening-sube vs relatividad-aprieta; splitting de
  anchuras 2p₁/₂ vs 2p₃/₂ de resonancias Stark (¡eco directo de nuestro
  Paper B!); fronteras (D,F) por componente de estructura fina.
- Validación regalada: escalera Ivanov-Ho → tablas Maltsev → límite
  D→∞ contra nuestro Paper D.
- Maquinaria nueva: acoplamiento dipolar relativista entre canales κ
  (z acopla κ→−κ±1; ~1 semana de álgebra angular), disciplina de
  balance cinético bajo escalado (nuestra teoría de polución ES el
  selling point).
- Objeción: relevancia física a Z alto + QED ignorada → framing de
  Maltsev verbatim ("Dirac-level benchmark values") + tabla del tamaño
  QED α(Zα)⁴ vs shifts calculados.

### #3 Floquet + screening (PARTIAL — mayor novedad/maquinaria)
- Existe: Floquet no-hermítico maduro SIN screening (escuela Chu);
  Lumb-Prasad 2015 "Floquet+Debye" pero hermítico/probabilidades, no
  cuasienergías complejas. Nadie publica anchuras Floquet apantalladas.
- Física: cierres de canal N-fotón inducidos por apantallamiento
  (discontinuidades de tasa sin análogo), cruces evitados de
  cuasienergías vs D, y **ω→0 reproduce el Paper D** (cross-check
  interno que nadie puede igualar).
- Estructura: canal = (ℓ, n_fotones) — nuestro block-tridiagonal tal
  cual. 2–4 semanas.
- Objeción SERIA: consistencia ω vs ω_p del modelo Debye estático →
  derivar la ventana de validez (ω<ω_p, plasmas densos/fríos) en la
  Sec. II, no en el rebuttal.

### #4 Muónico Dirac + screening (PARTIAL — follow-on barato tras #2)
Poszwa-Bahar-Soylu PoP 2016 (no-rel, perturbativo) sin sucesor. Full
Dirac + tamaño nuclear finito + resonancias certificadas: ~2 semanas
DESPUÉS del paper 2 (reusa el bloque Stark-κ). El contenido no-rel es
rescaling de hidrógeno — declararlo como lema, nunca paper standalone.

### #5 Screening + campo magnético (PARTIAL)
CERO cálculos de anchuras para átomo apantallado+B (solo ligados:
Ly/Hoang/Le PoP 28:063301 a 30 dígitos = benchmark de certificación
regalado). Δℓ=0,±2 diamagnético encaja en nuestros canales. PERO: sin
cross-check de shooting (no separable) y dilation-analyticity con B
requiere Avron-Herbst-Simon — la historia de certificación más débil
del top. Paper 4-5, no antes. Riesgo: Ly/Le a un paso de complex
rotation.

### #6 WDM / disolución de niveles (PARTIAL — capstone 2027)
El framing "proveedor de datos" está MUERTO (verificado: ningún código
de opacidad/line-shape ingiere tablas E_r/Γ externas — todos generan
interno). El sub-gap real: reemplazar la probabilidad de ocupación
Hummer-Mihalas / criterio Inglis-Teller con w_n(n_e,T) físico =
promedio de microcampo de NUESTROS mapas certificados, comparado contra
HM88/SP/EK en condiciones Hoarty/Ciricosta y fotosferas de enanas
blancas (Montgomery 2022 documenta fallas cuantitativas; Bonitz 2025
pide exactamente estos datos). Se fortalece tras papers 1–3. Matter
Radiat. Extremes (Q1, sin APC).

### #7 Nuclear Dirac Woods-Saxon (PARTIAL — opción con kill criterion)
Nadie certifica: Sun PRC 101 (2020) admite discrepancias inter-método
sin adjudicar en resonancias anchas; el zoo de métodos sigue creciendo
(Jost+CSM 2026) sin adjudicación. PERO la comunidad premia métodos, no
certificación, y si los métodos coinciden dentro de nuestras barras el
paper colapsa. Tratamiento: pre-check de 2-3 días (¿discrepan de verdad
en el régimen ancho?) antes de comprometer 4-6 semanas. Nota: polo del
Fermi en Im(r)=(2k+1)πa capa el ángulo de rotación (pieza analítica
nueva). PRC/EPJA. Estratégico por la relación Guo.

### #8 Dos electrones + campo (OPEN — pero NO para nosotros ahora)
El único OPEN formal es la peor jugada: 8–16 semanas, tira todo el
stack validado de un electrón + el cross-check independiente, y entra
a competir en precisión contra la escuela Hylleraas en su cancha
(Kar/Ghoshal/Jiao tooling activo 2024-2025). La restricción no es
apertura, es economía de capacidades. Proyecto deliberado de 3-4 meses
para más adelante, si acaso.

## Nota de competencia transversal

El grupo de Guo aparece en 3 direcciones (#1 extiende su paper, #2
tiene su maquinaria Dirac-Yukawa, #7 es su territorio). La relación de
dos etapas del Paper B se vuelve más valiosa: el mail etapa-2 puede
mencionar que extendimos su JPB 2019.
