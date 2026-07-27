# Plan Nivel 2 (2026-07-27): cola de vacíos reordenada + problemas abiertos

Reordenado por decisión de César (27-jul): G1 = B primero (rápido), WDM al
final (planificar, posible nube), dos-electrones antes del capstone.
Regla transversal (comentario de César sobre el Floquet): en plasma el
Hamiltoniano escalado/no-hermítico pierde la estructura autoadjunta usual
— cada paper lleva una sección explícita de CONDICIONES DE VALIDEZ
(ventanas de frecuencia/densidad, dominio del operador rotado, criterios
θ-estabilidad), como ya hicimos en D/E/F.

## G1 — Apantallado + campo magnético (PRIMERO, ~1.5-2 sem)
- Física: átomo H (luego ión hidrogenoide) en plasma Debye + B uniforme.
  B∥z conserva m (igual que Stark conserva m); término paramagnético es
  diagonal (m·B/2), el diamagnético B²ρ²/8 acopla ℓ con Δℓ=0,±2 → el
  motor de canales del Paper D se reusa casi entero (solo cambia el
  operador de acoplamiento: ρ²=r²sin²θ en vez de z=r·cosθ).
- Nuevo real: elementos ⟨ℓ'|sin²θ|ℓ⟩ (Gaunt, trivial) + r² radial (ya
  certificado). Resonancias por CS del término diamagnético (¡r²e^{2iθ}
  rota al revés que el Coulomb! — condición de sector específica, la
  "condición específica" de César aplica aquí también).
- Anclas: hidrógeno-en-B de alta precisión (Kravchenko exacto; tablas
  clásicas Rösner-Wunner-Ruder-Herold), límite Landau, λ=0.
- Gate previo (2-3 días): scout de literatura B+Debye (hay AIM/variacional
  para ligados; el hueco esperado = RESONANCIAS y B fuerte certificado).
  Kill-criterion: si alguien ya tiene el mapa (B,λ) con anchuras, pivotar.
- Riesgo conocido: sin cross-check de shooting (no separable); mitigación
  = dos bases/dos θ como en D + límite Landau analítico.

## G2 — Floquet no-hermítico apantallado (~2-4 sem)
- Canal = (ℓ, n_fotón); ω→0 debe reproducir Paper D (validación gratis).
- CONDICIONES ESPECÍFICAS (César): el operador de cuasi-energía
  K = H−i∂_t escalado NO es autoadjunto; teoría rigurosa disponible solo
  no-relativista (Howland-Yajima); ventana física ω vs ω_plasma explícita
  en Sec. II; θ-estabilidad + dos-tamaños-de-bloque-Floquet como criterio
  de aceptación.
- Enganche a problema abierto: conjetura de ionización de
  Costin-Lebowitz (¿forzamiento periódico genérico ioniza completamente?)
  — probada solo en modelos (interacciones puntuales); para
  potenciales 3D de corto alcance (¡Yukawa!) ABIERTA. Nuestro paper puede
  aportar la evidencia numérica controlada: ¿existen estados ligados de
  Floquet (Γ=0) en Yukawa+campo periódico o todo se ioniza? Sección
  dedicada + posible companion matemático.

## G3 — Muónico Dirac apantallado (~2 sem)
- Reuso casi total del motor F (masa reducida m_μ ≈ 206.77 m_e; escalas
  λ y F re-escaladas). Nivel Dirac como benchmark (VP de Uehling domina
  las correcciones reales en muónicos — declararlo como límite honesto,
  igual que Maltsev declara lo suyo).
- Física: los muónicos viven ~200× más cerca del núcleo → sensibilidad
  distinta al apantallamiento; mapa S_μ(Z,u) vs el S electrónico del F.

## G4 — Dos electrones apantallados (~4-8 sem, honesto)
- César: "solo experimentación y validaciones" — la validación sí es
  directa (tablas He/H⁻ en Debye de Ho y grupo Jiao), pero el motor NO es
  trivial: CI complejo-escalado con productos de nuestra base radial +
  integrales e-e (Slater/Condon en Laguerre; fórmula cerrada conocida).
  Scope quirúrgico: estados S (L=0) primero, H⁻ y He.
- Física objetivo: resonancias autoionizantes bajo apantallamiento
  (el clásico 2s2p ¹P de He en Debye) + estabilidad de H⁻ (λc de
  desaparición del único ligado) — único hueco OPEN del gap scan.

## G5 — WDM / Hummer-Mihalas (capstone, planificar aparte)
- Requiere: barridos masivos (Z,λ,F,T) → presupuesto de cómputo real.
  Opciones: servidor en la nube (spot GPU/CPU-hours), o el cluster de
  ESPOCH/Yachay si Manzanilla tiene acceso (preguntarle en el Zoom).
- Entregable: tablas de disolución/occupation probabilities con
  resonancias certificadas = insumo para modelos de opacidad.
- No empezar hasta que 2-3 de G1-G4 estén enviados.

## SUBIR DE NIVEL — problemas abiertos verificados (búsqueda 27-jul)

### N1 (RECOMENDADO, físico-matemático): carga crítica supercrítica EN PLASMA
- Problema abierto famoso: Z_cr≈173 (núcleo extendido) donde el 1s bucea
  al continuo inferior → emisión espontánea de positrones / vacío
  cargado. NUNCA observado; parámetros de resonancia supercrítica se
  calculan por extrapolación/continuación analítica (no directo).
- Nuestro ángulo ÚNICO (nadie lo tiene): (a) base exact-structure ya
  operativa a ν=0.9924; (b) con núcleo extendido cruzamos Z>137 sin
  ambigüedad y computamos la resonancia buceada DIRECTAMENTE con CS
  certificado (no extrapolación); (c) el knob nuevo = apantallamiento:
  Z_cr(λ) — ¿el plasma retrasa o adelanta el buceo? (hay precedente solo
  para B superfuerte, no para Debye); (d) conexión con Papers A/T: el
  hallazgo "el prefactor r^γ selecciona la extensión distinguida" toca la
  ambigüedad matemática ABIERTA de la extensión autoadjunta para ν>1.
- Paper físico (PRA/PRL-able si Z_cr(λ) da sorpresa) + posible teorema
  (track funcional César: ¿qué extensión selecciona la base cuando ν>1?).

### N2 (matemático puro, track César): teorema de rotación compleja Dirac–Stark
- Verificado: Šeba 1988 cubre Dirac + potenciales dilatación-analíticos
  (SIN campo); Herbst-Simon cubre Stark NO-relativista; el caso
  Dirac+lineal NO tiene teorema — nuestro Paper F lo declara
  explícitamente como fundamento pragmático. Probar (aunque sea parcial:
  espectro esencial del rotado, θ-independencia de autovalores discretos,
  límite F→0 → Šeba) = paper en J. Spectral Theory / LMP, patrón Paper T
  (la limitación honesta de un paper físico → teorema propio).
- Toda la evidencia numérica ya existe (mesetas θ del F).

### N3 (evidencia + companion): conjetura de ionización Costin-Lebowitz
- Va DENTRO de G2 (ver arriba). Entrada matemática realista porque
  Yukawa es de corto alcance (más cerca de lo que sus métodos manejan
  que Coulomb).

### N4 (contexto, no comprometerse): zero modes de Loss-Yau
- "¿Qué campos magnéticos admiten zero modes?" sigue abierto (criterio
  sharp parcial JEMS 2024). Relevante como CONTEXTO del G1 (estabilidad
  de la materia en B), pero B uniforme no toca zero modes L² — citar,
  no atacar.

## Secuencia propuesta
1. AHORA: consolidación (Zoom Manzanilla, feedback colaborador, César lee
   E/F, arXiv del A) — nada nuevo se abre hasta eso.
2. G1 (B) con su gate scout → paper 7.
3. N1 supercrítico en paralelo de fondo (gate: núcleo extendido en el
   motor A/B = ~3-4 días; luego Z_cr(λ) es un barrido corto). Si el gate
   pasa, N1 SUBE a flagship y G2-G4 esperan.
4. G2 (Floquet+N3), G3 (muónico), G4 (dos electrones), G5 (WDM, nube).
- César en paralelo permanente: T1-Coulomb (ensamblaje multi-escala) y
  N2 (Dirac-Stark riguroso) — sus dos teoremas.
