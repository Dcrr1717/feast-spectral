# PLAN PAPER C — El teorema completo: T1-Coulomb por compresiones continuas

**Título tentativo:** *A no-pollution theorem for cusp-adapted
kinetically balanced discretizations of the Coulomb–Dirac operator.*

**Tesis en una frase:** probar la Conjetura 7.2 del Paper T — la
protección de las familias balanceadas en base γ para V = −Z/r — usando
compresiones del operador libre CONTINUO (Lewin–Séré) y la desigualdad
de Hardy relativista en sentido de formas, en el rango de acoplamiento
donde medimos que la cota de forma continua converge y es < 1.

**Naturaleza:** track funcional de César (Claude: verificación numérica,
caza de contraejemplos bajo protocolo math_daniel, y formalización Lean
de los lemas que vayan cayendo). **Revista:** decidir al final según
alcance — Numerische Mathematik / IMA JNA si queda en ν≲0.75; algo
mayor si llega al régimen distinguido. Alternativa: fusionarlo como v2
del Paper T (compuerta G-C3). **Autoría: César solo.**

---

## Estado de partida (todo probado/medido ya)

- Herramientas P1/P2 (FINDINGS §13): ventana por resolvente e
  inclusión condicional por formas — correctas pero con hipótesis
  DIVERGENTES respecto al libre discreto (a_N ~ √N, α_N ~ log N).
- Diagnóstico de imposibilidad: (ΠH₀Π)⁻¹ ≠ ΠH₀⁻¹Π; el libre discreto
  tapa momento en σ_max ≈ 0.64βN. La ruta continua es OBLIGADA.
- E-T2: la cota de forma medida contra |H₀| de referencia grande
  CONVERGE (0.8633→0.8687 con incrementos a la mitad) y queda < 1
  hasta ν ≈ 0.75 a nb práctico.
- E-T3: mapa del impostor s=+1 (nace en el borde inferior, bucea con
  nb) = test de sharpness listo.
- Con Golub–Welsch log-espacio (FINDINGS §16) las referencias pueden
  crecer sin techo → mejores proxies del continuo para C2.

## Fase C0 — Formalización del marco (César, 1-2 semanas de fondo)

- Definir con precisión el objeto "compresión continua" en nuestro
  setting: base γ {f_j} ⊂ H^{1/2} (forma-dominio de |D₀|), Π_N el
  proyector, y las cantidades α(ν) := ‖|D₀|^{−1/2} V |D₀|^{−1/2}‖
  restringida al span creciente. Conectar con las condiciones (i)-(iii)
  de Lewin–Séré Proc. LMS 100:864.
- Recopilar las constantes sharp por canal de la desigualdad de Hardy
  relativista ⟨f,(Z/r)f⟩ ≤ C_κ⟨f,|D₀|f⟩: literatura DES / Tix /
  Burenkov–Evans; el techo de Talman ν≈0.9058 como valor crítico del
  canal |κ|=1. (E-T2 sugiere que la constante efectiva medida está
  entre Talman y Kato — entender POR QUÉ es parte de C0.)

## Fase C1 — El lema puente discreto-continuo (el crux técnico)

El obstáculo identificado: la fibra-diagonalidad de L2 vale para
|H₀ᴺ| = |ΠD₀Π| pero la cota de Hardy vive en Π|D₀|Π, y ambos difieren.
Dos rutas a explorar EN PARALELO (una teórica, una asistida):

- **Ruta (a) — desigualdad de operador:** x↦x² es operator-convex ⇒
  (ΠD₀Π)² ≤ ΠD₀²Π; con la monotonía de √ ⇒ |ΠD₀Π| ≤ (ΠD₀²Π)^{1/2}.
  Preguntar si (ΠD₀²Π)^{1/2} basta como intermediario entre la forma
  de Hardy y la estructura de fibras (D₀² es diagonal por canales y
  su compresión es computable: −d²/dr² + centrífugo + c⁴).
- **Ruta (b) — compresiones LS directas:** verificar las hipótesis del
  Teorema de Lewin–Séré con Λ± continuos y el defecto por modo de L1'
  (las fugas t³/t son exactamente esas compresiones cuantificadas);
  el trabajo es controlar el término ⟨Λ∓ψ, V Λ±ψ⟩ con Hardy + L1'.
- E-C1 (Claude): medir numéricamente la brecha |ΠD₀Π| vs (ΠD₀²Π)^{1/2}
  vs Π|D₀|Π en referencias grandes (ahora sin techo) — decide qué ruta
  tiene constantes viables ANTES de invertir semanas.

**COMPUERTA G-C1:** si ninguna ruta da constantes < 1 en ν ≤ 0.6
siquiera, replantear el enunciado (protección relativa a ventana más
pequeña) antes de seguir.

## Fase C2 — El teorema en el rango medido (ν ≲ 0.75)

- Enunciado objetivo: para la familia balanceada γ-basis sin truncar y
  V = −Z/r con ν < ν₀ (ν₀ el valor donde la cota de forma continua
  llega a 1), el Galerkin no produce autovalores bajo el centro de la
  brecha (versión cuantitativa con margen explícito).
- E-C2 (Claude): batería de verificación del enunciado exacto en
  discretizaciones reales (patrón t1_bounded_window.py) + búsqueda
  adversarial de contraejemplos en el borde ν→ν₀.

## Fase C3 — Sharpness y el impostor (cierre empírico)

- La trayectoria E-T3 del impostor anti-balanceado como test: el
  teorema debe EXCLUIR a las balanceadas exactamente donde el
  anti-balanceado SÍ produce el barrido. Cuantificar la brecha
  t_eff(nb) > t_max(nb) a Coulomb que quedó medida.
- **COMPUERTA G-C3 (decisión editorial):** ¿paper propio o v2 del
  Paper T? Criterio: si C2 cierra con ν₀ ≥ 0.75 y C1 produjo un lema
  puente citable por sí mismo → paper propio; si queda en ν₀ ~ 0.5 →
  sección nueva del Paper T en la revisión.

## Fase C4 — Régimen alto (estiramiento, NO bloqueante)

- Extender hacia el techo de Talman con las formas descompuestas de
  DES, y conjeturar/atacar el régimen distinguido ν→1 con la
  compatibilidad r^γ (la observación "el prefactor selecciona la
  extensión distinguida" de FINDINGS). Este es el premio Q1 gordo;
  se intenta solo si C2 cerró.

## Riesgos

| Riesgo | Mitigación |
|---|---|
| El lema puente no da constantes útiles | E-C1 lo detecta ANTES de invertir teoría; compuerta G-C1 |
| El enunciado útil requiere truncación controlada (nb finito) | formular con defecto por modo de L1' + tasa σ_max de L2 (ya exactos) |
| Tiempo de César (es el cuello real) | fases de 1-2 semanas con entregable mínimo cada una; Claude mantiene el andamiaje vivo |

**Presupuesto:** semanas-meses (ritmo de César). Depende de: nada
externo. Todo el andamiaje numérico está committeado y sin techos.
