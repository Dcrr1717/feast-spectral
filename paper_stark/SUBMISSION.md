# Paper D — notas de envío

## Revista objetivo

1. **Physics of Plasmas (AIP)** — fit temático PERFECTO: Yu–Ho 2005,
   Bahar–Soylu 2014 y Das–Das–Ghoshal 2025 (el registro previo del
   problema combinado) están todos ahí. Track subscription = $0 APC.
   La audiencia que usa estos números (espectroscopía de plasmas
   densos) lee esta revista.
2. **Physical Review A** — consistencia con Papers A y B (mismo
   formato REVTeX ya usado); el ángulo método-certificado encaja en
   la sección de física atómica computacional. Subscription $0.
3. J. Phys. B — alternativa Q2 (Lin–Ho 2011 está ahí).

Recomendación: PoP primero (audiencia natural del mapa), PRA como
plan B. SIEMPRE track subscription/non-OA.

## Verificación bibliográfica (2026-07-17, web)

- Kolosov 1987: ✓ JPB 20, 2359, DOI 10.1088/0022-3700/20/11/008.
- Telnov 1989: ✓ JPB 22, L399 — DOI corregido a .../22/14/002.
- Jentschura 2001: ✓ PRA 64, 013403 (título exacto confirmado).
- Yu–Ho 2005: CORREGIDO — título real "Stark shifts and widths of a
  hydrogen atom in Debye plasmas", autores A. C. H. Yu y Y. K. Ho.
- Das–Das–Ghoshal 2025: CORREGIDO — título real "Hydrogen atom under
  the joint influence of non-ideal plasma and a static electric
  field", autores B. Das, N. Das, A. Ghoshal; PoP 32, 043508 ✓.
  **PENDIENTE: confirmar DOI 10.1063/5.0250020 en la página del
  artículo antes de enviar.**
- Bahar–Soylu 2014: AÑADIDO (PoP 21, 092703, DOI 10.1063/1.4894684) —
  niveles ligados con campo débil (MGECSC + AIM), citado como parte
  del registro previo.
- Lin–Ho 2011: verificado por el scout del plan (JPB 44, 175001).
- Clásicos (Aguilar–Combes, Balslev–Combes, Herbst, Ho 1983,
  Moiseyev 1998, Rogers–Graboske–Harwood 1970): estándar, sin señales
  de error.

## Pendientes al enviar

- Reemplazar companionNR/companionA con arXiv IDs reales.
- Confirmar DOI de DasGhoshal2025 (arriba).
- Cross-check editorial de las anclas F=0.04/0.08 contra las tablas
  de Kolosov/Telnov dígito a dígito (nuestros valores son
  auto-convergidos con deriva ≤2e-9; F=0.10 ya coincide con el valor
  clásico a 10 dígitos).
- KarHo2007 está en refs.bib sin citar (bibtex la ignora) — borrar o
  citar según revisión final.

## Verificación de números (2026-07-17)

Script inline (scratchpad verify_paper_d.py): 70/73 chequeos PASS
contra los CSV; 2 falsas alarmas de redondeo en la 6ª cifra; 1 real
corregido (peor spread θ es 6.5e-6 → el texto ahora dice "below
7e-6... 73% of the grid"). Física verificada: PT2/PT4 Stark, exponente
de túnel, frontera de disolución interpolada (0.778/0.727/0.682),
cruce azul λ=0.1119, PT2 n=2 con desviaciones +0.0119/+0.0096,
pendientes iniciales ∈(0.9,1.0) los 5 F, F=5.12·F_c(n=2).
