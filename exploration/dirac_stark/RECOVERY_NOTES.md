# Recuperación del wipe de OneDrive (2026-07-26 ~19:00)

## Qué pasó
Entre ~18:50 y ~19:00, OneDrive vació el árbol local del proyecto
(`OneDrive\Desktop\metodo feast espectral\feast_spectral\` quedó solo
con lo escrito en los últimos minutos) Y borró los paquetes del
Desktop (`papers_overleaf\`, `papers_overleaf_profe.zip`,
`Respuesta_Manzanilla_26-7-2026\`). Los folders antiguos del Desktop
sobreviven — el patrón afecta lo creado/modificado HOY. `.git`
desapareció también → el working tree local está muerto.

## Qué está a salvo
- **GitHub tiene TODO hasta el commit F3e inclusive** (último push
  verificado: el mapa Z=54/82 + figura titular + FINDINGS F3e).
- Sobrevivientes locales asegurados en `C:\Users\dcrr1\feast_safe\survivors\`:
  f3d_triplet.csv (el triplete completo) + figs/figF_triplet.{pdf,png}.
- Los scripts post-F3e no comiteados (f3d_v3 con anclaje, f3f) se
  reconstruyen desde el contexto de la sesión (este directorio).
- La memoria de Claude (~/.claude) no estaba en OneDrive: intacta.

## Qué se perdió temporalmente (recuperable)
1. El working tree → re-clonado desde GitHub a
   `C:\Users\dcrr1\repos\feast-spectral` (FUERA de OneDrive).
2. Los paquetes del Desktop → 100% regenerables desde el repo
   (paper/ compila; los zips se rehacen con los mismos scripts).
3. Restauración TOTAL alternativa: onedrive.com → Papelera de
   reciclaje → Restaurar (recupera también los no-repo).

## Acciones para César
1. Entrar a onedrive.com → Papelera → restaurar lo de hoy (o
   verificar en la papelera de Windows).
2. Revisar por qué OneDrive hizo rollback (¿otro dispositivo?,
   ¿almacenamiento lleno?, ¿pausa/reset de sync?).
3. DECISIÓN RECOMENDADA: mover el repo de trabajo FUERA de OneDrive
   (C:\Users\dcrr1\repos\feast-spectral) y dejar OneDrive solo para
   documentos — git+GitHub ya es el respaldo del código.

## Estado de la física de hoy (todo seguro)
- F3c: primer mapa Stark relativista apantallado Z={10,20,40} ✓ en GitHub
- F3e: Z={54,82} — titular: screening des-hace supresión +62% a Z=82 ✓ en GitHub
- F3d v3: triplete n=2 completo (catálisis ×6000, cruce E_r=0 en u≈0.19,
  ratio FS no-monótono 1.02→0.22→0.33) — CSV+figs en feast_safe,
  scripts reconstruidos aquí, pendiente commit al clon.
