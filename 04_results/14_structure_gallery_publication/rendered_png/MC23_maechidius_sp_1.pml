
reinitialize
load /mnt/external/TC/alphafold3_revision/analysis_mature/01_models_pdb/maechidius_sp_1.pdb, peptide
hide everything
show cartoon, peptide
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 0
set cartoon_loop_radius, 0.25
set cartoon_oval_length, 1.20
set cartoon_oval_width, 0.30
set ray_opaque_background, off
bg_color white
orient peptide
zoom peptide, 3
set orthoscopic, on
set antialias, 2
set ray_trace_mode, 1
set ray_shadows, off
png /mnt/external/TC/alphafold3_revision/analysis_mature/14_structure_gallery_publication/rendered_png/MC23_maechidius_sp_1.png, width=1400, height=1100, dpi=250, ray=1
quit
