
reinitialize
load /mnt/external/TC/alphafold3_revision/analysis_mature/01_models_pdb/sitophilus_oryzae_2.pdb, peptide
hide everything
show cartoon, peptide
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set ray_opaque_background, off
bg_color white
orient peptide
zoom peptide, 6
set orthoscopic, on
set antialias, 2
set ray_trace_mode, 1
png /mnt/external/TC/alphafold3_revision/analysis_mature/13_structure_gallery/rendered_png/MC25_sitophilus_oryzae_2.png, width=1200, height=1000, dpi=200, ray=1
quit
