# Population only, short runs: what growth per division does each death rate
# actually give? Extrapolating this from one point is what cost the last probe.
source("load_origin.R")
st <- list(be_mutations_per_target_per_division=0,
  nuc_insertions_per_target_per_division=0, nuc_deletions_per_target_per_division=0,
  mt_mutations_per_genome_per_division=0, num_mt_genomes=0,
  bc_bg_insertion_prob_per_division=0, bc_bg_deletion_prob_per_division=0,
  bc_substitution_model="JC", bc_sub_model_params=list(0),
  mt_substitution_model="JC", mt_sub_model_params=list(0),
  mt_bg_insertion_prob_per_division=0, mt_bg_deletion_prob_per_division=0)
for (d in c(0.45, 0.50, 0.52, 0.55, 0.58)) {
  p <- list(num_init_cells=200, sim_length=list(8), random_seed=1,
    bc_length=40L, mito_genome_length=16600, max_bc_ints_per_cell=list(1),
    bc_nuc_composition=list(frac_a=1,frac_g=0,frac_c=0,frac_t=0),
    be_conversion_pattern="A --> G",
    be_targets=list(num_targets=2, config="S:1:10",
      edit_rate_class_fractions=list(high=1,medium=0,low=0),
      editing_window=list(size=0,decaying=FALSE,close_after_edit=FALSE)),
    editing_induction=list(timepoint=0,num_cells=NULL,frac_cells=1),
    differentiation_induction=list(timepoint=0,num_cells=NULL,frac_cells=1),
    cell_type_dict=list(founder_cell_type="progenitor",
      cell_type_params=list(progenitor=list(cell_cycle_length=1,
        death_per_cell_cycle_prob=d, bc_invariant_sites=0, mt_invariant_sites=0,
        induced_editing_params=st, uninduced_editing_params=st)),
      uninduced_transition_matrix=list(list(1)), induced_transition_matrix=list(list(1))))
  out <- file.path(tempdir(), sprintf("g%g", d)); dir.create(out, showWarnings=FALSE)
  ok <- tryCatch({
    run_gillespie_lineage_pipeline(p, params_path=NULL, output_dir=out,
      overrides=list(progress=FALSE, seed=1, modalities="barcode"))
    TRUE }, error=function(e) FALSE)
  if (!ok) { cat(sprintf("death %.2f -> failed/extinct\n", d)); next }
  n <- data.table::fread(file.path(out,"lineage_nodes.csv.gz"), showProgress=FALSE)
  a <- n[n$alive_at_end %in% c(TRUE,"TRUE"),]
  gen <- median(as.numeric(a$generation))
  cat(sprintf("death %.2f -> %6d alive from 200, generation %.0f, growth/division %.4f\n",
              d, nrow(a), gen, (nrow(a)/200)^(1/gen)))
}
