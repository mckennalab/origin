#!/bin/bash

set -o pipefail
 
base_path=$(pwd)

# initialize param_dir + scDesign3 options
param_dir=""
reference_path=""
celltype_col="cell_type"
celltype_map=""
use_pseudotime=""
sc_ncores=""
sc_max_cells=""

# Arg parsing:
#   -d <param_dir>           directory of JSON param files (required, prompted if absent)
#   -r <reference>           path to reference SCE/Seurat .rds or .h5ad for scDesign3
#   --celltype_col <name>    reference colData column for cell type (default: cell_type)
#   --celltype_map <json>    optional JSON mapping ct1/ct2/... to reference labels
#   --use_pseudotime         enable lineage-depth pseudotime as a covariate
#   --sc_ncores <n>          cores for scDesign3 fit/sim
#   --sc_max_cells <n>       cap simulated cells per stopping point
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d) param_dir="${2:-}"; shift 2 ;;
    -r) reference_path="${2:-}"; shift 2 ;;
    --celltype_col) celltype_col="${2:-}"; shift 2 ;;
    --celltype_map) celltype_map="${2:-}"; shift 2 ;;
    --use_pseudotime) use_pseudotime="--use_pseudotime"; shift 1 ;;
    --sc_ncores) sc_ncores="${2:-}"; shift 2 ;;
    --sc_max_cells) sc_max_cells="${2:-}"; shift 2 ;;
    *) echo "Unknown flag: $1"; exit 1 ;;
  esac
done

# if param_dir string is still empty (ie user did not provide param dir at execution)
if [[ -z "$param_dir" ]]; then
  read -p 'Enter path to directory with json param files: ' param_dir
fi

if [[ ! -d "$param_dir" ]]; then
  echo "Parameter directory does not exist: $param_dir" >&2
  exit 1
fi

shopt -s nullglob
param_files=("$param_dir"/*.json)
if [[ ${#param_files[@]} -eq 0 ]]; then
  echo "No JSON parameter files found in: $param_dir" >&2
  exit 1
fi

for filename in "${param_files[@]}"; do
  SECONDS=0

#  counter=$((counter + 1))
  echo "Now simulating with parameter file: $filename"
  simulation_log=$(mktemp "${TMPDIR:-/tmp}/remote_mito_sim.XXXXXX")
  if ! Rscript "${base_path}/legacy/sim5_code.R" -P "$filename" 2>&1 | tee "$simulation_log"; then
    echo "Simulation failed for $filename; skipping post-processing." >&2
    rm -f "$simulation_log"
    cd "$base_path" || exit 1
    continue
  fi
  most_recent_runid=$(sed -n 's/^UNIQUE_RUN_ID=//p' "$simulation_log" | tail -n 1)
  gillespie_output_dir=$(sed -n 's/^GILLESPIE_OUTPUT_DIR=//p' "$simulation_log" | tail -n 1)
  rm -f "$simulation_log"
  if [[ ! "$most_recent_runid" =~ ^[0-9]+$ ]]; then
    echo "Could not determine the completed simulation run id; skipping post-processing." >&2
    cd "$base_path" || exit 1
    continue
  fi

  echo "Simulation alone took $SECONDS seconds"

  if [[ -n "$gillespie_output_dir" ]]; then
    echo "Gillespie outputs are complete in: $gillespie_output_dir"
    echo "Skipping legacy FASTA/PHYLIP reconstruction; the Gillespie run already wrote lineage trees and character matrices."
    if [[ -n "$reference_path" ]]; then
      echo "Generating covariate-linked scDesign3 profiles for the Gillespie run..."
      sc_args=(--run-dir "$gillespie_output_dir" --reference "$reference_path" --celltype-col "$celltype_col")
      [[ -n "$celltype_map"   ]] && sc_args+=(--celltype-map "$celltype_map")
      [[ -n "$use_pseudotime" ]] && sc_args+=(--use-pseudotime)
      [[ -n "$sc_ncores"      ]] && sc_args+=(--ncores "$sc_ncores")
      [[ -n "$sc_max_cells"   ]] && sc_args+=(--max-cells "$sc_max_cells")
      Rscript "${base_path}/analysis/cli/generate_physicell_sc_profiles.R" "${sc_args[@]}"
    else
      echo "No -r reference provided; skipping scDesign3 single-cell profile generation."
    fi
    echo "Total sim param file took $SECONDS seconds"
    cd "$base_path" || exit 1
    continue
  fi
  
 #  num_param_files=$(ls | wc -l)
  # sleep 10s
   
  # The simulator emits its run id in a machine-readable line. Using that id
  # avoids races with concurrent runs writing under the same output directory.
  cd "${base_path}/output/processed_fastas/" || exit 1

  echo "completed run id: $most_recent_runid"

  # make a new directory specific to this runid that will house the trees generated
  mkdir -p "${base_path}/output/recon_trees/${most_recent_runid}"
  
  cd "$most_recent_runid" || { echo "Failed to enter $most_recent_runid"; exit 1; }
  
  # check if any fastas are present; if not, check for 10 more seconds
  fastas_present=false
  num_times_searched=0
  max_searches=3
  
  # if there's no fasta present after three searches, check to see if score mats exist
  while [[ $fastas_present = false ]] && [[ $num_times_searched -lt $max_searches ]]; do
    ((num_times_searched++))
    if [[ $(find . -maxdepth 1 -name "*.fasta" | wc -l) -gt 0 ]]; then
      fastas_present=true
    else
      echo "waiting for the first fasta file to appear ... (Search ${num_times_searched}/${max_searches})"
      sleep 3
    fi
  done

  if [[ "$fastas_present" == true ]]; then
    
  
    # get all fasta sequences ending in TERM (including varpos, if applicable)
    # exclude any _msa.fasta files
    terminal_fastas=(*_TERM.fasta)

    for terminal_fasta in "${terminal_fastas[@]}"; do
      echo "Processing $terminal_fasta"
	

      # muscle msa the sequences in the full fasta
      msa_fasta_name="${terminal_fasta%.fasta}_msa.fasta"
      muscle -align "$terminal_fasta" -output "$msa_fasta_name"

      treedir="${base_path}/output/recon_trees/${most_recent_runid}/fasta_${msa_fasta_name%.fasta}/"
      mkdir -p "$treedir"

      treefile_name="${treedir}fasta_${msa_fasta_name%.fasta}"
      echo "fasta treefile_name = ${treefile_name}"
  
      # use this new msa fasta as input to iqtree ...

      #  iqtree -s "$aligned_fasta" -m MFP -bb 1000 -alrt 1000 -bcor 0.9 -nstep 80
      #  iqtree -s "$full_fasta" -m MFP -bb 1000 -alrt 1000 -bcor 0.9 -nstep 80
      #  iqtree -s "$full_fasta" -m HKY -bb 1000 -alrt 1000 -nstep 80 -ntmax 20
      iqtree -s "$msa_fasta_name" -m HKY -bb 1000 -alrt 1000 -nstep 80 -nt 40 -pre "$treefile_name"


  
    done

  fi
  

  # unlike the processed_fasta dir, which will always exist because of reference being written to it,
  # there is no guarantee that the score mats subdirs exist
  # if it does exist, cd into it; if it doesn't, it means phylips will not be present
  
  phylip_dir="${base_path}/output/score_mats/${most_recent_runid}/phylips/"
  if [ -d "$phylip_dir" ]; then
    echo "phylip dir exists"
    cd "$phylip_dir" || exit 1

    # check if any phylips are present; if not, check for 10 more seconds
    phylips_present=false
    num_times_searched=0
    max_searches=3
  
    # if there's no fasta present after three searches, check to see if score mats exist
    while [[ $phylips_present = false ]] && [[ $num_times_searched -lt $max_searches ]]; do
      ((num_times_searched++))
      # if [[ $(find . -maxdepth 1 -name "*.phy" | wc -l) -gt 0 ]]; then
      if [[ $(find . -maxdepth 1 -name "*.fasta" | wc -l) -gt 0 ]]; then # temporary solution of saving fastas to phylip dir
        phylips_present=true
      else
        echo "waiting for the first phylip file to appear ... (Search ${num_times_searched}/${max_searches})"
        sleep 3
      fi
    done
  else # if the phylip dir does not exist, phylips_present can be set to false
    phylips_present=false
    echo "phylip dir does not exist"
  fi
    

  if [[ "$phylips_present" == true ]]; then

    echo "pwd here == $(pwd)"

    # all_phylips=$(ls *.phy)
    # temporary solution of saving fastas to phylip dir
    all_phylips=(*.fasta)

    for phy_path in "${all_phylips[@]}"; do
      echo "phy_path = ${phy_path}"
      # treedir="${base_path}/output/recon_trees/${most_recent_runid}/score_${phy_path%.phy}/"
      treedir="${base_path}/output/recon_trees/${most_recent_runid}/score_${phy_path%.fasta}/"
      mkdir -p "$treedir"

      treefile_name="${treedir}score_${phy_path%.fasta}"
      
      # # don't need this conversion anymore since score mats are saved to fastas and not phylips
      # # testing conversion from .phy to .fasta:
      # fasta_path="${base_path}/output/processed_fastas/${most_recent_runid}/${phy_path%.phy}.fasta"

      # python "${base_path}/phy_to_fasta.py" --phylip_path "$phy_path" --fasta_path "$fasta_path"
      # # don't need this conversion anymore since score mats are saved to fastas and not phylips

      
      # msa_fasta_path="${fasta_path%.fasta}_msa.fasta" # REMOVING THIS FOR NOW AS I AM ONLY WORKING WITH SCORE MATS
      # REINTRODUCE IF USING NUCLEOTIDE SEQUENCES
      # muscle -align "$fasta_path" -output "$msa_fasta_path"
      # iqtree -s ${msa_fasta_path} -bb 1000 -alrt 1000 -nstep 80 -nt 40 -pre $treefile_name
      iqtree -s "$phy_path" -st BIN -m MF -pre "$treefile_name"

    done

  fi
  
  # now look at all those reconstructed trees from both fasta and score methods
  cd "${base_path}/output/recon_trees/${most_recent_runid}/"  
  echo "pwd here below == $(pwd)"

  # now iterate through each directory in the urid subdir of recon_trees, read in the treefile in each subdir, and compare to ground truth
  # skip the
  all_recon_dirs=(*/)

  # if no recon dirs were found, can continue to processing for next param file
  if [[ ${#all_recon_dirs[@]} -eq 0 ]]; then
    echo "No tree reconstruction directories found. Param processing complete."
    echo "Total sim param file took $SECONDS seconds"
  fi
  for recon_tree_dir in "${all_recon_dirs[@]}"; do
    # it's possible that iqtree makes a results dir but not a treefile due to an error (e.g. too few sequences)
    # only try to copmare a treefile if it exists
    treefiles=("$recon_tree_dir"/*.treefile)
    if [[ ${#treefiles[@]} -gt 0 ]]; then
      echo "recon_tree_dir = ${recon_tree_dir}"
      cd "$recon_tree_dir" || exit 1
      treefile_name="$(basename "${treefiles[0]}")"
      echo "treefile_name = ${treefile_name}"
      this_savename="${treefile_name%.treefile}"
      global_path=$(pwd)
      Rscript "${base_path}/legacy/compare_trees_call_from_bash.r" -R "$treefile_name" -I "$most_recent_runid" -S "$this_savename" -P "$filename" -T "$global_path"
      cd .. 
    fi
  done



  cd "$base_path" || exit 1

  if [[ -d "${base_path}/output/rf_dist_files/${most_recent_runid}" ]]; then
    Rscript legacy/process_results_from_bash.r -I "$most_recent_runid"
  else
    echo "No RF-distance files found; skipping result aggregation."
  fi

  if [[ -d "${base_path}/output/score_mats/${most_recent_runid}/matrices/" ]]; then
    Rscript legacy/convert_scoremats_to_csvs.r --urid "$most_recent_runid" --score_mat_path "${base_path}/output/score_mats/${most_recent_runid}/matrices/"
  fi

  # Generate single-cell expression profiles for internal + terminal cells via
  # scDesign3, using the reference scRNA-seq dataset given by -r. Skipped if the
  # caller did not supply -r.
  if [[ -n "$reference_path" ]]; then
    echo "Generating single-cell profiles via scDesign3 (run id $most_recent_runid)..."
    sc_args=(-I "$most_recent_runid" -R "$reference_path" --celltype_col "$celltype_col")
    [[ -n "$celltype_map"   ]] && sc_args+=(--celltype_map "$celltype_map")
    [[ -n "$use_pseudotime" ]] && sc_args+=("$use_pseudotime")
    [[ -n "$sc_ncores"      ]] && sc_args+=(--ncores "$sc_ncores")
    [[ -n "$sc_max_cells"   ]] && sc_args+=(--max_cells_per_timepoint "$sc_max_cells")
    Rscript "${base_path}/legacy/generate_sc_profiles_from_bash.r" "${sc_args[@]}"
  else
    echo "No -r reference provided; skipping scDesign3 single-cell profile generation."
  fi

  echo "Total sim param file took $SECONDS seconds"

done
