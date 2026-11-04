#!/bin/bash
# make_transfer_bundle.sh -- run on CLUSTER 2, in the project root (the directory that contains run_1 ... run_10).
# Packs ONLY the files the new runs need into md_transfer_bundle.tgz (a few hundred MB at most) and prints sha256.
# Nothing is moved or deleted.
set -uo pipefail
REQ=(
  common.mod
  potential.mod
  potentials/library.meam
  potentials/CaZnMg.meam
  run_1/quenched_matrix_300K.data      # glass for Composition 1 (xCa 0.05): matched-geometry series, replicates, 0 % case, large cell
  run_5/quenched_matrix_300K.data      # glass for Composition 5 (xCa 0.03): the missing matched-geometry control
  run_6/quenched_matrix_300K.data      # glass for Composition 6 (xCa 0.07)
  run_1/AN_composite_relaxed_300K.data # temperature sweep
  run_3/AN_composite_relaxed_300K.data # strain-rate check, xZn = 0.40
  run_4/AN_composite_relaxed_300K.data # strain-rate check, xZn = 0.25
)
OPT=( pymatgen_offline_bundle.zip PYMATGEN_SETUP.md )     # only if cluster 1 cannot pip-install pymatgen
miss=0; list=()
for f in "${REQ[@]}"; do f=${f%% *}
  if [ -f "$f" ]; then list+=("$f"); printf "ok       %8s  %s\n" "$(du -h "$f" | cut -f1)" "$f"
  else printf "MISSING            %s\n" "$f"; miss=1; fi; done
for f in "${OPT[@]}"; do [ -f "$f" ] && { list+=("$f"); printf "optional %8s  %s\n" "$(du -h "$f" | cut -f1)" "$f"; }; done
[ $miss -eq 0 ] || { echo "Fix the MISSING files first (nothing packed)."; exit 1; }
tar czf md_transfer_bundle.tgz "${list[@]}" && sha256sum md_transfer_bundle.tgz | tee md_transfer_bundle.sha256
echo "On cluster 1:  tar xzf md_transfer_bundle.tgz   (in the new project root), then ./verify_inputs.sh"
