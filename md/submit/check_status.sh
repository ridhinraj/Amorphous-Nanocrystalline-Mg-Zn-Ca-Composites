#!/bin/bash
printf "%-14s %-6s %-6s %-10s %-9s %s\n" case build relax tension last_strain tension_wall_time
for d in $(cat cases/*.txt); do
  [ -d "$d" ] || continue
  b="-"; [ -f "$d/log.2_build" ] && b="yes"
  r="-"; ls $d/AN_composite_relaxed_300K.data $d/hcp_relaxed_300K.data >/dev/null 2>&1 && r="yes"
  t="-"; w="-"; ls=""
  if [ -f "$d/log.4_tension" ]; then
    w=$(grep -h "Total wall time" "$d/log.4_tension" | tail -1 | awk '{print $NF}')
    [ -n "$w" ] && t="DONE" || t="partial"
  fi
  [ -f "$d/tension_stress_strain.txt" ] && ls=$(tail -n 1 "$d/tension_stress_strain.txt" | awk '{printf "%.4f",$1}')
  printf "%-14s %-6s %-6s %-10s %-9s %s\n" "$d" "$b" "$r" "$t" "${ls:--}" "${w:--}"
done
