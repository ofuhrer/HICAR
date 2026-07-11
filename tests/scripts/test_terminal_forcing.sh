#!/usr/bin/env bash
# Verify the forcing-reader end-time contract with the bundled integration data.
#
# A run ending exactly at a forcing timestamp must consume that final record and
# finish without looking for another one.  In contrast, missing coverage and an
# end time between records must fail before the model starts.
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <hicar-repository> <hicar-executable>" >&2
    exit 2
fi

hicar_repo=$(cd "$1" && pwd)
hicar_exe=$(cd "$(dirname "$2")" && pwd)/$(basename "$2")
test_data="$hicar_repo/tests/Test_Cases"
input_dir="$test_data/input"
set_nml="$hicar_repo/helpers/example_namelists/set_nml_var.py"

[[ -x "$hicar_exe" ]] || { echo "HICAR executable not found: $hicar_exe" >&2; exit 2; }
[[ -d "$input_dir" ]] || { echo "HICAR test data not found: $input_dir" >&2; exit 2; }

mpiexec_path=""
IFS=':' read -r -a path_dirs <<< "$PATH"
for dir in "${path_dirs[@]}"; do
    if [[ -x "$dir/mpiexec" && ! "$dir" =~ python|conda ]]; then
        mpiexec_path="$dir/mpiexec"
        break
    fi
done
[[ -n "$mpiexec_path" ]] || { echo "mpiexec is required for terminal-forcing regression" >&2; exit 2; }

work_dir=$(mktemp -d "${TMPDIR:-/tmp}/hicar_terminal_forcing.XXXXXX")
cleanup() { rm -rf "$work_dir"; }
trap cleanup EXIT
mkdir -p "$work_dir/output" "$work_dir/restart"

write_list() {
    local list_path=$1
    shift
    : > "$list_path"
    local stamp
    for stamp in "$@"; do
        printf '"%s/forcing/laf%s.nc"\n' "$test_data" "$stamp" >> "$list_path"
    done
}

make_namelist() {
    local name=$1
    local end_date=$2
    local list_path=$3
    local nml="$work_dir/$name.nml"

    cp "$hicar_repo/helpers/example_namelists/alpine_realdata.nml" "$nml"
    python3 "$set_nml" "$nml" start_date "'2017-02-14 00:00:00'" --group general --insert
    python3 "$set_nml" "$nml" end_date "'$end_date'" --group general --insert
    python3 "$set_nml" "$nml" forcing_file_list "'$list_path'" --group forcing --insert
    python3 "$set_nml" "$nml" wait_for_ready_file ".False." --group forcing --insert
    python3 "$set_nml" "$nml" output_folder "'$work_dir/output/'" --group output --insert
    python3 "$set_nml" "$nml" restart_folder "'$work_dir/restart/'" --group restart --insert
    # Keep this I/O-contract regression small and self-contained; the forcing
    # reader is exercised while optional physics support tables are not needed.
    python3 "$set_nml" "$nml" pbl "'none'" --group physics --insert
    python3 "$set_nml" "$nml" lsm "'none'" --group physics --insert
    python3 "$set_nml" "$nml" sfc "'none'" --group physics --insert
    python3 "$set_nml" "$nml" water "'none'" --group physics --insert
    python3 "$set_nml" "$nml" mp "'none'" --group physics --insert
    python3 "$set_nml" "$nml" rad "'none'" --group physics --insert
    python3 "$set_nml" "$nml" terrain_shading ".False." --group rad_parameters --insert
    printf '%s\n' "$nml"
}

run_case() {
    local nml=$1
    local log=$2
    (cd "$input_dir" && OMP_NUM_THREADS=1 "$mpiexec_path" -np 2 "$hicar_exe" "$nml") >"$log" 2>&1
}

# The final available record is exactly the requested end time: this must pass.
write_list "$work_dir/exact.list" 2017021400 2017021401
exact_nml=$(make_namelist exact "2017-02-14 01:00:00" "$work_dir/exact.list")
if ! run_case "$exact_nml" "$work_dir/exact.log"; then
    cat "$work_dir/exact.log" >&2
    exit 1
fi
grep -q "Simulation completed successfully" "$work_dir/exact.log"

# A final timestamp before end_time remains a hard coverage failure.
write_list "$work_dir/missing.list" 2017021400
missing_nml=$(make_namelist missing "2017-02-14 01:00:00" "$work_dir/missing.list")
if run_case "$missing_nml" "$work_dir/missing.log"; then
    echo "Missing terminal forcing unexpectedly succeeded" >&2
    exit 1
fi
grep -q "forcing data does not cover the requested simulation end" "$work_dir/missing.log"

# An end between forcing times requires the next real interpolation endpoint.
write_list "$work_dir/non_aligned.list" 2017021400 2017021401
non_aligned_nml=$(make_namelist non_aligned "2017-02-14 01:30:00" "$work_dir/non_aligned.list")
if run_case "$non_aligned_nml" "$work_dir/non_aligned.log"; then
    echo "Non-aligned end without its interpolation endpoint unexpectedly succeeded" >&2
    exit 1
fi
grep -q "forcing data does not cover the requested simulation end" "$work_dir/non_aligned.log"

echo "Terminal forcing regression passed"
