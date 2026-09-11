#!/usr/bin/env bash
set -euo pipefail
task_fixture=$(mktemp -d "${TMPDIR:-/tmp}/lheh5-test.XXXXXXXX")
echo "LHEH5 test artifacts: $task_fixture"
python3 "${SCRAM_TEST_PATH}/makeLHEH5Fixtures.py" "$task_fixture"
testLHEH5Reader "$task_fixture"
for task_case in basic skip runs shower loss; do
  task_args=()
  task_files=events3
  task_count=3
  case "$task_case" in
    skip) task_files=events0,events3,events0,events3; task_args=(skip=4 maxEvents=1); task_count=1 ;;
    runs) task_files=events3,different_run; task_count=6 ;;
    shower) task_args=(shower=True) ;;
    loss) task_files=metadata; task_args=(allowLoss=True); task_count=1 ;;
  esac
  for task_format in xml hdf5; do
    task_suffix=lhe
    if [[ $task_format == hdf5 ]]; then task_suffix=h5; fi
    task_inputs=()
    IFS=',' read -ra task_names <<< "$task_files"
    for task_name in "${task_names[@]}"; do
      task_inputs+=("file:$task_fixture/$task_name.$task_suffix")
    done
    task_input_arg=$(IFS=,; echo "${task_inputs[*]}")
    cmsRun "${SCRAM_TEST_PATH}/lheInputBackend_cfg.py" encoding="$task_format" \
      inputFiles="$task_input_arg" outputFile="$task_fixture/$task_case-$task_format.root" \
      "${task_args[@]}" > "$task_fixture/$task_case-$task_format.log" 2>&1
  done
  task_check=()
  if [[ $task_case == shower ]]; then task_check=(--shower); fi
  if [[ $task_case == loss ]]; then task_check=(--loss); fi
  python3 "${SCRAM_TEST_PATH}/checkLHEInputEDM.py" "$task_fixture/$task_case-xml.root" \
    "$task_fixture/$task_case-hdf5.root" "$task_count" "${task_check[@]}"
done
