#!/bin/bash
# Run the mhdFoam Hartmann sweep over Ha in {1,5,10,20,50}.
# With rho=nu=sigma=mu=a=1 the Hartmann number equals the applied field: Ha = B0.
# Each case: copy template -> set B0=(0,Ha,0) -> blockMesh -> mhdFoam -> sample U(y) at x=10.
# Usage: bash run_openfoam_sweep.sh ["1 5 10 50"]   (default skips 20 if already run)
set -u
cd "$(dirname "$0")/openfoam" || exit 1
source /usr/share/openfoam/etc/bashrc 2>/dev/null

HALIST="${1:-1 5 10 50}"
for HA in $HALIST; do
  CASE="Ha_${HA}"
  echo "=================== Ha = ${HA} (B0=${HA}) ==================="
  rm -rf "$CASE"
  cp -r template "$CASE"
  cd "$CASE" || exit 1
  foamDictionary -entry internalField -set "uniform (0 $HA 0)" 0/B >/dev/null 2>&1
  foamDictionary -entry "boundaryField.lowerWall.value" -set "uniform (0 $HA 0)" 0/B >/dev/null 2>&1
  foamDictionary -entry "boundaryField.upperWall.value" -set "uniform (0 $HA 0)" 0/B >/dev/null 2>&1
  blockMesh > log.blockMesh 2>&1
  echo "  blockMesh exit=$?  ($(grep -m1 nCells log.blockMesh | tr -s ' '))"
  mhdFoam > log.mhdFoam 2>&1
  rc=$?
  echo "  mhdFoam exit=${rc}  last: $(grep '^Time = ' log.mhdFoam | tail -1)"
  # sample U(y) at the latest converged time
  postProcess -func sample -latestTime > log.sample 2>&1
  prof=$(ls postProcessing/*/*/line_centreProfile_U.xy 2>/dev/null | tail -1)
  echo "  profile: ${prof:-NONE}"
  cd ..
done
echo "SWEEP DONE"
