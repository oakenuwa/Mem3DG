#!/bin/sh
# Regenerate the coupled test dataset.
#
# Prerequisites (build these first):
#   libmem3dg      cmake --build <mem3dg-build>  --target mem3dg
#   cytosim libs   cmake --build <cyto-build>    --target cytosim cytomath cytobase
#   the driver     cmake --build <driver-build>
#
# Point DRIVER at the built coupled_driver binary, then run this from inside
# this directory.
set -e

DRIVER=${DRIVER:-coupled_driver}
MESH=${MESH:-../frame9.ply}

# args: <config> <mesh> <nSteps> <cutoff> <stiffness> <strang> <saveEvery>
#
#   cutoff    0.05  physical interaction range; the membrane and filaments
#                   genuinely overlap at this distance (see README)
#   strang    1     Strang splitting, error ~dt^2
#   saveEvery 20    one exported frame per 20 coupled steps -> 11 frames
"$DRIVER" config.cym "$MESH" 200 0.05 20.0 1 20

echo
echo "membrane frames: $(ls frames/membrane_*.obj 2>/dev/null | wc -l | tr -d ' ')"
echo "fiber frames:    $(ls frames/fibers_*.txt   2>/dev/null | wc -l | tr -d ' ')"
echo "(the two counts must match -- they are written in the same loop)"
