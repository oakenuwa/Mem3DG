#!/bin/sh
# Regenerate the CLOSED-membrane (icosphere) coupled dataset.
# See run.sh for the open-cap variant and for build prerequisites.
set -e

DRIVER=${DRIVER:-coupled_driver}
MESH=${MESH:-../icosphere.obj}

# args: <config> <mesh> <nSteps> <shell> <stiffness> <strang> <saveEvery> <outDir> <mode>
#
#   mesh       icosphere.obj: R=0.8, 642 verts, 1280 faces, CLOSED (Euler 2)
#   cutoff     0.05
#   stiffness  20.0  large enough to matter; see README
#   saveEvery  10    10000/10 -> 1001 frames (~63 MB of OBJ; reduce if too big)
#   mode       1     one-sided wall: confine filaments inside, push membrane out
#                    (0 = symmetric repulsion, which does NOT confine)
"$DRIVER" config-icosphere.cym "$MESH" 10000 0.05 20.0 1 10 frames-icosphere 1

echo
echo "membrane frames: $(ls frames-icosphere/membrane_*.obj 2>/dev/null | wc -l | tr -d ' ')"
echo "fiber frames:    $(ls frames-icosphere/fibers_*.txt   2>/dev/null | wc -l | tr -d ' ')"
