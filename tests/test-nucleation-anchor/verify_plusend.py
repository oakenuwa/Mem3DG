# Verify the plus-end-held nucleation geometry and the steric contrast.
#
# Run with Blender's python (it has numpy):
#   blender --background --python-expr "import sys,subprocess; \
#     subprocess.run([sys.executable,'verify_plusend.py'])"
#
# The polyline in fibers_%04d.txt runs MINUS -> PLUS, so the held (plus) end is
# the LAST point.  Body direction is measured from the held end toward the free
# minus end, which is the direction that should be normal to the membrane.
#
# Check the launch angle at frame 0000, NOT frame 0001: saveEvery is 400 steps,
# so by frame 0001 the filaments have had 0.02 s to reorient under load and read
# ~57 deg even when the launch was exact.
import numpy as np
def load_fibers(p):
    fibs=[];cur=None
    for L in open(p):
        L=L.strip()
        if not L or L.startswith("#"): continue
        if L.startswith("fiber"): cur=[]; fibs.append(cur)
        elif cur is not None: cur.append([float(v) for v in L.split()[:3]])
    return [np.array(f) for f in fibs if len(f)>1]
for d in ("frames-steric","frames-no-steric"):
    F=load_fibers(f"{d}/fibers_0000.txt")
    plus=np.array([f[-1] for f in F]); minus=np.array([f[0] for f in F])
    body=minus-plus; body/=np.linalg.norm(body,axis=1,keepdims=True)
    nin=-plus/np.linalg.norm(plus,axis=1,keepdims=True)
    ang=np.degrees(np.arccos(np.clip((body*nin).sum(1),-1,1)))
    L=np.array([np.linalg.norm(np.diff(f,axis=0),axis=1).sum() for f in F])
    print(f"{d}: n={len(F)} at frame 0000 (just nucleated, len {L.mean():.4f})")
    print(f"   |plus end| {np.linalg.norm(plus,axis=1).mean():.4f}   "
          f"body angle to inward normal: {ang.round(2).tolist()} deg")

print()

import numpy as np, sys
def load_fibers(p):
    fibs=[];cur=None
    for L in open(p):
        L=L.strip()
        if not L or L.startswith("#"): continue
        if L.startswith("fiber"): cur=[]; fibs.append(cur)
        elif cur is not None: cur.append([float(v) for v in L.split()[:3]])
    return [np.array(f) for f in fibs if len(f)>1]

for d in ("frames-steric","frames-no-steric"):
    print("="*60); print(d)
    for fr in ("0001","0010","0030"):
        F=load_fibers(f"{d}/fibers_{fr}.txt")
        if not F: continue
        # polyline order is MINUS -> PLUS.  plus end is the LAST point.
        plus  = np.array([f[-1] for f in F])
        minus = np.array([f[0]  for f in F])
        rp = np.linalg.norm(plus,axis=1); rm = np.linalg.norm(minus,axis=1)
        # body direction = from the HELD (plus) end toward the free (minus) end
        body = minus-plus
        body /= np.linalg.norm(body,axis=1,keepdims=True)
        nin = -plus/np.linalg.norm(plus,axis=1,keepdims=True)   # inward normal at the held end
        ang = np.degrees(np.arccos(np.clip((body*nin).sum(1),-1,1)))
        L=np.array([np.linalg.norm(np.diff(f,axis=0),axis=1).sum() for f in F])
        print(f"  frame {fr}: n={len(F):3d}  len mean {L.mean():.4f}")
        print(f"     |plus end|  (held, should stay ~0.80): mean {rp.mean():.4f}  min {rp.min():.4f}  max {rp.max():.4f}")
        print(f"     |minus end| (free, should move inward): mean {rm.mean():.4f}  min {rm.min():.4f}")
        print(f"     body angle to INWARD normal at held end: mean {ang.mean():.1f} deg  max {ang.max():.1f} deg")

print()

import numpy as np, sys
def load_fibers(p):
    fibs=[];cur=None
    for L in open(p):
        L=L.strip()
        if not L or L.startswith("#"): continue
        if L.startswith("fiber"):
            cur=[]; fibs.append(cur)
        elif cur is not None:
            cur.append([float(v) for v in L.split()[:3]])
    return [np.array(f) for f in fibs if len(f)>1]

RAD=0.01; RANGE=2*RAD
for d in ("frames-steric","frames-no-steric"):
    F=load_fibers(f"{d}/fibers_0030.txt")
    n=len(F)
    mins=[]
    for i in range(n):
        for j in range(i+1,n):
            D=np.sqrt(((F[i][:,None,:]-F[j][None,:,:])**2).sum(-1))
            mins.append(D.min())
    mins=np.array(mins)
    npair=len(mins)
    print(f"{d}")
    print(f"   filament pairs                        : {npair}")
    print(f"   pairs closer than {RANGE:.3f} (steric range) : {(mins<RANGE).sum()}  ({100*(mins<RANGE).mean():.1f}%)")
    print(f"   pairs closer than {RAD:.3f}                : {(mins<RAD).sum()}  ({100*(mins<RAD).mean():.1f}%)")
    print(f"   min / 1st pct / 5th pct               : {mins.min():.5f} / {np.percentile(mins,1):.5f} / {np.percentile(mins,5):.5f}")
    print(f"   median pair distance                  : {np.median(mins):.5f}")
