# Verify the nucleation-patch setup and contrast the two steric conditions.
import numpy as np, sys, os

def load_obj(p):
    V=[];F=[]
    for L in open(p):
        if L.startswith("v "): V.append([float(t) for t in L.split()[1:4]])
        elif L.startswith("f "): F.append([int(t.split("/")[0])-1 for t in L.split()[1:4]])
    return np.array(V), np.array(F)

def load_fibers(p):
    fib=[];cur=None
    for L in open(p):
        L=L.strip()
        if not L or L.startswith("#"): continue
        if L.startswith("fiber"):
            cur=[]; fib.append(cur)
        elif cur is not None:
            cur.append([float(v) for v in L.split()[:3]])
    return [np.array(f) for f in fib if len(f)>=2]

def load_points(p):
    out={};cur=None
    for L in open(p):
        L=L.strip()
        if not L or L.startswith("#"): continue
        if L.startswith("class"):
            t=L.split(); cur=" ".join(t[1:4]); out[cur]=[]
        elif cur: out[cur].append([float(v) for v in L.split()[:3]])
    return {k:np.array(v) for k,v in out.items() if v}

def pt_tri(P,A,B,C):
    AB,AC,AP=B-A,C-A,P-A
    d1=(AB*AP).sum(-1);d2=(AC*AP).sum(-1)
    BP=P-B;d3=(AB*BP).sum(-1);d4=(AC*BP).sum(-1)
    CP=P-C;d5=(AB*CP).sum(-1);d6=(AC*CP).sum(-1)
    va=d3*d6-d5*d4;vb=d5*d2-d1*d6;vc=d1*d4-d3*d2
    den=va+vb+vc;den=np.where(np.abs(den)<1e-300,1e-300,den)
    v=np.clip(vb/den,0,1);w=np.clip(vc/den,0,1)
    s=v+w;f=s>1;v=np.where(f,v/np.where(s==0,1,s),v);w=np.where(f,w/np.where(s==0,1,s),w)
    Q=A+v[...,None]*AB+w[...,None]*AC
    return np.linalg.norm(P-Q,axis=-1)

def report(d):
    print("="*66); print(d); print("="*66)
    V0,F0=load_obj(f"{d}/membrane_0000.obj")
    A,B,C=V0[F0[:,0]][None],V0[F0[:,1]][None],V0[F0[:,2]][None]

    pts=load_points(f"{d}/points_0000.txt")
    nuc=np.vstack([v for v in pts.values()])
    dsurf=pt_tri(nuc[:,None,:],A,B,C).min(1)
    r=np.linalg.norm(nuc,axis=1)
    print(f"nucleators: n={len(nuc)}  dist-to-membrane max {dsurf.max():.2e}")
    print(f"            |pos| {r.min():.4f}..{r.max():.4f}   z {nuc[:,2].min():.4f}..{nuc[:,2].max():.4f}")
    lat=np.linalg.norm(nuc[:,:2],axis=1)
    print(f"            lateral radius max {lat.max():.4f}  -> patch diameter {2*lat.max():.4f} um"
          f"   (target 1.0, z>0.6245: {'OK' if nuc[:,2].min()>0.6245-1e-6 else 'VIOLATED'})")

    # launch direction vs the outward normal at the base, frame 0
    fib=load_fibers(f"{d}/fibers_0000.txt")
    if fib:
        cosines=[]
        for f in fib:
            base=f[0]; tip=f[-1]
            axis=tip-base
            n=np.linalg.norm(axis)
            if n<1e-12: continue
            axis/=n
            rad=base/np.linalg.norm(base)      # outward radial ~ surface normal
            cosines.append(np.dot(axis,-rad))  # +1 == exactly inward normal
        cosines=np.array(cosines)
        ang=np.degrees(np.arccos(np.clip(cosines,-1,1)))
        print(f"filaments @frame0: n={len(fib)}  angle to INWARD normal: "
              f"mean {ang.mean():.2f} deg  max {ang.max():.2f} deg")

    # final frame: spread and length
    last=sorted(f for f in os.listdir(d) if f.startswith("fibers_"))[-1]
    fibL=load_fibers(f"{d}/{last}")
    L=np.array([np.linalg.norm(np.diff(f,axis=0),axis=1).sum() for f in fibL])
    tips=np.array([f[-1] for f in fibL])
    tipr=np.linalg.norm(tips[:,:2],axis=1)     # lateral spread of the tips
    print(f"filaments @{last}: n={len(fibL)}  length mean {L.mean():.4f} max {L.max():.4f}")
    print(f"            tip lateral spread: mean {tipr.mean():.4f}  p90 {np.percentile(tipr,90):.4f}  max {tipr.max():.4f}")
    print(f"            tip z: mean {tips[:,2].mean():.4f}  min {tips[:,2].min():.4f}")
    # min inter-filament point distance, a direct probe of fiber-fiber steric
    allp=np.vstack(fibL); owner=np.concatenate([[i]*len(f) for i,f in enumerate(fibL)])
    sub=slice(None,None,3)
    P=allp[sub]; O=owner[sub]
    D=np.linalg.norm(P[:,None,:]-P[None,:,:],axis=-1)
    same=O[:,None]==O[None,:]
    np.fill_diagonal(D,np.inf); D[same]=np.inf
    print(f"            min distance between DIFFERENT filaments: {D.min():.5f}")
    return dict(n0=len(fib) if fib else 0, nL=len(fibL), Lmean=L.mean(),
                spread=tipr.mean(), dmin=D.min())

a=report(sys.argv[1]); b=report(sys.argv[2])
print("="*66)
print(f"{'':22}{'steric ON':>14}{'steric OFF':>14}")
for k,lab in [("nL","filaments"),("Lmean","mean length"),
              ("spread","tip lateral spread"),("dmin","min inter-fil dist")]:
    print(f"{lab:>22}{a[k]:>14.5f}{b[k]:>14.5f}")
