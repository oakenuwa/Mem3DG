"""Is the deformation SHARED with the neighbourhood, or concentrated?"""
import glob,math,sys
def rd(p):
    V=[]
    for l in open(p):
        if l.startswith('v '): V.append(tuple(float(x) for x in l.split()[1:4]))
    return V
def rdf(p):
    out=[]
    for l in open(p):
        if l.startswith('#') or l.startswith('fiber'): continue
        t=l.split()
        if len(t)==3: out.append(tuple(float(x) for x in t))
    return out
d=sys.argv[1]
mems=sorted(glob.glob(d+'/membrane_*.obj')); fibs=sorted(glob.glob(d+'/fibers_*.txt'))
V0=rd(mems[0]); V=rd(mems[-1]); Fp=rdf(fibs[-1])
rows=[(min(math.dist(v,q) for q in Fp), math.dist(v,v0)) for v,v0 in zip(V,V0)]
disp=[x for _,x in rows]
nz=sum(1 for x in disp if x>0.01)
print('  verts moved >0.01: %3d/%d (%4.1f%%)   max %.4f   mean %.4f'%(
    nz,len(disp),100*nz/len(disp),max(disp),sum(disp)/len(disp)))
for lo,hi in [(0,0.15),(0.15,0.3),(0.3,0.6),(0.6,9)]:
    sel=[x for dd,x in rows if lo<=dd<hi]
    if sel: print('    dist[%.2f,%.2f)  %3d verts  mean|disp| %.4f  max %.4f'%(lo,hi,len(sel),sum(sel)/len(sel),max(sel)))
