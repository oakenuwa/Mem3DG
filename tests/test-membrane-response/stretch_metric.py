"""How needle-like is the deformation?  Reports edge stretch vs the initial mesh."""
import glob,math,sys
def rd(p):
    V=[];F=[]
    for l in open(p):
        if l.startswith('v '): V.append(tuple(float(x) for x in l.split()[1:4]))
        elif l.startswith('f '): F.append([int(t.split('/')[0])-1 for t in l.split()[1:]])
    return V,F
def edges(F):
    E=set()
    for f in F:
        for i in range(3):
            a,b=f[i],f[(i+1)%3]; E.add((min(a,b),max(a,b)))
    return sorted(E)
d=sys.argv[1]
mems=sorted(glob.glob(d+'/membrane_*.obj'))
V0,F0=rd(mems[0]); E=edges(F0)
L0=[math.dist(V0[a],V0[b]) for a,b in E]
print('%-14s %8s %8s %8s %8s'%('frame','maxStretch','p99','meanStr','r_max'))
for i in [0,len(mems)//4,len(mems)//2,3*len(mems)//4,len(mems)-1]:
    V,_=rd(mems[i])
    st=sorted((math.dist(V[a],V[b])/l0) for (a,b),l0 in zip(E,L0))
    rmax=max(math.sqrt(sum(c*c for c in q)) for q in V)
    print('%-14d %8.2f %8.2f %8.3f %8.3f'%(i,st[-1],st[int(0.99*len(st))],sum(st)/len(st),rmax))
