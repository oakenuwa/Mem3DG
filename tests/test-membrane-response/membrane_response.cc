// =============================================================================
// membrane_response.cc
//
// Does a local deformation propagate through the membrane?
//
// Pull ONE vertex outward, then read the membrane's own internal forcing and
// ask: do the pulled vertex's NEIGHBOURS feel anything?  A real membrane
// resists locally (bending) and transmits the strain outward (tension), so the
// 1-ring, 2-ring, ... should all carry force that decays smoothly with
// topological distance.  If only the pulled vertex responds, the surface is
// behaving as a bag of independent points and any protrusion will come out as a
// needle with hugely stretched triangles rather than a smooth bulge.
// =============================================================================
#include <mem3dg/capi/mem3dg_c.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <queue>
#include <map>
#include <set>

static const char* MESH = "/Users/hakenuwa/Desktop/Mem3DG/tests/icosphere.obj";

struct Case { const char* name; double Kb, Ksg, areaFrac, Kv, Kse; };

int main() {
  // topological rings around vertex 0, from the connectivity
  mem3dg_system* probe = mem3dg_system_create_from_file(MESH,1e-5,1e-5,0,2.5e-5);
  size_t nV = mem3dg_system_num_vertices(probe);
  size_t nF = mem3dg_system_num_faces(probe);
  std::vector<size_t> FI(nF*3); mem3dg_system_get_faces(probe,&FI[0]);
  std::vector<double> P0(nV*3); mem3dg_system_get_vertex_positions(probe,&P0[0]);
  mem3dg_system_destroy(probe);

  std::vector<std::set<size_t>> adj(nV);
  for (size_t f=0; f<nF; ++f)
    for (int i=0;i<3;++i){ size_t a=FI[3*f+i], b=FI[3*f+(i+1)%3];
                           adj[a].insert(b); adj[b].insert(a); }
  std::vector<int> ring(nV,-1); ring[0]=0;
  std::queue<size_t> q; q.push(0);
  while(!q.empty()){ size_t v=q.front(); q.pop();
    for(size_t w: adj[v]) if(ring[w]<0){ ring[w]=ring[v]+1; q.push(w); } }

  const double PULL = 0.15;   // how far vertex 0 is dragged outward

  Case cases[] = {
    {"bending only",                     8.22e-6, 0,     1.00, 0,    0},
    {"+ tension Ksg=100, At=A0",         8.22e-6, 100,   1.00, 0,    0},
    {"+ tension Ksg=100, PRE-TENSED .98",8.22e-6, 100,   0.98, 0,    0},
    {"+ Ksg=100 pre + Kv=100",           8.22e-6, 100,   0.98, 100,  0},
    {"+ ... + EDGE SPRINGS Kse=0.1",     8.22e-6, 100,   0.98, 100,  0.1},
    {"+ ... + EDGE SPRINGS Kse=1.0",     8.22e-6, 100,   0.98, 100,  1.0},
  };

  std::printf("Pull vertex 0 outward by %.2f, then read the membrane's own forces.\n",PULL);
  std::printf("Reported: mean |force| on each topological ring around vertex 0.\n");
  std::printf("A membrane that transmits strain shows force on rings 1,2,3 - not only ring 0.\n\n");

  for (const Case& c : cases) {
    mem3dg_system* m = mem3dg_system_create_from_file(MESH,c.Kb,c.Kb,0,2.5e-5);
    if (c.Ksg>0) mem3dg_system_set_tension(m,c.Ksg,c.areaFrac*mem3dg_system_area(m));
    if (c.Kv >0) mem3dg_system_set_osmotic(m,c.Kv,-1);
    if (c.Kse>0) mem3dg_system_set_regularization(m,c.Kse,0,0);

    // impose the deformation
    std::vector<double> P = P0;
    double r = std::sqrt(P0[0]*P0[0]+P0[1]*P0[1]+P0[2]*P0[2]);
    for (int j=0;j<3;++j) P[j] = P0[j]*(1.0+PULL/r);
    mem3dg_system_set_vertex_positions(m,&P[0]);

    std::vector<double> F(nV*3);
    mem3dg_system_get_mechanical_forces(m,&F[0]);

    std::map<int,std::pair<double,int>> byRing;
    for (size_t v=0; v<nV; ++v) {
      double mg = std::sqrt(F[3*v]*F[3*v]+F[3*v+1]*F[3*v+1]+F[3*v+2]*F[3*v+2]);
      auto& e = byRing[ring[v]]; e.first += mg; e.second += 1;
    }
    std::printf("%-36s", c.name);
    for (int rg=0; rg<=4; ++rg) {
      auto it = byRing.find(rg);
      if (it==byRing.end()) { std::printf("      -    "); continue; }
      std::printf(" %9.3e", it->second.first/it->second.second);
    }
    // ratio ring1/ring0 tells how much is transmitted
    double r0=byRing[0].first/byRing[0].second, r1=byRing[1].first/byRing[1].second;
    std::printf("   ring1/ring0 = %.4f\n", r0>0? r1/r0 : 0.0);
    mem3dg_system_destroy(m);
  }
  std::printf("\n%-36s %9s %9s %9s %9s %9s\n","","ring0","ring1","ring2","ring3","ring4");
  return 0;
}
