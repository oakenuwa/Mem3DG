// Verify Cytosim placement/confinement against the Mem3DG membrane.
#include <mem3dg/capi/mem3dg_c.h>
#include "simul.h"
#include "parser.h"
#include "single.h"
#include "single_set.h"
#include "couple.h"
#include "couple_set.h"
#include "spaces/space_mem3dg.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <string>
static double surfDist(mem3dg_system* m, Vector p, int* insideOut){
  double pt[3];
  mem3dg_system_project(m,p.XX,p.YY,p.ZZ,pt,0);
  *insideOut = mem3dg_system_inside(m,p.XX,p.YY,p.ZZ);
  double dx=p.XX-pt[0],dy=p.YY-pt[1],dz=p.ZZ-pt[2];
  return std::sqrt(dx*dx+dy*dy+dz*dz);
}
static void report(const char* label, std::vector<double>& d, std::vector<int>& in){
  if(d.empty()){ std::printf("  %-22s (none)\n",label); return; }
  std::sort(d.begin(),d.end());
  double mean=0; for(double x:d) mean+=x; mean/=d.size();
  int nin=0; for(int x:in) nin+= (x==1);
  std::printf("  %-22s n=%-4zu  dist mean %.3e  max %.3e   inside %d/%zu\n",
              label,d.size(),mean,d.back(),nin,in.size());
}
int main(int argc,char**argv){
  mem3dg_system* m=mem3dg_system_create_from_file(argv[2],8.22e-6,8.22e-6,0,2.5e-5);
  SpaceMem3DG::setHandle(m);
  Simul simul; simul.initCytosim();
  Parser(&simul,1,1,1,0,0).readConfig(argv[1]);
  simul.prepare();
  std::printf("membrane R=0.8, %zu faces\n\n", mem3dg_system_num_faces(m));

  std::vector<double> dSurfA,dSurfD,dInner; std::vector<int> iSurfA,iSurfD,iInner;
  for(Single* s=simul.singles.firstF(); s; s=s->next()){
    int in; double d=surfDist(m,s->position(),&in);
    std::string n=s->prop->name();
    if(n=="surfaceAnchor"){dSurfA.push_back(d);iSurfA.push_back(in);}
    else if(n=="surfaceDiffuser"){dSurfD.push_back(d);iSurfD.push_back(in);}
    else {dInner.push_back(d);iInner.push_back(in);}
  }
  std::printf("SINGLES\n");
  report("surfaceAnchor",dSurfA,iSurfA);
  report("surfaceDiffuser",dSurfD,iSurfD);
  report("innerDiffuser",dInner,iInner);

  std::vector<double> cS,cI; std::vector<int> jS,jI;
  for(Couple* c=simul.couples.firstFF(); c; c=c->next()){
    int in; double d=surfDist(m,c->position(),&in);
    std::string n=c->prop->name();
    if(n=="surfaceCross"){cS.push_back(d);jS.push_back(in);} else {cI.push_back(d);jI.push_back(in);}
  }
  std::printf("\nCOUPLES\n");
  report("surfaceCross",cS,jS);
  report("innerCross",cI,jI);

  // --- now step, and check the confinement still holds -------------------
  int nSteps = (argc>3)? atoi(argv[3]) : 2000;
  for(int i=0;i<nSteps;++i){ simul.solve_meca(); simul.steps(); }
  std::printf("\n=== after %d steps ===\n", nSteps);

  dSurfA.clear();dSurfD.clear();dInner.clear();iSurfA.clear();iSurfD.clear();iInner.clear();
  for(Single* s=simul.singles.firstF(); s; s=s->next()){
    int in; double d=surfDist(m,s->position(),&in);
    std::string n=s->prop->name();
    if(n=="surfaceAnchor"){dSurfA.push_back(d);iSurfA.push_back(in);}
    else if(n=="surfaceDiffuser"){dSurfD.push_back(d);iSurfD.push_back(in);}
    else {dInner.push_back(d);iInner.push_back(in);}
  }
  std::printf("SINGLES\n");
  report("surfaceAnchor",dSurfA,iSurfA);
  report("surfaceDiffuser",dSurfD,iSurfD);
  report("innerDiffuser",dInner,iInner);
  cS.clear();cI.clear();jS.clear();jI.clear();
  for(Couple* c=simul.couples.firstFF(); c; c=c->next()){
    int in; double d=surfDist(m,c->position(),&in);
    std::string n=c->prop->name();
    if(n=="surfaceCross"){cS.push_back(d);jS.push_back(in);} else {cI.push_back(d);jI.push_back(in);}
  }
  std::printf("\nCOUPLES\n");
  report("surfaceCross",cS,jS);
  report("innerCross",cI,jI);

  mem3dg_system_destroy(m);
  return 0;
}
