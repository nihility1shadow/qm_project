"""Compare compressed/fallback CSR to the v1.13 release, including exact signs."""
from pathlib import Path
import argparse,subprocess,tempfile
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--cxx',default='g++');args=parser.parse_args()
old=subprocess.check_output(['git','show','v1.13-two-buffer-cache:ahm-mb-sep.cpp'],cwd=root,text=True)
new=(root/'ahm-mb-sep.cpp').read_text()
def block(s):return s[s.index('struct SepmbReferenceTransitions {'):s.index('void sepmb_fock_rhs(')]
prefix='''#include <complex>
#include <vector>
#include <cstring>
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <stdexcept>
using namespace std;
using dcomplex=complex<double>;
const dcomplex I(0,1);
'''
main=r'''
bool same_bits(double a,double b) {return memcmp(&a,&b,sizeof(double))==0;}
int main() {
  mt19937_64 rng(114);uniform_real_distribution<double> sample(-1,1);
  int graphs=0,checks=0;size_t decoded=0;double maxerror=0;
  for(int ns : {1,2,7,31,193}) for(int mode=0;mode<9;mode++) {
    vector<vector<pair<int,double>>> edges(ns);
    for(int source=0;source<ns;source++) for(int target=0;target<ns;target++) {
      if(mode==0 || (source+3*target)%5!=0) continue;
      double coupling=.03125;
      if(mode==1) coupling=0.0;
      if(mode==3) coupling=sample(rng);
      if(mode==4 && (source+target)%2) coupling=nextafter(coupling,1.0);
      if(mode==5) coupling=1.e-310;
      if(mode==6) coupling=3.7e-200;
      if(mode==7) coupling=3.7e200;
      if(mode==8) coupling=0.0;
      if((source+target)%2 || mode==8) coupling=-coupling;
      edges[source].emplace_back(target,coupling);
      if(target%3==0) edges[source].emplace_back(target,-coupling); // duplicate edges
    }
    auto visit=[&](int source,const auto& callback) {
      for(auto e:edges[source]) callback(e.first,e.second);
    };
    auto a=previous::sepmb_build_reference_transitions(ns,visit);
    auto b=candidate::sepmb_build_reference_transitions(ns,visit);
    auto fallback=candidate::sepmb_build_reference_transitions(ns,visit,false);
    if(fallback.uniform_signed || a.offsets!=b.offsets || a.offsets!=fallback.offsets)
      throw runtime_error("offset/fallback mismatch");
    if(b.uniform_signed && !b.couplings.empty()) throw runtime_error("extra coupling storage");
    if(a.sources.empty() && b.uniform_signed) throw runtime_error("empty graph uniform flag");
    for(size_t k=0;k<a.sources.size();k++) {
      int encoded=b.sources[k];int source=b.uniform_signed && encoded<0?-(encoded+1):encoded;
      double coupling=b.uniform_signed?(encoded<0?-b.common_magnitude:b.common_magnitude):b.couplings[k];
      if(source!=a.sources[k] || !same_bits(coupling,a.couplings[k]) ||
          fallback.sources[k]!=a.sources[k] || !same_bits(fallback.couplings[k],a.couplings[k]))
        throw runtime_error("decoded signed source/coupling mismatch");
      decoded++;
    }
    previous::SepmbReferenceContext ca;candidate::SepmbReferenceContext cb,cf;
    ca.transitions=&a;cb.transitions=&b;cf.transitions=&fallback;
    vector<dcomplex> state(ns);
    for(int trial=0;trial<7;trial++) {
      for(int i=0;i<ns;i++) {
        state[i]=dcomplex(sample(rng),sample(rng));
        if(mode==7) state[i]*=1.e-100;
        if(trial==1) state[i]=0;
        if(trial==2) state[i]*=1.e-310;
        if(trial==3) state[i]=dcomplex(-0.0,0.0);
      }
      for(int target=0;target<ns;target++) {
        auto x=previous::sepmb_reference_hopping(ca,state.data(),target);
        auto y=candidate::sepmb_reference_hopping(cb,state.data(),target);
        auto z=candidate::sepmb_reference_hopping(cf,state.data(),target);
        maxerror=max(maxerror,abs(x-y));
        if(!same_bits(real(x),real(y)) || !same_bits(imag(x),imag(y)) ||
            !same_bits(real(x),real(z)) || !same_bits(imag(x),imag(z)))
          throw runtime_error("hopping bitwise mismatch");
        checks++;
      }
    }
    graphs++;
  }
  // A nextafter-sized difference must select the real-valued fallback.
  auto nonuniform=candidate::sepmb_build_reference_transitions(2,[](int source,const auto& v) {
    v(0,source?nextafter(.25,1.0):.25);
  });
  if(nonuniform.uniform_signed) throw runtime_error("rounded nonuniform coupling");
  cout<<graphs<<" graph cases, "<<decoded<<" decoded edges, "<<checks
      <<" bitwise hopping checks, max error "<<maxerror<<" passed\n";
}
'''
with tempfile.TemporaryDirectory(prefix='qm-signed-csr-') as directory:
    directory=Path(directory);cpp=directory/'test.cpp';exe=directory/'test.exe'
    cpp.write_text(prefix+'namespace previous {\n'+block(old)+'\n}\nnamespace candidate {\n'+block(new)+'\n}\n'+main)
    subprocess.run([args.cxx,'-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True,timeout=120)
    subprocess.run([str(exe)],check=True,timeout=120)
