"""Compare every reference operator to the tagged v1.12 implementation."""
from pathlib import Path
import argparse,subprocess,tempfile
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--cxx',default='g++');args=parser.parse_args()
old=subprocess.check_output(['git','show','v1.12-compact-reference:ahm-mb-sep.cpp'],cwd=root,text=True)
new=(root/'ahm-mb-sep.cpp').read_text()
def block(source,marker):
    return source[source.index(marker):source.index('void sepmb_fock_rhs(')]
prefix='''#include <complex>
#include <vector>
#include <cstring>
#include <cmath>
#include <random>
#include <iostream>
#include <algorithm>
using namespace std;
using dcomplex=complex<double>;
const dcomplex I(0,1);
#define bzero(p,n) memset(p,0,n)
'''
main=r'''
template<class Context>
void config(Context& c,int ns,int nf,double* energy,int* occupied,dcomplex* work) {
  c.nstate=ns;c.nfock=nf;c.fock_begin=0;c.total_fock=nf;c.distributed=false;
  c.electronic_energy=energy;c.molecule_occupied=occupied;c.frequency=.037;
  c.half_displacement=1.73;c.constant_shift=.019;c.energy_origin=-.124;
}
int main() {
  mt19937_64 rng(113);uniform_real_distribution<double> sample(-1,1);
  double maximum=0;int checks=0;
  for(int ns : {1,7,31}) for(int nf : {2,7,16}) {
    int dim=ns*nf;
    vector<vector<pair<int,dcomplex>>> old_graph(ns);
    candidate::SepmbReferenceTransitions graph;
    for(int s=0;s<ns;s++) {
      for(int target=0;target<ns;target++) if(target!=s && (s+target)%3==0) {
        double weight=sample(rng)*.03;
        old_graph[s].emplace_back(target,weight);
      }
    }
    graph=candidate::sepmb_build_reference_transitions(ns,[&](int source,const auto& visit) {
      for(auto edge:old_graph[source]) visit(edge.first,real(edge.second));
    });
    vector<double> energy(ns);vector<int> occupied(ns);
    for(int s=0;s<ns;s++) {energy[s]=sample(rng);occupied[s]=s%2;}
    vector<dcomplex> state(dim),a(dim),b(dim),wa(dim),wb(dim);
    previous::SepmbReferenceContext ca;candidate::SepmbReferenceContext cb;
    config(ca,ns,nf,energy.data(),occupied.data(),wa.data());
    config(cb,ns,nf,energy.data(),occupied.data(),wb.data());
    ca.hamiltonian_work=wa.data();
    ca.transitions=&old_graph;cb.transitions=&graph;
    for(int trial=0;trial<5;trial++) for(int part=0;part<5;part++) {
      for(int i=0;i<dim;i++) {
        state[i]=dcomplex(sample(rng),sample(rng));
        if(i%9==0)state[i]=0;
        if(i%11==0)state[i]*=1e-310;
        a[i]=b[i]=dcomplex(sample(rng),sample(rng));
      }
      vector<dcomplex> initial_derivative=a;
      ca.operator_part=cb.operator_part=part;
      double af=sample(rng),df=sample(rng);
      previous::sepmb_reference_rhs(dim,0,af,df,&ca,state.data(),a.data());
      candidate::sepmb_reference_rhs(dim,0,af,df,&cb,state.data(),b.data());
      for(int i=0;i<dim;i++)maximum=max(maximum,abs(a[i]-b[i]));
      checks++;
      if(nf>=4) {
        for(int half=0;half<2;half++) {
          int begin=half?nf/2:0,end=half?nf:nf/2;
          cb.nfock=end-begin;cb.fock_begin=begin;cb.halo.assign(4*ns,0.0);
          for(int h=0;h<4;h++) {
            int row=h<2?begin-2+h:end+h-2;
            if(row>=0 && row<nf)
              copy(state.begin()+row*ns,state.begin()+(row+1)*ns,cb.halo.begin()+h*ns);
          }
          vector<dcomplex> local(initial_derivative.begin()+begin*ns,initial_derivative.begin()+end*ns);
          candidate::sepmb_reference_rhs((end-begin)*ns,0,af,df,&cb,state.data()+begin*ns,local.data());
          for(size_t i=0;i<local.size();i++)maximum=max(maximum,abs(local[i]-a[begin*ns+i]));
          checks++;
        }
        cb.nfock=nf;cb.fock_begin=0;
      }
    }
  }
  cout<<checks<<" reference operator comparisons; max difference="<<maximum<<"\n";
  return maximum>2e-15?1:0;
}
'''
with tempfile.TemporaryDirectory(prefix='qm-rhs-') as directory:
    directory=Path(directory);cpp=directory/'test.cpp';exe=directory/'test.exe'
    cpp.write_text(prefix+'namespace previous {\n'+block(old,'struct SepmbReferenceContext {')+'}\nnamespace candidate {\n'+block(new,'struct SepmbReferenceTransitions {')+'}\n'+main)
    subprocess.run([args.cxx,'-O2','-std=c++17',str(cpp),'-o',str(exe)],check=True,timeout=120)
    subprocess.run([str(exe)],check=True,timeout=120)
