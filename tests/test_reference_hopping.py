"""Validate real mask edges against creation/annihilation on ordered sets."""
from pathlib import Path
import argparse, subprocess, tempfile
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--cxx',default='g++');args=parser.parse_args()
source=(root/'ahm-mb-sep.cpp').read_text()
start=source.index('double sepmb_binom(')
stop=source.index('double sepmb_kondo_degeneracy(',start)
prefix='''#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <stdexcept>
using namespace std;
'''
main=r'''
void check_graph(int norb, const set<int>& initial, int distance, size_t& edge_checks) {
  auto states=sepmb_generate_reference_states(norb, initial, distance);
  auto masks=sepmb_generate_reference_masks(norb, initial, distance);
  map<set<int>,int> set_index;
  unordered_map<unsigned long long,int> mask_index;
  vector<double> coupling(norb);
  for(int i=0;i<norb;i++) coupling[i]=((i%3)-1.25)/(1.0+i);
  for(size_t i=0;i<states.size();i++) {set_index[states[i]]=i;mask_index[masks[i]]=i;}
  vector<map<int,double>> matrix(states.size());
  for(size_t source=0;source<states.size();source++) {
    vector<pair<int,double>> expected,actual;
    for(int orbital=1;orbital<norb;orbital++) {
      int annihilate=states[source].count(0)?0:orbital;
      int create=annihilate==0?orbital:0;
      auto target=states[source];
      auto a=target.find(annihilate);
      if(a==target.end()) continue;
      int sign=0;
      sign=(std::distance(target.begin(),a)%2)?-1:1;
      target.erase(a);
      if(target.count(create)) continue;
      if(std::distance(target.begin(),target.lower_bound(create))%2) sign=-sign;
      target.insert(create);
      auto found=set_index.find(target);
      if(found!=set_index.end()) expected.emplace_back(found->second,sign*coupling[orbital]);
    }
    sepmb_visit_mask_reference_edges(norb,masks[source],mask_index,coupling.data(),
      [&](int target,double value) {actual.emplace_back(target,value);matrix[source][target]=value;});
    if(actual!=expected) throw runtime_error("signed edge/order mismatch");
    edge_checks+=actual.size();
  }
  for(size_t source=0;source<states.size();source++)
    for(auto edge:matrix[source]) {
      auto reverse=matrix[edge.first].find(source);
      if(reverse==matrix[edge.first].end() || reverse->second!=edge.second)
        throw runtime_error("non-Hermitian reference graph");
    }
}
int main() {
  size_t checks=0;int graphs=0;
  for(int norb=3;norb<=10;norb++) for(int nel=1;nel<norb;nel++)
    for(int mode=0;mode<2;mode++) for(int depth=0;depth<=3;depth++) {
      set<int> initial;
      for(int i=0;i<nel;i++) initial.insert(mode?norb-1-i:i);
      check_graph(norb,initial,depth,checks);graphs++;
    }
  for(const set<int>& initial : {set<int>{0,63},set<int>{62,63}}) {
    check_graph(64,initial,1,checks);graphs++;
  }
  cout<<graphs<<" signed reference graphs, "<<checks<<" directed edges passed\n";
}
'''
with tempfile.TemporaryDirectory(prefix='qm-hopping-') as directory:
    directory=Path(directory);cpp=directory/'test.cpp';exe=directory/'test.exe'
    cpp.write_text(prefix+source[start:stop]+main)
    subprocess.run([args.cxx,'-O2','-std=c++17',str(cpp),'-o',str(exe)],check=True,timeout=120)
    subprocess.run([str(exe)],check=True,timeout=120)
