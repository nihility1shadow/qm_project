#include "./Johnson-path-sampler.h"

#include <cstdlib>
#include <iostream>
#include <iterator>

double *log_table = NULL;

namespace {

int choose_from_set(const set<int>& s) {
  if(s.empty()) return -1;
  int idx = static_cast<int>(drand48()*s.size());
  if(idx >= static_cast<int>(s.size())) idx = static_cast<int>(s.size()) - 1;
  auto it = s.begin();
  advance(it, idx);
  return *it;
}

void split_by_target(const int N, const set<int>& current, const set<int>& target,
    set<int>& both, set<int>& out_only, set<int>& in_only, set<int>& neither) {
  set_intersection(current.begin(), current.end(), target.begin(), target.end(),
                   inserter(both, both.begin()));
  set_difference(current.begin(), current.end(), target.begin(), target.end(),
                 inserter(out_only, out_only.begin()));
  set_difference(target.begin(), target.end(), current.begin(), current.end(),
                 inserter(in_only, in_only.begin()));

  set<int> current_or_target;
  set_union(current.begin(), current.end(), target.begin(), target.end(),
            inserter(current_or_target, current_or_target.begin()));
  for(int j=1; j<=N; j++) {
    if(current_or_target.find(j) == current_or_target.end()) neither.insert(j);
  }
}

}

void JohnsonPathSampler::compute_dp_table() {
  Ptd = array2d<double>(Kmax+1, Dmax+1);
  Ptd[0][0] = 1.0;
  if(M == 0 || M == N) return;

  const double denom = M*(N-M);
  for(int t=1; t<=Kmax; t++) {
    for(int d=0; d<=Dmax; d++) {
      double p = 0.0;
      if(d > 0) {
        p += d*d/denom*Ptd[t-1][d-1];
      }
      p += d*(N-2*d)/denom*Ptd[t-1][d];
      if(d < Dmax) {
        p += (M-d)*(N-M-d)/denom*Ptd[t-1][d+1];
      }
      Ptd[t][d] = p;
    }
  }

  return;
}

int JohnsonPathSampler::sample_path(const int k, const set<int>& S0, const set<int>& S1,
    vector<pair<int, int>> &path) const {
  path.clear();
  if(k < 0 || k > Kmax) return -1;
  if(static_cast<int>(S0.size()) != M || static_cast<int>(S1.size()) != M) return -2;

  set<int> diff;
  set_symmetric_difference(S0.begin(), S0.end(), S1.begin(), S1.end(),
                           inserter(diff, diff.begin()));
  int distance = static_cast<int>(diff.size())/2;
  if(distance > Dmax || Ptd[k][distance] <= 0.0) return -3;

  set<int> current = S0;
  for(int remaining=k; remaining>0; remaining--) {
    set<int> both, out_only, in_only, neither;
    split_by_target(N, current, S1, both, out_only, in_only, neither);
    distance = static_cast<int>(out_only.size());

    double w_down = 0.0, w_same_keep = 0.0, w_same_miss = 0.0, w_up = 0.0;
    if(distance > 0) {
      w_down = out_only.size()*in_only.size()*Ptd[remaining-1][distance-1];
    }
    w_same_keep = both.size()*in_only.size()*Ptd[remaining-1][distance];
    w_same_miss = out_only.size()*neither.size()*Ptd[remaining-1][distance];
    if(distance < Dmax) {
      w_up = both.size()*neither.size()*Ptd[remaining-1][distance+1];
    }

    const double total = w_down + w_same_keep + w_same_miss + w_up;
    if(total <= 0.0) return -4;

    double r = drand48()*total;
    int out = -1, in = -1;
    if(r < w_down) {
      out = choose_from_set(out_only);
      in  = choose_from_set(in_only);
    } else if((r -= w_down) < w_same_keep) {
      out = choose_from_set(both);
      in  = choose_from_set(in_only);
    } else if((r -= w_same_keep) < w_same_miss) {
      out = choose_from_set(out_only);
      in  = choose_from_set(neither);
    } else {
      out = choose_from_set(both);
      in  = choose_from_set(neither);
    }

    current.erase(out);
    current.insert(in);
    path.push_back({out, in});
  }

  return current == S1 ? 0 : -5;
}

int get_random_element(const set<int>& s) {
  return choose_from_set(s);
}

void print_path(const vector<pair<int, int>>& path, const set<int>& S0) {
  set<int> current = S0;
  print_set("start: ", current, "\n");
  for(size_t j=0; j<path.size(); j++) {
    current.erase(path[j].first);
    current.insert(path[j].second);
    std::cout<<"step "<<j+1<<": "<<path[j].first<<" -> "<<path[j].second<<" ";
    print_set(NULL, current, "\n");
  }

  return;
}

bool verify_path(const set<int> &S0, const set<int> &S1,
    const vector<pair<int, int>>& path) {
  set<int> current = S0;
  const size_t size0 = S0.size();

  for(size_t j=0; j<path.size(); j++) {
    if(current.find(path[j].first) == current.end()) return false;
    if(current.find(path[j].second) != current.end()) return false;
    current.erase(path[j].first);
    current.insert(path[j].second);
    if(current.size() != size0) return false;
  }

  return current == S1;
}
