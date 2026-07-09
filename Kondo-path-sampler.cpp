#include "./Kondo-path-sampler.h"

namespace {

double kondo_binom(const int n, const int k) {
  if(k < 0 || k > n) return 0.0;
  int kk = k < n-k ? k : n-k;
  double r = 1.0;
  for(int j=1; j<=kk; j++) r *= (n-kk+j)*1.0/j;
  return r;
}

int diff_size_no_zero(const set<int> &a, const set<int> &b) {
  int n = 0;
  for(int x : a) if(x != 0 && b.find(x) == b.end()) n++;
  return n;
}

}

/*
 * this function works now.
 */
void KondoPathSampler::compute_dp_table() {
  int Nvac = N - M;
  if(Ptrans) {
    free3d(Ptrans);
    Ptd = NULL;
    Qtd = NULL;
  }
  Ptrans = array3d<double>(2, Kmax+1, Dmax+1);
  Ptd = Ptrans[0];
  Qtd = Ptrans[1];

  Ptd[0][0] = 1.0;
  Qtd[0][0] = 1.0;

  if(M == 0 || Nvac == 0) return;

  const double inv_M = 1.0/M, inv_Nvac = 1.0/Nvac;

  for(int t=0; t<Kmax; t++) {
    if(t%2) {
      // P: A-start odd B-distance -> even A-distance.
      // B_d returns to A_d by removing a newly added orbital,
      // or to A_{d+1} by removing one of the original occupied orbitals.
      for(int d=0; d<=Dmax; d++) {
        if(d + 1 <= M) Ptd[t+1][d] += (d+1)*inv_M*Ptd[t][d];
        if(d > 0 && d <= M) Ptd[t+1][d] += (M-d)*inv_M*Ptd[t][d-1];

        // Q: B-start odd A-distance -> even B-distance.
        if(d + 1 <= Nvac) Qtd[t+1][d] += (d+1)*inv_Nvac*Qtd[t][d];
        if(d > 0 && d <= Nvac) Qtd[t+1][d] += (Nvac-d)*inv_Nvac*Qtd[t][d-1];
      }
    } else {
      // P: A-start even A-distance -> odd B-distance.
      for(int d=0; d<=Dmax; d++) {
        if(d <= Nvac) Ptd[t+1][d] += (Nvac-d)*inv_Nvac*Ptd[t][d];
        if(d < Dmax) Ptd[t+1][d] += (d+1)*inv_Nvac*Ptd[t][d+1];

        // Q: B-start even B-distance -> odd A-distance.
        if(d <= M) Qtd[t+1][d] += (M-d)*inv_M*Qtd[t][d];
        if(d < Dmax) Qtd[t+1][d] += (d+1)*inv_M*Qtd[t][d+1];
      }
    }
  }

  // Convert distance-class probabilities into probabilities for one concrete endpoint.
  // This is what the separated projection needs for a fixed target determinant.
  const int NA = M - 1;
  const int NBvac = N - M - 1;
  for(int t=0; t<=Kmax; t++) {
    for(int d=0; d<=Dmax; d++) {
      double pdeg = (t%2 == 0)
        ? kondo_binom(NA, d)*kondo_binom(Nvac, d)
        : kondo_binom(NA, d)*kondo_binom(Nvac, d+1);
      double qdeg = (t%2 == 0)
        ? kondo_binom(M, d)*kondo_binom(NBvac, d)
        : kondo_binom(M, d+1)*kondo_binom(NBvac, d);
      Ptd[t][d] = pdeg > 0.0 ? Ptd[t][d]/pdeg : 0.0;
      Qtd[t][d] = qdeg > 0.0 ? Qtd[t][d]/qdeg : 0.0;
    }
  }

  return;
}

//This function, coming from deepseek, works perfectly. 
//What it actually calculates is the sum of 
//probability of the nodes whose distance to the initial node is d after k moves.
#if 0
void KondoPathSampler::compute_dp_table_sum() {
  // 桶外元素总数
  int Nvac = N - M;
  double inv_NM = 1.0/Nvac, inv_M = 1.0/M;

  Ptrans = array3d<double>(2, Kmax+1, Dmax+1);
  Ptd = Ptrans[0];
  Qtd = Ptrans[1];

  // 边界条件: t=0 时 d=0 的概率为1
  Ptd[0][0] = 1.0;
  // 填充DP表
  for (int t = 1; t < Kmax; t++) {
    for (int d = 0; d <= Dmax; d++) {
      if (t%2) { // 奇数步
        if (d > 0)  Ptd[t][d] += (d+1)*inv_NM*Ptd[t-1][d+1];
                    Ptd[t][d] += (Nvac-d)*inv_NM*Ptd[t-1][d];
      } else { // 偶数步
        if (d < Dmax) Ptd[t][d] += (M-d+1)*inv_M*Ptd[t-1][d-1];
                      Ptd[t][d] += d*inv_M*Ptd[t-1][d];
      }
    }
  }

  // 边界条件: t=0 时 d=0 的概率为1
  Qtd[0][0] = 1.0;
  for (int n = 1; n <= Kmax; n++) {
      if (n % 2 == 1) { // 奇数步：到达A类
          for (int d = 1; d <= M; d++) {
              double term1 = d * Qtd[n-1][d];
              double term2 = (d > 0) ? (M - d + 1) * Qtd[n-1][d-1] : 0;
              Qtd[n][d] = (term1 + term2) / M;
          }
      } else { // 偶数步：到达B类
          for (int d = 0; d <= M; d++) {
              double term1 = (d < M) ? (d+1) * Qtd[n-1][d+1] : 0;
              double term2 = (N-M-d) * Qtd[n-1][d];
              Qtd[n][d] = (term1 + term2) / (N-M);
          }
      }
  }

  return;
}
#endif

void KondoPathSampler::compute_dp_table_sum() {
  compute_dp_table();
  return;
}

double KondoPathSampler::endpoint_probability(const set<int> &current,
    const set<int> &target, const int remaining) const {
  if(remaining < 0 || remaining > Kmax) return 0.0;
  if((int)current.size() != M || (int)target.size() != M) return 0.0;

  const bool current_A = current.find(0) != current.end();
  const bool target_A  = target.find(0)  != target.end();
  if(((remaining%2) == 0) != (current_A == target_A)) return 0.0;

  const int current_only = diff_size_no_zero(current, target);
  const int target_only  = diff_size_no_zero(target, current);
  int d = -1;

  if(current_A) {
    if(target_A) {
      if(current_only != target_only) return 0.0;
      d = current_only;
    } else {
      if(target_only != current_only + 1) return 0.0;
      d = current_only;
    }
    return (d >= 0 && d <= Dmax) ? Ptd[remaining][d] : 0.0;
  }

  if(target_A) {
    if(current_only != target_only + 1) return 0.0;
    d = target_only;
  } else {
    if(current_only != target_only) return 0.0;
    d = current_only;
  }
  return (d >= 0 && d <= Dmax) ? Qtd[remaining][d] : 0.0;
}

int KondoPathSampler::sample_path_direct(const int ksteps, const set<int> &S0,
    const set<int> &S1, vector<pair<int, int> > &path) const {
  path.clear();
  if(ksteps < 0 || ksteps > Kmax) return -1;
  if((int)S0.size() != M || (int)S1.size() != M) return -2;
  if(endpoint_probability(S0, S1, ksteps) <= 0.0) return -3;

  set<int> current = S0;
  for(int remaining=ksteps; remaining>0; remaining--) {
    struct Candidate { int first, second; double weight; };
    vector<Candidate> candidates;
    double total = 0.0;

    if(current.find(0) != current.end()) {
      for(int y=1; y<N; y++) {
        if(current.find(y) != current.end()) continue;
        set<int> next = current;
        next.erase(0);
        next.insert(y);
        double w = endpoint_probability(next, S1, remaining-1);
        if(w > 0.0) {
          candidates.push_back({0, y, w});
          total += w;
        }
      }
    } else {
      for(int x : current) {
        if(x == 0) continue;
        set<int> next = current;
        next.erase(x);
        next.insert(0);
        double w = endpoint_probability(next, S1, remaining-1);
        if(w > 0.0) {
          candidates.push_back({x, 0, w});
          total += w;
        }
      }
    }

    if(total <= 0.0 || candidates.empty()) return -4;
    double r = drand48()*total;
    const Candidate *choice = &candidates.back();
    for(const Candidate &cand : candidates) {
      if(r < cand.weight) {
        choice = &cand;
        break;
      }
      r -= cand.weight;
    }

    if(current.find(choice->first) == current.end() ||
       current.find(choice->second) != current.end()) return -6;
    current.erase(choice->first);
    current.insert(choice->second);
    path.push_back({choice->first, choice->second});
  }

  return current == S1 ? 0 : -5;
}

int KondoPathSampler::sample_path_A2A(const int ksteps, const set<int> &kS0, const set<int> &kS1,
        vector<pair<int, int> > &path) const {
  if(kS0.find(0) == kS0.end() || kS1.find(0) == kS1.end() || ksteps%2) return -1;
  return sample_path_direct(ksteps, kS0, kS1, path);
}

int KondoPathSampler::sample_path_A2B(const int ksteps, const set<int> &kS0, const set<int> &kS1,
        vector<pair<int, int> > &path) const {
  if(kS0.find(0) == kS0.end() || kS1.find(0) != kS1.end() || ksteps%2==0) return -1;
  return sample_path_direct(ksteps, kS0, kS1, path);
}
 
int KondoPathSampler::sample_path_B2A(const int ksteps, const set<int> &kS0, const set<int> &kS1,
        vector<pair<int, int> > &path) const {
  if(kS0.find(0) != kS0.end() || kS1.find(0) == kS1.end() || ksteps%2==0) return -1;
  return sample_path_direct(ksteps, kS0, kS1, path);
}

int KondoPathSampler::sample_path_B2B(const int ksteps, const set<int> &kS0, const set<int> &kS1,
        vector<pair<int, int> > &path) const {
  if(kS0.find(0) != kS0.end() || kS1.find(0) != kS1.end() || ksteps%2) return -1;
  return sample_path_direct(ksteps, kS0, kS1, path);
}

#if 1
int exchange_A2B_for(set<int> &Sk, 
    int &size_A, set<int> &A, 
    int &size_B, set<int> &B, 
    int &size_C, set<int> &C, 
    int &size_D, set<int> &D,
    multiset<int> &opts, vector<pair<int, int> > &path) {
    //at this moment, current set Sk is an A
    //we need choose one element from C or D and exchange it with 0 in Sk
    int idx = drand48()*(size_C + size_D);
    if(idx<size_C) {
       auto it = C.begin();
       advance(it, idx);
       idx = *it;
       auto pos = opts.find(*it);
       opts.erase(pos);
       C.erase(idx);
       size_C--;
       if(opts.find(idx)!=opts.end()) {
         A.insert(idx);
         size_A++;
       }
    } else {
       auto it = D.begin();
       advance(it, idx-size_C);
       auto pos = opts.find(*it);
       opts.erase(pos);
       idx = *it;
       D.erase(idx);
       size_D--;
       if(opts.find(idx)!=opts.end()) {
         B.insert(idx);
         size_B++;
       }
    }
#ifdef __ENABLE_CHECK_
    if(Sk.find(idx) != Sk.end()) {
      cerr<<"the element "<<idx<<" to be added is already in S."<<endl;
    }
#endif
    Sk.erase(0);
    Sk.insert(idx);
    path.push_back({0, idx});

    return idx;
}

int exchange_B2A_for(set<int> &Sk, 
    int &size_A, set<int> &A, 
    int &size_B, set<int> &B, 
    int &size_C, set<int> &C, 
    int &size_D, set<int> &D,
    multiset<int> &opts, vector<pair<int, int> > &path) {
    //Now, current set Sk is an B
    //we need choose one element from B or A, which is automatically in Sk,
    //to exchange with 0 outside
    int idx = drand48()*(size_A + size_B);
    if(idx<size_A) {
       auto it = A.begin();
       advance(it, idx);
       idx = *it;
       auto pos = opts.find(*it);
       opts.erase(pos);
       A.erase(idx);
       size_A--;
       if(opts.find(idx)!=opts.end()) {
         C.insert(idx);
         size_C++;
       }
    } else {
       auto it = B.begin();
       advance(it, idx-size_A);
       auto pos = opts.find(*it);
       opts.erase(pos);
       idx = *it;
       B.erase(idx);
       size_B--;
       if(opts.find(idx)!=opts.end()) {
         D.insert(idx);
         size_D++;
       }
    }
    Sk.erase(idx);
    Sk.insert(0);
    path.push_back({idx, 0});

    return idx;
}

int exchange_A2B_back(set<int> &Sk, 
    int &size_A, set<int> &A, 
    int &size_B, set<int> &B, 
    int &size_C, set<int> &C, 
    int &size_D, set<int> &D,
    multiset<int> &opts, vector<pair<int, int> > &path) {
    int idx = drand48()*(size_B + size_D);
    if(idx<size_B) {
       auto it = B.begin();
       advance(it, idx);
       idx = *it;
       auto pos = opts.find(*it);
       opts.erase(pos);
       B.erase(idx);
       size_B--;
       if(opts.find(idx)!=opts.end()) {
         A.insert(idx);
         size_A++;
       }
    } else {
       auto it = D.begin();
       advance(it, idx-size_B);
       auto pos = opts.find(*it);
       opts.erase(pos);
       idx = *it;
       D.erase(idx);
       size_D--;
       if(opts.find(idx)!=opts.end()) {
         C.insert(idx);
         size_C++;
       }
    }
#ifdef __ENABLE_CHECK_
    if(Sk.find(idx) != Sk.end()) {
      cerr<<"the element "<<idx<<" to be added is already in S."<<endl;
    }
#endif
    Sk.erase(0);
    Sk.insert(idx);
    path.push_back({0, idx});

    return idx;
}

int exchange_B2A_back(set<int> &Sk, 
    int &size_A, set<int> &A, 
    int &size_B, set<int> &B, 
    int &size_C, set<int> &C, 
    int &size_D, set<int> &D,
    multiset<int> &opts, vector<pair<int, int> > &path) {
    int idx = drand48()*(size_A + size_C);
    if(idx<size_A) {
       auto it = A.begin();
       advance(it, idx);
       idx = *it;
       auto pos = opts.find(*it);
       opts.erase(pos);
       A.erase(idx);
       size_A--;
       if(opts.find(idx)!=opts.end()) {
         B.insert(idx);
         size_B++;
       }
    } else {
       auto it = C.begin();
       advance(it, idx-size_A);
       auto pos = opts.find(*it);
       opts.erase(pos);
       idx = *it;
       C.erase(idx);
       size_C--;
       if(opts.find(idx)!=opts.end()) {
         D.insert(idx);
         size_D++;
       }
    }

    Sk.erase(idx);
    Sk.insert(0);
    path.push_back({idx, 0});

    return idx;
}

/*
 * S0 and S1 should be A nodes and ksteps should be even.
 * forward and backward search for the path
 */
int KondoPathSampler::sample_path_A2A_fb(const int ksteps, const set<int> &S0, 
    const set<int> &S1, vector<pair<int, int> > &path) const {
#ifdef __ENABLE_CHECK_
   if(S0.find(0) == S0.end() || 
      S1.find(0) == S1.end() || 
      ksteps%2) {
     cerr<<"KondoPathSampler::sample_path_A2A_uniform : \n";
     cerr<<"    S0 and S1 should be A nodes and ksteps should be even"<<endl;
     abort();
   }
#endif 
  vector<pair<int, int> > path_bck;
  path.clear();
  // 初始化四个集合
  set<int> A, B, C, D;
  multiset<int> opts;
  // A =  S0 ∩ S1  ∩ opts
  // B = (S0 - S1) ∩ opts : elements needed to move out
  // C = (S1 - S0) ∩ opts : elments needed to move in
  // D = opts - S0 - S1
  //     those elemnts needed for exchange lie neither in S0 nor in S1
  set_difference(S0.begin(), S0.end(), S1.begin(), S1.end(), 
                 inserter(B, B.begin()));
  set_difference(S1.begin(), S1.end(), S0.begin(), S0.end(), 
                 inserter(C, C.begin()));
  int current_d = C.size();

  {
    set<int> reduced_opts;
    set_union(B.begin(), B.end(), C.begin(), C.end(), 
            inserter(reduced_opts, reduced_opts.begin()));
    for(int j : reduced_opts) opts.insert(j);
  }

#ifdef __ENABLE_CHECK_
  if(ksteps<opts.size()) {
     cerr<<"KondoPathSampler::sample_path_A2A_uniform : \n";
     cerr<<"    ksteps should be at least 2*d\n";
     cerr<<"    But now kstep = "<<ksteps<<" d = "<<current_d<<endl;
     abort();
  }
#endif
  const int hsteps = ksteps/2;
  for(int j=0; j<hsteps-current_d; j++) {
    int idx = drand48()*(N-1)+1;
    opts.insert(idx);
    opts.insert(idx); //duplicte the element
    //D = opts - S0 - S1
    //A = S0 ∩ S1 ∩ opts
    bool in_S0 = S0.find(idx) != S0.end(),
         in_S1 = S1.find(idx) != S1.end();
    if((!in_S0) && (!in_S1)) D.insert(idx); 
    if(  in_S0  &&   in_S1 ) A.insert(idx); 
  }

#ifdef __ENABLE_CHECK_
  cout<<"Before exchange: \n";
  print_set("S0: ",S0, "\n");
  print_set("S1: ",S1, "\n");
  print_set("A : ", A, "\n");
  print_set("B : ", B, "\n");
  print_set("C : ", C, "\n");
  print_set("D : ", D, "\n");
  print_set("op: ", opts, "\n");
#endif

  int size_A = A.size(),
      size_B = B.size(),
      size_C = C.size(), 
      size_D = D.size();
  set<int> Sk_for = S0,
           Sk_bck = S1;
  for(int k=0; k<hsteps/2; k++) {
    if(size_B + size_D==0) return -1;
    int idx = exchange_A2B_back(Sk_bck, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path_bck);
#ifdef __ENABLE_CHECK_
    cout<<"The "<<ksteps-2*k<<" exchange: 0 -> "<<idx<<"\n";
    print_set("    S0: ",Sk_for, "\n");
    print_set("    S1: ",Sk_bck, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif

    if(size_C + size_D==0) return -1;
    idx = exchange_A2B_for(Sk_for, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path);
#ifdef __ENABLE_CHECK_
    cout<<"The "<<2*k+1<<" exchange: 0 -> "<<idx<<"\n";
    print_set("    S0: ",Sk_for, "\n");
    print_set("    S1: ",Sk_bck, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif

    if(size_A + size_B==0) return -1;
    idx = exchange_B2A_for(Sk_for, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path);
#ifdef __ENABLE_CHECK_
    cout<<"The "<<2*k+2<<" exchange: "<<idx<<" -> 0\n";
    print_set("    S0: ",Sk_for, "\n");
    print_set("    S1: ",Sk_bck, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif

    if(size_A + size_C==0) return -1;
    idx = exchange_B2A_back(Sk_bck, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path_bck);
#ifdef __ENABLE_CHECK_
    cout<<"The "<<ksteps-2*k-1<<" exchange: "<<idx<<" -> 0\n";
    print_set("    S0: ",Sk_for, "\n");
    print_set("    S1: ",Sk_bck, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif
  }

  for(int k=0; k<hsteps%2; k++) {
    if(size_B + size_D==0) return -1;
    int idx = exchange_A2B_back(Sk_bck, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path_bck);
#ifdef __ENABLE_CHECK_
    cout<<"The "<<hsteps+2<<" exchange: 0 -> "<<idx<<"\n";
    print_set("    S0: ",Sk_for, "\n");
    print_set("    S1: ",Sk_bck, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif

    if(size_C + size_D==0) return -1;
    idx = exchange_A2B_for(Sk_for, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path);
#ifdef __ENABLE_CHECK_
    cout<<"The "<<hsteps+1<<" exchange: 0 -> "<<idx<<"\n";
    print_set("    S0: ",Sk_for, "\n");
    print_set("    S1: ",Sk_bck, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif
  }

  for(int k=hsteps-1; k>-1; k--) {
    path.push_back({path_bck[k].second, path_bck[k].first});
  }

  return 0;
}

/*
 * S0 and S1 should be A nodes and ksteps should be even
 * This function does not work at this moment.
 */
int KondoPathSampler::sample_path_A2A_uniform(const int ksteps, const set<int> &S0, 
    const set<int> &S1, vector<pair<int, int> > &path) const {
#ifdef __ENABLE_CHECK_
   if(S0.find(0) == S0.end() || 
      S1.find(0) == S1.end() || 
      ksteps%2) {
     cerr<<"KondoPathSampler::sample_path_A2A_uniform : \n";
     cerr<<"    S0 and S1 should be A nodes and ksteps should be even"<<endl;
     abort();
   }
#endif 

  path.clear();
  // 初始化四个集合
  set<int> A, B, C, D;
  multiset<int> opts;
  // A =  S0 ∩ S1  ∩ opts
  // B = (S0 - S1) ∩ opts : elements needed to move out
  // C = (S1 - S0) ∩ opts : elments needed to move in
  // D = opts - S0 - S1
  //     those elemnts needed for exchange lie neither in S0 nor in S1
  set_difference(S0.begin(), S0.end(), S1.begin(), S1.end(), 
                 inserter(B, B.begin()));
  set_difference(S1.begin(), S1.end(), S0.begin(), S0.end(), 
                 inserter(C, C.begin()));
  int current_d = C.size();

  {
    set<int> reduced_opts;
    set_union(B.begin(), B.end(), C.begin(), C.end(), 
            inserter(reduced_opts, reduced_opts.begin()));
    for(int j : reduced_opts) opts.insert(j);
  }

  print_set("op = B+C: ", opts, "\n");
#ifdef __ENABLE_CHECK_
  if(ksteps<opts.size()) {
     cerr<<"KondoPathSampler::sample_path_A2A_uniform : \n";
     cerr<<"    ksteps should be at least 2*d\n";
     cerr<<"    But now kstep = "<<ksteps<<" d = "<<current_d<<endl;
     abort();
  }
#endif
  const int hsteps = ksteps/2;
  for(int j=0; j<hsteps-current_d; j++) {
    int idx = drand48()*(N-1)+1;
    opts.insert(idx);
    print_set("op 1 : ", opts, "\n");
    opts.insert(idx); //duplicte the element
    print_set("op 2 : ", opts, "\n");
    //D = opts - S0 - S1
    //A = S0 ∩ S1 ∩ opts
    bool in_S0 = S0.find(idx) != S0.end(),
         in_S1 = S1.find(idx) != S1.end();
    if((!in_S0) && (!in_S1)) D.insert(idx); 
    if(  in_S0  &&   in_S1 ) A.insert(idx); 
  }

#ifdef __ENABLE_CHECK_
  cout<<"Before exchange: \n";
  print_set("S0: ",S0, "\n");
  print_set("S1: ",S1, "\n");
  print_set("A : ", A, "\n");
  print_set("B : ", B, "\n");
  print_set("C : ", C, "\n");
  print_set("D : ", D, "\n");
  print_set("op: ", opts, "\n");
#endif

  int size_A = A.size(),
      size_B = B.size(),
      size_C = C.size(), 
      size_D = D.size();
  set<int> Sk = S0;
  //forward path
  for(int k=0; k<hsteps; k++) {
    if(size_C + size_D==0) return -1;
#if 0
    //at this moment, current set Sk is an A
    //we need choose one element from C or D and exchange it with 0 in Sk
    int idx = drand48()*(size_C + size_D);
    if(idx<size_C) {
       auto it = C.begin();
       advance(it, idx);
       idx = *it;
       auto pos = opts.find(*it);
       opts.erase(pos);
       C.erase(idx);
       size_C--;
       if(opts.find(idx)!=opts.end()) {
         A.insert(idx);
         size_A++;
       }
    } else {
       auto it = D.begin();
       advance(it, idx-size_C);
       auto pos = opts.find(*it);
       opts.erase(pos);
       idx = *it;
       D.erase(idx);
       size_D--;
       if(opts.find(idx)!=opts.end()) {
         B.insert(idx);
         size_B++;
       }
    }
#ifdef __ENABLE_CHECK_
    if(Sk.find(idx) != Sk.end()) {
      cerr<<"the element "<<idx<<" to be added is already in S."<<endl;
    }
#endif
    Sk.erase(0);
    Sk.insert(idx);
    path.push_back({0, idx});
#else
    int idx = exchange_A2B_for(Sk, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path);
#endif

#ifdef __ENABLE_CHECK_
    cout<<"The "<<2*k<<" exchange: 0 -> "<<idx<<"\n";
    print_set("    Sk: ",Sk, "\n");
    print_set("    S1: ",S1, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif

    //Now, current set Sk is an B
    //we need choose one element from B or A, which is automatically in Sk,
    //to exchange with 0 outside
    if(size_A + size_B==0) return -1;
#if 0
    idx = drand48()*(size_A + size_B);
    if(idx<size_A) {
       auto it = A.begin();
       advance(it, idx);
       idx = *it;
       auto pos = opts.find(*it);
       opts.erase(pos);
       A.erase(idx);
       size_A--;
       if(opts.find(idx)!=opts.end()) {
         C.insert(idx);
         size_C++;
       }
    } else {
       auto it = B.begin();
       advance(it, idx-size_A);
       auto pos = opts.find(*it);
       opts.erase(pos);
       idx = *it;
       B.erase(idx);
       size_B--;
       if(opts.find(idx)!=opts.end()) {
         D.insert(idx);
         size_D++;
       }
    }
    Sk.erase(idx);
    Sk.insert(0);
    path.push_back({idx, 0});
#else
    idx = exchange_B2A_for(Sk, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path);
#endif

#ifdef __ENABLE_CHECK_
    cout<<"The "<<2*k+1<<" exchange: "<<idx<<" -> 0\n";
    print_set("    Sk: ",Sk, "\n");
    print_set("    S1: ",S1, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif
  }

  return 0;
}

/*
 * S0 and S1 should be A nodes and ksteps should be even
 * This function does not work at this moment.
 */
int KondoPathSampler::sample_path_A2A_backward(const int ksteps, const set<int> &S0, 
    const set<int> &S1, vector<pair<int, int> > &forpath) const {
#ifdef __ENABLE_CHECK_
   if(S0.find(0) == S0.end() || 
      S1.find(0) == S1.end() || 
      ksteps%2) {
     cerr<<"KondoPathSampler::sample_path_A2A_uniform : \n";
     cerr<<"    S0 and S1 should be A nodes and ksteps should be even"<<endl;
     abort();
   }
#endif 

  forpath.clear();
  // 初始化四个集合
  set<int> A, B, C, D;
  multiset<int> opts;
  // A =  S0 ∩ S1  ∩ opts
  // B = (S0 - S1) ∩ opts : elements needed to move out
  // C = (S1 - S0) ∩ opts : elments needed to move in
  // D = opts - S0 - S1
  //     those elemnts needed for exchange lie neither in S0 nor in S1
  set_difference(S0.begin(), S0.end(), S1.begin(), S1.end(), 
                 inserter(B, B.begin()));
  set_difference(S1.begin(), S1.end(), S0.begin(), S0.end(), 
                 inserter(C, C.begin()));
  int current_d = C.size();

  {
    set<int> reduced_opts;
    set_union(B.begin(), B.end(), C.begin(), C.end(), 
            inserter(reduced_opts, reduced_opts.begin()));
    for(int j : reduced_opts) opts.insert(j);
  }

  print_set("op = B+C: ", opts, "\n");
#ifdef __ENABLE_CHECK_
  if(ksteps<opts.size()) {
     cerr<<"KondoPathSampler::sample_path_A2A_uniform : \n";
     cerr<<"    ksteps should be at least 2*d\n";
     cerr<<"    But now kstep = "<<ksteps<<" d = "<<current_d<<endl;
     abort();
  }
#endif
  const int hsteps = ksteps/2;
  for(int j=0; j<hsteps-current_d; j++) {
    int idx = drand48()*(N-1)+1;
    opts.insert(idx);
    print_set("op 1 : ", opts, "\n");
    opts.insert(idx); //duplicte the element
    print_set("op 2 : ", opts, "\n");
    //D = opts - S0 - S1
    //A = S0 ∩ S1 ∩ opts
    bool in_S0 = S0.find(idx) != S0.end(),
         in_S1 = S1.find(idx) != S1.end();
    if((!in_S0) && (!in_S1)) D.insert(idx); 
    if(  in_S0  &&   in_S1 ) A.insert(idx); 
  }

#ifdef __ENABLE_CHECK_
  cout<<"Before exchange: \n";
  print_set("S0: ",S0, "\n");
  print_set("S1: ",S1, "\n");
  print_set("A : ", A, "\n");
  print_set("B : ", B, "\n");
  print_set("C : ", C, "\n");
  print_set("D : ", D, "\n");
  print_set("op: ", opts, "\n");
#endif

  int size_A = A.size(),
      size_B = B.size(),
      size_C = C.size(), 
      size_D = D.size();
  set<int> Sk = S1;
  vector<pair<int, int> > path;
  //backward path
  for(int k=0; k<hsteps; k++) {
    //at this moment, current set Sk is an A
    //we need choose one element from C or D and exchange it with 0 in Sk
    if(size_B + size_D==0) return -1;
#if 0
    int idx = drand48()*(size_B + size_D);
    if(idx<size_B) {
       auto it = B.begin();
       advance(it, idx);
       idx = *it;
       auto pos = opts.find(*it);
       opts.erase(pos);
       B.erase(idx);
       size_B--;
       if(opts.find(idx)!=opts.end()) {
         A.insert(idx);
         size_A++;
       }
    } else {
       auto it = D.begin();
       advance(it, idx-size_B);
       auto pos = opts.find(*it);
       opts.erase(pos);
       idx = *it;
       D.erase(idx);
       size_D--;
       if(opts.find(idx)!=opts.end()) {
         C.insert(idx);
         size_C++;
       }
    }
#ifdef __ENABLE_CHECK_
    if(Sk.find(idx) != Sk.end()) {
      cerr<<"the element "<<idx<<" to be added is already in S."<<endl;
    }
#endif
    Sk.erase(0);
    Sk.insert(idx);
    path.push_back({0, idx});
#else
    int idx = exchange_A2B_back(Sk, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path);
#endif

#ifdef __ENABLE_CHECK_
    cout<<"The "<<2*k<<" exchange: 0 -> "<<idx<<"\n";
    print_set("    S0: ",S0, "\n");
    print_set("    Sk: ",Sk, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif

    //Now, current set Sk is an B
    //we need choose one element from A or C, which is automatically in Sk,
    //to exchange with 0 outside
    if(size_A + size_C==0) return -1;
#if 0
    idx = drand48()*(size_A + size_C);
    if(idx<size_A) {
       auto it = A.begin();
       advance(it, idx);
       idx = *it;
       auto pos = opts.find(*it);
       opts.erase(pos);
       A.erase(idx);
       size_A--;
       if(opts.find(idx)!=opts.end()) {
         B.insert(idx);
         size_B++;
       }
    } else {
       auto it = C.begin();
       advance(it, idx-size_A);
       auto pos = opts.find(*it);
       opts.erase(pos);
       idx = *it;
       C.erase(idx);
       size_C--;
       if(opts.find(idx)!=opts.end()) {
         D.insert(idx);
         size_D++;
       }
    }
    Sk.erase(idx);
    Sk.insert(0);
    path.push_back({idx, 0});
#else
    idx = exchange_B2A_back(Sk, size_A, A, size_B, B, size_C, C, 
                     size_D, D, opts,  path);
#endif

#ifdef __ENABLE_CHECK_
    cout<<"The "<<2*k+1<<" exchange: "<<idx<<" -> 0\n";
    print_set("    S0: ",S0, "\n");
    print_set("    Sk: ",Sk, "\n");
    print_set("    A : ", A, "\n");
    print_set("    B : ", B, "\n");
    print_set("    C : ", C, "\n");
    print_set("    D : ", D, "\n");
    print_set("    op: ", opts, "\n");
#endif
  }

  for(int k=ksteps-1; k>-1; k--) {
    forpath.push_back({path[k].second, path[k].first});
  }

  return 0;
}
#endif

#if 0
/*
 * Direct sampler based on the Kondo DP table. But it does not work at this moment.
 */
void KondoPathSampler::sample_path_DP(const int k, const set<int>& S0, const set<int>& S1, 
         const int k, vector<pair<int, int>>& path, set<int>& final_set) {
  const double inv_M_NM = 1.0/(M*(N-M));
  // 初始化四个集合
  set<int> A, B, C, D;
        
  // A = S0 ∩ S1
  set_intersection(S0.begin(), S0.end(), S1.begin(), S1.end(), 
                   inserter(A, A.begin()));
  // B = S0 - S1
  set_difference(S0.begin(), S0.end(), S1.begin(), S1.end(), 
                 inserter(B, B.begin()));
  // C = S1 - S0
  set_difference(S1.begin(), S1.end(), S0.begin(), S0.end(), 
                 inserter(C, C.begin()));
  // D = [N] - (S0 ∪ S1)
  set<int> all_elements;
  for (int i = 1; i < N; i++) all_elements.insert(i);
  set<int> S0_union_S1;
  set_union(S0.begin(), S0.end(), S1.begin(), S1.end(), 
            inserter(S0_union_S1, S0_union_S1.begin()));
  set_difference(all_elements.begin(), all_elements.end(), 
                 S0_union_S1.begin(), S0_union_S1.end(), 
                 inserter(D, D.begin()));
        
  if(A.find(0)!=A.end()) A.erase(0);
  if(B.find(0)!=B.end()) B.erase(0);
  if(C.find(0)!=C.end()) C.erase(0);

  // 当前状态
  // 计算初始距离
  set<int> sym_diff;
  set_symmetric_difference(S0.begin(), S0.end(), S1.begin(), S1.end(), 
                                 inserter(sym_diff, sym_diff.begin()));
  int d0 = sym_diff.size() / 2;
  int current_d = d0;
  set<int> current_set = S0;
  bool is_A = true;   // 起始为A类

  // 从最后一步开始反向采样
  for (int step = k; step >= 1; step--) {
    if (is_A) {
      // 从A类移动到B类：移除0，添加一个新元素
      double weight_same = 0.0, weight_increase = 0.0;

      // 保持距离的权重 (d -> d)
      if (!B.empty()) {
        weight_same = B.size() * Qtd[step-1][current_d];
      }
      // 远离目标的权重 (d -> d+1)
      if (!D.empty() && current_d < M) {
        weight_increase = D.size() * Qtd[step-1][current_d+1];
      }

      double total_weight = weight_same + weight_increase;
      if (total_weight <= 0) {
        break;  // 没有可行移动
      }

      // 随机选择移动类型
      double r = drand48() * total_weight;
      int e_out = 0;  // 总是移除0
      int e_in;

      if (r < weight_same) {
        // 保持距离：从B中随机选择元素添加
        e_in = get_random_element(B);
                      
        // 更新集合
        B.erase(e_in);
        A.insert(e_in);
      } else {
        // 增加距离：从D中随机选择元素添加
         e_in = get_random_element(D);
                      
        // 更新集合
        D.erase(e_in);
        C.insert(e_in);
        current_d++;  // 距离增加
      }

      // 记录交换
      path.push_back({e_out, e_in});

      // 更新当前集合
      current_set.erase(e_out);
      current_set.insert(e_in);
      is_A = false;
    } else {
      // 从B类移动到A类：添加0，移除一个元素
      double weight_same = 0.0, weight_decrease = 0.0;

      // 保持距离的权重 (d -> d)
      if (!A.empty()) {
         weight_same = A.size() * Ptd[step-1][current_d];
      }

      // 靠近目标的权重 (d -> d-1)
      if (!C.empty() && current_d > 0) {
          weight_decrease = C.size() * Ptd[step-1][current_d-1];
      }

      double total_weight = weight_same + weight_decrease;
      if (total_weight <= 0) {
        break;  // 没有可行移动
      }
                
      // 随机选择移动类型
      double r = drand48() * total_weight;
      int e_in = 0;  // 总是添加0
      int e_out;
                
      if (r < weight_same) {
        // 保持距离：从A中随机选择元素移除
        e_out = get_random_element(A);
        // 更新集合
        A.erase(e_out);
        B.insert(e_out);
      } else {
        // 减少距离：从C中随机选择元素移除
        e_out = get_random_element(C);
        // 更新集合
        C.erase(e_out);
        D.insert(e_out);
        current_d--;  // 距离减少
      }
                
      // 记录交换
      path.push_back({e_out, e_in});
                
      // 更新当前集合
      current_set.erase(e_out);
      current_set.insert(e_in);
      is_A = true;
    }
  }
        
  final_set = current_set;
  return;
}
#endif

#if 0
int main() {
    // 参数设置
    int N = 5;   // 元素总数 {0,1,2,3,4}
    int M = 2;   // 每个集合大小
    int k = 4;   // 总步数
    
    // 定义起始节点 S0 (必须包含0)
    set<int> S0 = {0, 1};
    
    // 定义目标节点 S1
    set<int> S1 = {0, 3};
    
    // 计算目标距离
    set<int> sym_diff;
    set_symmetric_difference(S0.begin(), S0.end(), S1.begin(), S1.end(), 
                             inserter(sym_diff, sym_diff.begin()));
    int d_target = sym_diff.size() / 2;
    
    // 目标节点类型 (由步数奇偶决定)
    bool target_is_A = (k % 2 == 0);  // 偶数步终点为A类
    
    cout << "Kondo图路径采样\n";
    cout << "N=" << N << ", M=" << M << ", 步数=" << k << "\n";
    print_set("起始: ", S0, "\n");
    print_set("目标: ", S1, "\n");
    cout << "目标距离: " << d_target << ", 目标类型: " 
         << (target_is_A ? "A类" : "B类") << "\n\n";
    
    // 创建采样器
    KondoPathSampler sampler(N, M, k, d_target, target_is_A);
    
    // 采样路径
    vector<pair<int, int>> path;
    set<int> final_set;
    sampler.sample_path(S0, S1, path, final_set);
    
    // 打印结果
    cout << "===== 采样路径 =====\n";
    print_path(path, S0);
    
    // 检查是否到达目标
    print_set("最终状态: ", final_set);
    if (final_set == S1) {
        cout << "  (成功到达目标)" << endl;
    } else {
        cout << "  (未到达目标)" << endl;
    }
    
    return 0;
}
#endif

#if 0
int main() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  srand48(tv.tv_sec-800000+9317);

  // 参数设置
  int N = 200;
  int M = 100;
  int k = 36;
  int d = 9;
    
  // 创建采样器
  KondoPathSampler sampler(N, M, 40);

  // 生成初始状态和目标状态
  set<int> S0; 
  for(int j=0; j<M-d; j++) S0.insert(j);
  set<int> S1=S0; 
  for(int j=M-d; j<M; j++) {
    S0.insert(j);
    S1.insert(j+d);
  }

  // 采样路径
  clock_t start = clock();
  vector<pair<int, int>> path;
  int count=0, bad=0, status;
  while(count<1) {
    status = sampler.sample_path_A2A_fb(k, S0, S1, path);
    count += 1+status;
    bad   -= status;
  }

  double sample_time = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
    
  // 验证结果
  print_set("S0 : ", S0, "\n");
  print_set("S1 : ", S1, "\n");
  cout << "10000 A->A 采样路径，失败次数 ："<<bad<<". 耗时: " << fixed << setprecision(6) << sample_time << "秒\n";
  cout << "路径长度: " << path.size() << endl;
  //print_path(path, S0); 
  // 打印交换操作
  count = static_cast<int>(path.size());
  cout << "\n前"<<count<<"个交换操作:\n";
  for (int i = 0; i < count; i++) {
    cout << "步骤 " << (i+1) << ": 移出 " << path[i].first 
         << ", 移入 " << path[i].second << "\n";
  }
  verify_path(S0, S1, path); 

 
#if 0
  S1.erase(0);
  S1.insert(M+d+1);
  print_set("S0 : ", S0, "\n");
  print_set("S1 : ", S1, "\n");
  start = clock();
  status = sampler.sample_path_A2B(k+1, S0, S1, path);
  sample_time = static_cast<double>(clock() - start) / CLOCKS_PER_SEC;
  cout << "A->B 采样路径，耗时: " << fixed << setprecision(6) << sample_time << "秒\n";
  cout << "路径长度: " << path.size() << endl;
  print_path(path, S0); 
  verify_path(S0, S1, path); 
#endif
  return 0;
}
#endif
