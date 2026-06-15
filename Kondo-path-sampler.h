#ifndef _YYY_KONDO_PATH_SAMPLER_
#define _YYY_KONDO_PATH_SAMPLER_

#include <vector>
#include <set>
#include <sys/time.h>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <random>
#include <cmath>

#include <yyy_inlines.h>

#include "./Johnson-path-sampler.h"

using namespace std;

/*
 * 现有一个图，图上每个节点由M个不同自然数构成的集合代表，
 * 其中每个数有N个可能，即0，1，2，..., N-1。
 * 两个节点当且仅当它们的集合的对称差集有2个数，且其中一个为0时，
 * 它们之间有边。这时节点可分两类：A类（包含0）和B类（不包含0）。
 * 我们可将这类图称为Kondo图，因为它与凝聚态物理中的Kondo问题密切相关。
 * Kondo图上单次移动总会从A类中的一个节点到B类中的一个，或者反过来。
 * 如果偶数次移动，A类中的节点总会到A类中的另一个，B类中的节点也不会跑到A类中。
 * 经过分析，从Kondo图上A类节点S0开始移动时，到达另一个节点S1的概率只与
 * S0和S1对称差的元素个数的一半d和移动次数n有关，可写为P(n,d)。
 * 移动次数为奇或偶时情况有些不同，递推公式为：
 * strat from A      
 * P(2k, d)   =   (M-d+1)/M P(2k-1, d-1) + d/M P(2k-1, d)
 * P(2k+1, d) =   (d+1)/(N-M) P(2k,   d+1) + (N-M-d)/(N-M) P(2k,   d),
 * strat from B
 * Q(2k+1  , d) = (M-d+1)/M Q(2k, d-1) + d/M Q(2k, d),
 * Q(2k  , d) =   (d+1)/(N-M) Q(2k-1, d+1) + (N-M-d)/(N-M) Q(2k-1, d)
 * 其中 P(n,d') 前的系数是转移概率，只依赖N, M, d 的函数，跟k无关。
 *
 * 以下，我们将从A类节点出发，k 次交换后获得一个距离为d的节点的概率
 * 记为 Ptd[k][d]，从B类节点出发得到 的概率记为 Qtd[k][d].
 * Ptrans[0] = Ptd, Ptrans[1] = Qtd;
 *
 * 在Kondo图上的游走，
 * 如果是从A类开始 2 k 次游走，问题等价于在 Jonhnson(N-1,M-1)图上的k次游走；
 * 如果是从B类开始 2 k 次游走，问题等价于在 Jonhnson(N-1,M)图上的k次游走。
 *
 * 高效起见，Kondo 图上的路径搜索被划归为Johnson图的问题。
 */
class KondoPathSampler {
  private :
    int N, M, Kmax, Dmax;
    double **Ptd; 
    double **Qtd;
    double ***Ptrans;
    JohnsonPathSampler JPS_N1_M1; //N-1,M-1,Kmax, for A2A
    JohnsonPathSampler JPS_N1_M; //N-1,M,Kmax, for B2B
    void compute_dp_table();
    void compute_dp_table_sum();
  public :
    KondoPathSampler(const int n0, const int m0, const int km):
      N(n0), M(m0), Kmax(km), Dmax(min(m0, n0-m0)),
      Ptd(NULL), Qtd(NULL), Ptrans(NULL),
      JPS_N1_M1(n0-1, m0-1, km), JPS_N1_M(n0-1, m0, km) {
      compute_dp_table();
    }
    ~KondoPathSampler() { free3d(Ptrans); Ptd=NULL; Qtd=NULL; }
    double get_Ptd(const int k, const int d) {return Ptd[k][d];}
    double get_Qtd(const int k, const int d) {return Qtd[k][d];}
    int sample_path_A2A(const int k, const set<int> &S0, const set<int> &S1, 
        vector<pair<int, int> > &path) const;
    int sample_path_A2B(const int k, const set<int> &S0, const set<int> &S1, 
        vector<pair<int, int> > &path) const;
    int sample_path_B2A(const int k, const set<int> &S0, const set<int> &S1, 
        vector<pair<int, int> > &path) const;
    int sample_path_B2B(const int k, const set<int> &S0, const set<int> &S1, 
        vector<pair<int, int> > &path) const;
#if 1
    // this function doest not work yet.
    int sample_path_A2A_uniform(const int k, const set<int> &S0, 
        const set<int> &S1, vector<pair<int, int> > &path) const;
    int sample_path_A2A_backward(const int k, const set<int> &S0, 
        const set<int> &S1, vector<pair<int, int> > &path) const;
    int sample_path_A2A_fb(const int ksteps, const set<int> &S0, 
        const set<int> &S1, vector<pair<int, int> > &path) const;
#endif
};

#endif
