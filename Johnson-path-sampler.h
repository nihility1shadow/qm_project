#ifndef _YYY_JOHNSON_PATH_SAMPLER_
#define _YYY_JOHNSON_PATH_SAMPLER_

#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <iomanip>
#include <limits>

#include <yyy_inlines.h>
using namespace std;

inline int MIN(int x, int y) {return x>y? y :x; }
extern double *log_table;

/*
 * NOTICE：The N elements in the full set start from 1 stead of 0
 */
// The class for Johnson Graph Path sampling
class JohnsonPathSampler {
public:
    double **Ptd;
    JohnsonPathSampler(int N0, int M0, int kmax) 
        : Ptd(NULL), N(N0), M(M0), Kmax(kmax), Dmax(MIN(M0, N0 - M0)) {
        if(!log_table) {
          log_table = array1d<double>(2*N0+1);
          log_table[0] = -1.e-100;
          for(int j=1; j<=2*N0; j++) log_table[j] = log(j);
        }
        compute_dp_table() ;
    }
    ~JohnsonPathSampler() { free2d(Ptd); }
    int sample_path(const int k, const set<int>& S0, const set<int>& S1,
        vector<pair<int, int>> &path) const;
    double get_Ptd(const int t, const int d) const { return Ptd[t][d]; }

private:
    int N, M, Kmax, Dmax;
    void compute_dp_table() ;
};
 
// 从集合中随机选择一个元素
int get_random_element(const set<int>& s) ;

// 打印集合
template<class T>
void print_set(const char *heading, const T& s, const char *ending=NULL) {
    if(heading) printf("%s{", heading);
    else        printf("{");
    for (auto it = s.begin(); it != s.end(); ) {
        printf("%d", *it);
        if (++it != s.end()) printf(", ");
    }
    if(ending) printf("}%s", ending);
    else       printf("}");
}

void print_path(const vector<pair<int, int>>& path, const set<int>& S0) ;
bool verify_path(const set<int> &S0, const set<int> &S1, 
            const vector<pair<int, int>>& path);

#endif
