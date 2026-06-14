#ifndef _YYY_AHM_
#define _YYY_AHM_

#include <set>

#include <yyy_inlines.h>
#include <cxxsparse_double.h>

/*The Anderson-Holstein model
 * H = \hat{p}^2/2 m + U0(\hat{x}) + (U1-U0)\hat{d}^\dagger\hat{d}
 *    + \sum_k (En[k] \hat{c}_k^\dagger\hat{c}_k + cpl[k] \hat{c}_k^\dagger \hat{d} + cpl[k] \hat{d}^\dagger\hat{c}_k)
 *    U0 = 1/2 m freq^2 x^2
 *    U1 = 1/2 m freq^2 (x-delx)^2 + delE
 */
class AHM {
  private:
    int    Nel;   // number of electrons in the systems, including the one in the molecule
    int    Norb;  // nubmer of electronic orbitals, being the number in the continuum plus 1 (the excited state)
    int    Nhs;   // the size of the Hilbert space
    int    npt;   // number of grids for data analysis
    int    Nex;   // Nex is the first integer to satisfy occ[Nex][0]==1
    double *En;   // energies of the electronic orbitals in the continuum
    double *cpl;  // the electron-vibration coupling in the continuum
    double *Eocc; // the energy of the state occ
    int   **occ;  // the occupied orbitals with Norb electrons
    int  **virt;  // virt[j] is an array storing the virtual orbitals corresponding to the occupied state occ[j].
                  // The combination of virt[j] occ[j] is {0, ..., Norb}.
    double mass;  // the mass of the vibration
    double freq;  // the frequency of the vibration
    double delx;  // the displacement between the ground and excited state PES.
    double delE;  // the energy difference between bottoms of the ground and excited state PES
    double xmin;  // for data analysis
    double xmax;  // for data analysis
    DSPMatrix exc;  // exc[j][k] is the matrix element of the j<->k transition; 
                  //<k| \sum_l g_l(\hat{c}_l^\dagger \hat{d} + \hat{d}^\dagger\hat{c}_l) |j> = exc[k][j] |k>
  public:
    AHM() {Nel=0; Norb = 0; En=0; cpl=0; occ=0; virt=0; Eocc=0; mass=0; freq=0; delx=0; delE=0;}
    ~AHM() {free1d(En); free1d(cpl); free1d(Eocc); free2d(occ); free2d(virt);}
    void set_el(const int n, double *const e, double *const c, const int mel=1) {
      Nel  = mel;
      Norb = n;
      En  = array1d<double>(Norb);
      cpl = array1d<double>(Norb);
      memcpy(En,  e, sizeof(double)*n);
      memcpy(cpl, c, sizeof(double)*n);
      return;
    }
    void set_basis();
    void set_exc();
    CMatrix calc_U1t(const double dt) const;
    void set_Nel (const int ne, const int nlevels)    { Nel  = ne; Norb = nlevels; return;}
    void set_mass(const double m) { mass = m; return;}
    void set_freq(const double f) { freq = f; return;}
    void set_delx(const double x) { delx = x; return;}
    void set_delE(const double e) { delE = e; En[0] = e; return;}
    void calc_Eocc();
    void set_grids(const int n, const double x0, const double x1) { 
      npt =n; 
      xmin=x0; 
      xmax=x1; 
      return;
    }
    void set_Nex() {
      Nex = 1;
      int den=1;
      for(int n=1; n<Nel; n++) {
        Nex *= Norb-n;
        den *= n;
      }
      Nex /= den;
      if(occ && (occ[Nex-1][0] != 0 || occ[Nex][0] != 1)) {
        printf("something goes wrong or the structure of basis set.\n");
        abort();
      }
      return;
    }

    void U1(const double dt, dcomplex **psi) const ;
    void exwf(const int npt, const double Fa, const double Fdt, 
       dcomplex *const psi, dcomplex *phi) const ;
    void MBexwf(const int npt, const double Fa, const double Fdt, 
       dcomplex *const psi, dcomplex *phi) const ;
    void qm(const int nstep, const double dt, const dcomplex alp0) const ;
    void qm2(const int nstep, const double dt, const dcomplex alp0) const ;
    void CSpoisson  (const int ntraj, const int nstep, const double dt, const dcomplex alp0) const ;
    void Semipoisson(const int ntraj, const int nstep, const double dt, const dcomplex alp0) const ;
/*
 * CSpoisson2 is almost the same as CSpoisson except it 
 * calculated the average with respect to the electronic level #1
 */
    void CSpoisson2(const int ntraj, const int nstep, const double dt, const dcomplex alp0) const ;
    // NHpoisson is working now.
    void NHpoisson (const int ntraj, const int nstep, const double dt, const dcomplex alp0) const ;
    // Sepoisson is not successful
    //void Sepoisson (const int ntraj, const int nstep, const double dt, const dcomplex alp0) const ;
    void SepCSpoisson (const int ntraj, const int nstep, const double dt, const dcomplex alp0) const ;
    void Inhomopoisson (const int ntraj, const int nstep, const double dt, const dcomplex alp0) const ;
    // NHSepCSpoisson works but its performance is worse than that of SepCSpoisson
    void NHSepCSpoisson (const int ntraj, const int nstep, const double dt, const dcomplex alp0) const ;
    void MBpoisson(const int ntraj, const int nstep, const double dt, 
              const dcomplex alp_init, const dcomplex amp, const int state_init) const ;
    void SepMBpoisson(const int ntraj, const int nstep, const double dt, 
              const dcomplex alp_init, const dcomplex amp, const set<int> &state) const ;
    void wfanalysis(FILE *FL, const double t, dcomplex **const wf) const ;
    void csproj(const dcomplex wgt_for, const dcomplex alp_for, 
                const dcomplex wgt_bck, const dcomplex alp_bck, 
                const int excited, const set<int> &state, double *rlt) const ;
    void discretize(const int Nstate, const double eta, const double wc);
    void diseven (const int Nstate, const double eps, const double Emax);
    void diseven2(const int Nstate, const double eps, const double Emax);
    void disrandom(const int Nstate, const double eta, const double wc);
};

int binary_search3(const void *key, void *base, int nel, int width,
        int *found, int (*compar)(const void *, const void *, size_t));
int intcmp (const void *a, const void *b) ;
int _mycmp3(const void *a, const void *b, size_t __cmp_num);
void complement(const int Nall, int *aux, const int Nx, int *const x, int *y) ;

#endif
