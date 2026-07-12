#include <sys/times.h>
#include <sys/time.h>
#include <cstdlib>

#include <yyy_inlines.h>
#include "./qj.h"
#include "./ahm.h"
#include "./semi.h"
#include "./adt.h"
//#include "./eigen-rep.h"

int diabatic(int argc, char **argv) {
  char tags[12][256] = {"tls", "Tully1", "Tully2", "Tully3", "Tanh", "Rabi", "LZ", "AHM", "AHMZ", "AHMR", "AHAU", "LUZH"};
  char fnm[256];
  FILE *FL; 

  double kappa=0,  // coupling strength
         lambda = 0.0, // the jumping rate of the Poisson process, to be modified later.
         xmin   ,  // Spatial grid lower limits
         xmax   ,  // Spatial grid upper limits
         mass = 1.0,  // mass of the particle
         dt   = 0.01,   // step size of time propagation
         // the parameters for the initial wavepacket
         xinit,    // average position
         kinit,    // average wave vector
         gmm=1.0,  // Gaussian width of the HK propagator
         width;    //width of the initial wavepacket

    /******************************************** 
     * problem = 0 : two-level system;
     * problem = 1 : the Tully model 1
     * problem = 2 : the Tully model 2
     * problem = 3 : the Tully model 3
     * problem = 4 : the tanh potential;
     * problem = 5 : the Rabi model;
     * problem = 6 : the Landau-Zener model;
     * problem = 7 : the Rabi model implemented in the Anderson-Holstein model;
     * problem = 8 : the Anderson-Holstein model;
     * problem = 9 : the Anderson-Holstein model with Zhou's parameters;
     * problem = 10: the Anderson-Holstein model with random parameters;
     * problem = 11: the Example 3 used by LU & Zhou in Math. Comput. 87, 2189, (2017);
     * job = 0, exact quantum mechanical.
     * job = 1, exact Poissonization.
     * job = 2, Poissonization + semiclassical surface hopping method.
     * job = 3, Poissonization + Poor Person's version of the semiclassic hopping.
     * job = 4, Poissonization for Heisenberg equation.
     *******************************************/
  int problem = strtol(argv[1], NULL, 10),
      job     = strtol(argv[2], NULL, 10),
      ntraj   = strtol(argv[3], NULL, 10), //Number of stochastic samples
      nstep ,  // number of time steps
      nwf;     // number of stored wavefunctions
  // The rate factor. jumping rate is rescaled by this variable 
  //   if the offdiagonal matrix elements are constant. 
  // Otherwise, the rate is set to this value.
  double rtfct = 1.0;
  if(argc > 4) rtfct = strtod(argv[4], NULL); 

  _diab_sys_t sys;

  int npt   = 1024; // Grid points
  switch (problem) {
    case 0: // the TLS problem
      nstep  = 2000;
      kappa  = 0.5;
      //BP_const (ntraj, rtfct, nstep, dt);
      BP_tls   (ntraj, 1.0, kappa, rtfct, nstep, dt);
      PJ_tls   (ntraj, 1.0, kappa, rtfct, nstep, dt);
      //PJ_tls_rr(ntraj, 1.0, kappa,        nstep, dt);
      return (0);
    case 1: // Tully model 1
      mass    = 2000.0;
      xinit   = -1;
      kinit   = 20.0;
      width   = 0.5*M_SQRT1_2; // 20*M_SQRT1_2/kinit;
      gmm     = 12.0;
      xmin    = -6;  
      xmax    =  8; 
      npt     = 4096*4;
      nstep   = 200;
      dt      = 0.5;
      nwf     = 200;
      sys.set_mass(mass);
      sys.set_pesup  (&tully1up,   &der_tully1up,   &d2_tully1up,   &d3_tully1up); 
      sys.set_pesdown(&tully1down, &der_tully1down, &d2_tully1down, &d3_tully1down); 
      sys.set_hint   (&tully1int,  &der_tully1int,  &d2_tully1int,  &d3_tully1int ); 
      sys.set_cohint (&coh_tully1int);
      break;
    case 2: // Tully model 2
      mass    = 2000.0;
      xinit   = -5;
      kinit   = 28.5;
      width   = 0.25; //20*M_SQRT1_2/kinit;
      xmin    = -10;  
      xmax    =  10; 
      npt     = 8096;
      nstep   = 400;
      dt      = 0.5;
      nwf     =  80;

      sys.set_mass(mass);
      sys.set_pesup  (&tully2up,   &der_tully2up,   &d2_tully2up,   &d3_tully2up); 
      sys.set_pesdown(&tully2down, &der_tully2down, &d2_tully2down, &d3_tully2down); 
      sys.set_hint   (&tully2int,  &der_tully2int,  &d2_tully2int,  &d3_tully2int ); 
      sys.set_cohint (&coh_tully2int);
      break;
    case 3: // Tully model 3
      mass    = 2000.0;
      xinit   = -4;
      kinit   = 40;
      width   = 0.25; //20*M_SQRT1_2/kinit;
      xmin    = -15;  
      xmax    =  15; 
      npt     = 8096;
      nstep   = 800;
      dt      = 0.5;
      nwf     = 200;

      sys.set_mass(mass);
      sys.set_pesup  (&tully3up,   &der_tully3up,   &d2_tully3up,   &d3_tully3up); 
      sys.set_pesdown(&tully3down, &der_tully3down, &d2_tully3down, &d3_tully3down); 
      sys.set_hint   (&tully3int,  &der_tully3int,  &d2_tully3int,  &d3_tully3int ); 
      sys.set_cohint (&coh_tully3int);
      break;
    case 4: // tanh potentials
      mass    = 1000.0;
      dt      = 0.1;
      nstep   = 250;
      nwf     = 250;
      kappa   = 0.1;  
      xinit   = -3.0;
      kinit   = 800.0;
      width   = 0.5; // 20*M_SQRT1_2/kinit;
      xmin    = -15;  
      xmax    =  15; 
      npt     = 4096*4;

      sys.set_mass(mass);
      sys.set_pesup  (&tanhup,   &der_tanhup,   &d2_tanhup,   &d3_tanhup); 
      sys.set_pesdown(&tanhdown, &der_tanhdown, &d2_tanhdown, &d3_tanhdown); 
      sys.set_hint   (&tanhint,  &der_tanhint,  &d2_tanhint,  &d3_tanhint ); 
      sys.set_cohint (&coh_tanhint);
      sys.set_param(kappa);
      break;
    case 5: // the Rabi problem
      kappa   = 0.05;  
      xmin    = -20;  
      xmax    =  20; 
      nstep   = 2000;
      nwf     = 500;
      npt     = 256; // Grid points
      xinit   = -6.0;
      kinit   = 0.0;
      width   = 1.0;

      sys.set_mass(mass);
      sys.set_pesup  (&rabiup,    &der_rabiup,   &d2_rabiup,   &d3_rabiup); 
      sys.set_pesdown(&rabidown,  &der_rabidown, &d2_rabidown, &d3_rabidown); 
      sys.set_hint   (&Hint_const,&Fzero2p,      &Fzero2p,     &Fzero2p ); 
      sys.set_cohint (&coh_constint);
      sys.set_param(kappa);
      break;
    case 6: // the Landau-Zener model
      kappa   = 0.2;  
      xmin    = -13;  
      xmax    =  17; 
      nstep   = 500;
      nwf     = 100;
      xinit   = -8.0;
      kinit   = 0.0;
      width   = 0.5;

      sys.set_mass(mass);
      sys.set_pesup  (&LZup,    &der_LZup,   &d2_LZup,   &d3_LZup); 
      sys.set_pesdown(&LZdown,  &der_LZdown, &d2_LZdown, &d3_LZdown); 
      sys.set_hint   (&Hint_const,&Fzero2p,      &Fzero2p,     &Fzero2p ); 
      sys.set_cohint (&coh_constint);
      sys.set_param(kappa);
      break;
    case 7: // the Rabi model implemented as a specific case of the Anderson-Holstein problem
      {
        xinit       = -8.0;
        npt         = 512;
        double En   = 0.0,
               delx = 2.0,
               xmin = -20,
               xmax =  20;
        nstep  = 2000;
        AHM ahm;
        ahm.set_mass(1.0);
        ahm.set_freq(1.0);
        ahm.set_delx(delx);
        ahm.set_delE(En);
        kappa = 0.2;  
        En   *=  -1;
        ahm.set_el(1, &En, &kappa);
        ahm.set_grids(npt, xmin, xmax);

        dcomplex alp0 = xinit*M_SQRT1_2;
        ahm.CSpoisson2(ntraj, nstep, dt, alp0);
      }
      return (0);
    case 8: // The Anderson-Holstein problem
      {
        if(argc != 3 && argc != 4) {
          printf("Syntax: %s task=8 ntraj Norb\n", argv[0]);
          abort();
        }
        double wc      = 3.6749323758566216e-01, // 10 eV
               eta     = 0.005,
               delx    = 0.5,
               xmin    =-8.0,
               xmax    =+8.0,
               pinit   = 0.0,
               mass    = 14583.1067146087,       // reduced mass of O2 molecule
               freq    = 3.6749323758566211e-03; // 0.1 eV, O2
        int    Norb    = strtol(argv[4], NULL, 10),
               Nel     = 1;
        npt    = 1024;
        nstep  = 4000;
        dt     = 0.5;
        xinit  = delx;
        AHM ahm;
        ahm.set_mass(mass);
        ahm.set_freq(freq);
        ahm.set_delx(delx);
        ahm.set_grids(npt, xmin, xmax);
        //ahm.discretize(Norb, eta, wc);
        ahm.diseven(Norb, eta, wc);
        ahm.set_Nel(Nel, Norb);
        ahm.set_delE(-9.1873309396415540e-02); //-2.5 eV 
        ahm.set_basis();
        ahm.calc_Eocc();
        ahm.set_exc();
        ahm.set_Nex();

        dcomplex alp0 = sqrt(mass*freq/2)*(xinit+(pinit/(mass*freq))*I);
        if(job==0) ahm.qm(nstep, dt, alp0);
        else ahm.SepCSpoisson(ntraj, nstep, dt, alp0);
      }
      return (0);
    case 9: // The Anderson-Holstein problem with Huang-Zhao parameters
      {
        if(argc != 5) {
          printf("Syntax: %s 9 ntraj Norb 1/epsilon\n", argv[0]);
          abort();
        }
        double epsilon = 1.0/strtod(argv[4], NULL),
               xmin    = -5,
               xmax    = -xmin,
               pinit   = 0.0,
               delx    = -0.2,
               mass    = 1/(epsilon*epsilon),
               freq    = epsilon;
        int    Norb    = 5,
               Nel     = 1; 
        if(argc > 3) Norb = strtol(argv[3], NULL, 10); 
        npt    = 1024;
        nstep  = 4000;
        dt     = 0.004/epsilon;
        AHM ahm;
        ahm.set_mass(mass);
        ahm.set_freq(freq);
        ahm.set_delx(delx);
        ahm.diseven(Norb, epsilon, 1.0);
        ahm.set_Nel(Nel, Norb);
        ahm.set_delE(0.095);
        ahm.set_basis();
        ahm.calc_Eocc();
        ahm.set_grids(npt, xmin, xmax);
        ahm.set_Nex();
        xinit = delx;

        dcomplex alp0 = sqrt(mass*freq/2)*(xinit+(pinit/(mass*freq))*I);
        if(job==0) ahm.qm(nstep, dt, alp0);
        else {
          ahm.SepCSpoisson(ntraj, nstep, dt, alp0);
          //ahm.NHpoisson(ntraj, nstep, dt, alp0);
        }
      }
      return (0);
    case 10: // The Anderson-Holstein problem in atomic units
      {
        if(argc != 6) {
          printf("Syntax: %s task=10 job ntraj Norb Nel\n", argv[0]);
          abort();
        }
        double wc      = 3.6749323758566216e-01, // 10 eV
               eta     = 0.01,
               delx    = 2.0,
               xmin    =-8.0,
               xmax    =+8.0,
               pinit   = 0.0,
               mass    = 14583.1067146087,       // reduced mass of O2 molecule
               freq    = 3.6749323758566211e-03; // 0.1 eV, O2
        int    Norb    = strtol(argv[4], NULL, 10),
               Nel     = strtol(argv[5], NULL, 10);
        npt    = 1024;
        nstep  = 8000;
        dt     = 0.5;
        {
          const char *env_wc_ev = getenv("AHM_WC_EV");
          const char *env_eta   = getenv("AHM_ETA");
          const char *env_nstep = getenv("AHM_NSTEP");
          if(env_wc_ev) wc = atof(env_wc_ev)/27.211386245988;
          if(env_eta) eta = atof(env_eta);
          if(env_nstep) nstep = atoi(env_nstep);
          if(nstep <= 0) {
            printf("AHM_NSTEP must be positive.\n");
            abort();
          }
          printf("#AHAU_PARAMS wc_au=%1.16e wc_eV=%1.16e eta=%1.16e nstep=%d dt=%g\n",
                 wc, wc*27.211386245988, eta, nstep, dt);
        }
        xinit  = delx;
        AHM ahm;
        ahm.set_mass(mass);
        ahm.set_freq(freq);
        ahm.set_delx(delx);
        ahm.set_grids(npt, xmin, xmax);
        //ahm.discretize(Norb, eta, wc);
        ahm.diseven(Norb, eta, wc);
        ahm.set_Nel(Nel, Norb);
        ahm.set_delE(0.375*wc); // aligned with the Fermi level
        ahm.set_basis();
        ahm.calc_Eocc();
        ahm.set_exc();
        ahm.set_Nex();

        dcomplex alp0 = sqrt(mass*freq/2)*(xinit+(pinit/(mass*freq))*I),
                 wgt  = 1.0;
        if(job==0) ahm.qm(nstep, dt, alp0);
        else {
#if 0
          ahm.SepCSpoisson(ntraj, nstep, dt, alp0);
          //ahm.MBpoisson(ntraj, nstep, dt, alp0, wgt, 0);
#else
          set<int> state;
          for(int j=0; j<Nel; j++) state.insert(j);
          ahm.SepMBpoisson(ntraj, nstep, dt, alp0, wgt, state);
#endif
        }
      }
      return (0);
    case 11: // modified tanh potentials
      _LU_eps = 1.0/16;
      mass    = 1.0/(_LU_eps*_LU_eps);
      dt      = 0.005/_LU_eps;
      kappa   = 0.1;  
      xmin    = -7;  
      xmax    = 15; 
      npt     = 2048*4;
      nstep   = 800;
      nwf     = 100;
      xinit   = -1.0;
      kinit   = 2.0/_LU_eps;
      width   = 1.0/sqrt(32.0);
      gmm     = 12.0;

      sys.set_mass(mass);
      sys.set_pesup  (&LUZHup,   &der_LUZHup,   &d2_LUZHup,   &d3_LUZHup); 
      sys.set_pesdown(&LUZHdown, &der_LUZHdown, &d2_LUZHdown, &d3_LUZHdown); 
      sys.set_hint   (&LUZHint,  &der_LUZHint,  &d2_LUZHint,  &d3_LUZHint ); 
      sys.set_cohint (&coh_LUZHint);
      sys.set_param(kappa);
      break;
    default:
      printf("problem should be 0...11.\n");
      exit(1);
  }


  if(problem==5 && job==1) {
    RabiPoisson(ntraj, nstep, dt, npt, xmin, xmax, 2.0, 0.0, kappa, xinit);
    //RabiNHPP(ntraj, nstep, dt, npt, xmin, xmax, 2.0, 0.0, kappa, xinit);
    return 0;
  }

  double dx      = (xmax - xmin)/npt, // spatial grid size
         **sigma = array2d<double>(2, nwf+1); 
       // sigma : the deviation in psi1 (sigma[0]) and psi2 (sigma[1])

  dcomplex *psi1 = array1d<dcomplex>(npt), 
           *kap1 = array1d<dcomplex>(npt), 
           *psi2 = array1d<dcomplex>(npt), 
           *kap2 = array1d<dcomplex>(npt),
           *qwvf = array1d<dcomplex>(npt);
  double   *v1   = array1d<double>(npt),
           *v2   = array1d<double>(npt),
           *vint = array1d<double>(npt);

  double   *Eg   = array1d<double>(npt),
           *Ee   = array1d<double>(npt),
           *Vna  = array1d<double>(npt),
           *cst  = array1d<double>(npt),
           *snt  = array1d<double>(npt),
           *Kge = array1d<double>(npt);

  if(myid==master) {
    sprintf(fnm, "%s-pes.dat", tags[problem]);
    FL = fopen(fnm, "w");
    fprintf(FL, "#%11s %16s %16s %16s %16s %16s %16s\n",
        "t", "cst", "snt", "Eg", "Ee", "Kge", "Vna");
  }

  //initializing the wavepacket
  dcomplex norx = 0.0;
  for(int i=0; i<npt; i++) {
    double x = xmin + i*dx;
    v1[i]   = sys.calc_pesup(x);
    v2[i]   = sys.calc_pesdown(x);
    vint[i] = sys.calc_Hint(x, kappa);
    lambda += abs(vint[i]);
    qwvf[i] = init_gaussian(x, xinit, kinit, width); 
    norx   += qwvf[i]*conj(qwvf[i]);
  }
  if(myid==master) fclose(FL);

  norx = 1./csqrt(norx*dx);
  for(int i=0; i<npt; i++)  *(qwvf + i) *= norx ; 
  lambda *= dx*rtfct/(xmax-xmin);
  //if(myid==master)printf("kappa = %1.16e, lambda = %1.16e\n", kappa, lambda);

  fftw_plan fftpln_for1, fftpln_bck1, fftpln_for2, fftpln_bck2;
  fftpln_for1 = fftw_plan_dft_1d(npt, 
                  (fftw_complex *)psi1, (fftw_complex *)kap1, 
                  FFTW_FORWARD,  FFTW_MEASURE);
  fftpln_bck1 = fftw_plan_dft_1d(npt, 
                  (fftw_complex *)kap1, (fftw_complex *)psi1, 
                  FFTW_BACKWARD, FFTW_MEASURE);     
  fftpln_for2 = fftw_plan_dft_1d(npt, 
                  (fftw_complex *)psi2, (fftw_complex *)kap2, 
                  FFTW_FORWARD,  FFTW_MEASURE);
  fftpln_bck2 = fftw_plan_dft_1d(npt, 
                  (fftw_complex *)kap2, (fftw_complex *)psi2, 
                  FFTW_BACKWARD, FFTW_MEASURE);     
  memcpy(psi1, qwvf, sizeof(dcomplex)*npt);
  bzero(psi2, sizeof(dcomplex)*npt);

//  memcpy(psi2, qwvf, sizeof(dcomplex)*npt);
  double ***rlt = array3d<double>(2, nstep+1, 4),
         ***sig = array3d<double>(2, nstep+1, 3);
#if 0
        sprintf(fnm, "wf-qm-0.dat");
        FL = fopen(fnm, "w");
        fprintf(FL, "#%12s %25s %25s %25s %25s\n", "x", 
            "Re[psi1]", "Im[psi1]", "Re[psi2]", "Im[psi2]");
        for(int j=0; j<npt; j++) {
          double x = xmin + j*dx;
          fprintf(FL, "%+12.8f %+25.16e %+25.16e %+25.16e %+25.16e\n", x, 
              real(psi1[j]), imag(psi1[j]),
              real(psi2[j]), imag(psi2[j]));
        }
        fclose(FL);
#endif

  struct tms st_cpu, en_cpu, pr_cpu;
  clock_t st_time = times(&st_cpu);

#if 1
  if(myid==master) {
#if 0
    pnpsi(1, xmin, xmax, npt, psi1, kap1, fftpln_for1, fftpln_bck1);
    double p = 0.0, prob=0.0;
    for(int n=0; n<npt; n++) {
      p += real(psi1[n]*conj(qwvf[n]))*dx;
      prob += real(qwvf[n]*conj(qwvf[n]))*dx;
    }
    printf("%e\n", p/prob);
    memcpy(psi1, qwvf, sizeof(dcomplex)*npt);
#endif

    qav(xmin, xmax, npt, psi1, 0.0, 1.0, rlt[0][0]);
    qav(xmin, xmax, npt, psi2, 0.0, 1.0, rlt[0][0]+2);
    correlation(xmin, xmax, npt, psi1, psi2, 0.0, 1.0, sig[0][0]);
  }
  for(int l=0; l<nstep; l++) {
    QJqm(xmin, xmax, npt, l*dt, dt, mass, v1, v2, vint,
      psi1, kap1, fftpln_for1, fftpln_bck1,
      psi2, kap2, fftpln_for2, fftpln_bck2);
    if(myid==master) {
      qav(xmin, xmax, npt, psi1, 0.0, 1.0, rlt[0][l+1]);
      qav(xmin, xmax, npt, psi2, 0.0, 1.0, rlt[0][l+1]+2);
      correlation(xmin, xmax, npt, psi1, psi2, 0.0, 1.0, sig[0][l+1]);
#if 0
      if((l+1)%100==0) {
        sprintf(fnm, "wf-qm-%05d.dat", l+1);
        FL = fopen(fnm, "w");
        fprintf(FL, "#%12s %25s %25s %25s %25s\n", "x", 
            "Re[psi1]", "Im[psi1]", "Re[psi2]", "Im[psi2]");
        for(int j=0; j<npt; j++) {
          double x = xmin + j*dx;
          fprintf(FL, "%+12.8f %+25.16e %+25.16e %+25.16e %+25.16e\n", x, 
              real(psi1[j]), imag(psi1[j]),
              real(psi2[j]), imag(psi2[j]));
        }
        fclose(FL);
      }
#endif
    }
  }
#endif
  clock_t pr_time = times(&pr_cpu);

#if 0
    memcpy(psi2, psi1, sizeof(dcomplex)*npt);
    pnpsi(1, xmin, xmax, npt, psi1, kap1, fftpln_for1, fftpln_bck1);
    double p = 0.0, prob=0.0;
    for(int n=0; n<npt; n++) {
      p    += real(psi1[n]*conj(psi2[n]))*dx;
      prob += real(psi2[n]*conj(psi2[n]))*dx;
    }
    printf("%e\n", p/prob);
#endif

  memcpy(psi1, qwvf, sizeof(dcomplex)*npt);
  bzero(psi2, sizeof(dcomplex)*npt);
  dcomplex **wf1 = array2d<dcomplex>(nwf+1, npt),
           **wf2 = array2d<dcomplex>(nwf+1, npt);
  if(job==0) {
    sprintf(fnm,  "avqm-%s-%d.dat", tags[problem],  ntraj);
  } else if(job==1) {
    sprintf(fnm, "avpp-%s-%d.dat", tags[problem],   ntraj);
    NHpoisson(ntraj, vint, nstep, dt, 
      nwf, wf1, wf2, sigma,
      xmin, xmax, npt, mass, v1, v2, vint,
      psi1, kap1, fftpln_for1, fftpln_bck1,
      psi2, kap2, fftpln_for2, fftpln_bck2);
  } else if(job==2) {
    sprintf(fnm, "avsemi-%s-%d.dat", tags[problem], ntraj);
    sys.hksh(ntraj, nstep, nstep*dt, xinit, kinit, gmm, 
        npt, xmin, xmax, nwf, wf1, wf2) ;
  } else if(job==3) {
    sprintf(fnm, "avpphk-%s-%d.dat", tags[problem], ntraj);
    sys.hksh_ppv(ntraj, nstep, nstep*dt, xinit, kinit, gmm, 
        npt, xmin, xmax, nwf, wf1, wf2) ;
  } else if(job==4){
    sprintf(fnm, "avlvl-%s-%d.dat", tags[problem], ntraj);
    sys.HeisenbergHK(ntraj, nstep, nstep*dt, xinit, kinit, gmm, fnm);
  } else if(job!=0){
    fprintf(stderr, "job should be 0, 1, 2, 3. But now job = %d\n", job);
    abort();
  }

  clock_t en_time = times(&en_cpu);

  if(myid == master) {
    char str[256];
    sprintf(str, "wf-%s-qj-%d.dat", tags[problem], abs(ntraj));
    FL = fopen(str, "w");
    fprintf(FL, "%12s %25s %25s %25s %25s\n", "x", 
        "Re[psi1]", "Im[psi1]", "Re[psi2]", "Im[psi2]");
    for(int j=0; j<npt; j++) {
      double x = xmin + j*dx;
      fprintf(FL, "%+12.8f %+25.16e %+25.16e %+25.16e %+25.16e\n", x, 
          real(wf1[nwf][j]), imag(wf1[nwf][j]),
          real(wf2[nwf][j]), imag(wf2[nwf][j]));
    }
    fclose(FL);

    if(job ==4 ) { 
      FL = fopen(fnm, "a");
    } else {
      FL = fopen(fnm, "w");
    }
  }

  print_timing(st_time, pr_time, en_time,
               st_cpu,  pr_cpu,  en_cpu, FL);

  if(myid == master && job == 4) { fclose(FL);}

  if(myid == master && job != 4) {
    int njunk = nstep/nwf;
    fprintf(FL, "#%11s %25s %25s %25s %25s %25s %25s %25s %25s %25s %25s %25s %25s %25s %25s %25s %25s\n", 
        "t", 
        "QM <psi1|psi1>", "QM <psi1|x|psi1>", 
        "QM <psi2|psi2>", "QM <psi2|x|psi2>",
        "QM <psi|sigma_x|psi>", 
        "QM <psi|sigma_y|psi>", 
        "QM <psi|sigma_z|psi>", 
        "PJ <psi1|psi1>", "PJ <psi1|x|psi1>", 
        "PJ <psi2|psi2>", "PJ <psi2|x|psi2>",
        "PJ <psi|sigma_x|psi>", 
        "PJ <psi|sigma_y|psi>", 
        "PJ <psi|sigma_z|psi>",
        "dev<psi1||psi1>",
        "dev<psi2||psi2>"
    );
    for(int k=0; k<=nwf; k++){
      int j = k*njunk;
      double sd1=0.0, sd2=0.0;
      for(int l=0; l<npt; l++) {
        sd1 += real(wf1[k][l]*conj(wf1[k][l]));
        sd2 += real(wf2[k][l]*conj(wf2[k][l]));
      }
      sd1 = sqrt(abs(sigma[0][k] - sd1*dx)/ntraj);
      sd2 = sqrt(abs(sigma[1][k] - sd2*dx)/ntraj);
      qav(xmin, xmax, npt, wf1[k], 0.0, 1.0, rlt[1][j]);
      qav(xmin, xmax, npt, wf2[k], 0.0, 1.0, rlt[1][j]+2);
      correlation(xmin, xmax, npt, wf1[k], wf2[k], 0.0, 1.0, sig[1][j]);
      fprintf(FL, "%+12.8f %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e %+25.16e\n", 
               j*dt, 
               rlt[0][j][0], rlt[0][j][1], rlt[0][j][2], rlt[0][j][3],
               sig[0][j][0], sig[0][j][1], sig[0][j][2],
               rlt[1][j][0], rlt[1][j][1], rlt[1][j][2], rlt[1][j][3],
               sig[1][j][0], sig[1][j][1], sig[1][j][2],
               sd1, sd2
      );
#if 0
    double pold = 0;
      if(j) {
        double rate = (rlt[0][j][3]-pold)/(dt*pold);
        printf("%+12.8f %+25.16e\n", j*dt, rate);
      }
      pold = rlt[0][j][3];
#endif
    }
    fclose(FL);
  }

  fftw_destroy_plan(fftpln_for1);
  fftw_destroy_plan(fftpln_bck1);
  fftw_destroy_plan(fftpln_for2);
  fftw_destroy_plan(fftpln_bck2);
  free(psi1);
  free(kap1);
  free(psi2);
  free(kap2);
  free2d(wf1);
  free2d(wf2);
  free(v1);
  free(v2);
  free2d(sigma);
  free3d(rlt);
  free3d(sig);

  return 0;
}
