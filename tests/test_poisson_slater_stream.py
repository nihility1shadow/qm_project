"""Independent checks for the streaming Poisson prototype; no fitted QM input."""
import sys,json
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'))
from poisson_slater_stream import Ensemble,model,centered_pairs,counts_and_times,transform_schedule,category_probability
from coherent_poisson_kernel import CoherentSlaterPath,pair_observables
from slater_path_kernel import SlaterPath
from scipy.stats import poisson
from scipy.special import roots_legendre

def validate():
 rng=np.random.default_rng(11601);stateerr=obserr=0.;checks=0
 sizes=[(3,1),(4,2),(10,5),(30,15),(50,25),(50,1),(50,49)]
 for norb,nel in sizes:
  m=model(norb,nel);rate=np.linalg.norm(m['coupling']);batch=4
  ens=[Ensemble(m,batch),Ensemble(m,batch)]
  old=[[CoherentSlaterPath(SlaterPath.basis(norb,m['occupied']),m['displacement']) for _ in range(batch)] for _ in range(2)]
  for step in range(5):
   for side in range(2):
    ens[side].advance(np.arange(batch),float(step+1)*7)
    for i in range(batch):old[side][i]=old[side][i].advance(m['energies'],m['frequency'],m['displacement'],7.)
    ids=np.flatnonzero(rng.random(batch)<.7);ens[side].jump(ids,rate)
    for i in ids:old[side][i]=old[side][i].jump(m['coupling'],-1j/rate)
    for i in range(batch):
     err=max(float(np.max(abs(ens[side].u[i]-old[side][i].electronic.orbitals))),float(abs(ens[side].w[i]-old[side][i].electronic.weight)),float(abs(ens[side].alpha[i]-old[side][i].alpha)))
     stateerr=max(stateerr,err);checks+=1
   obs=centered_pairs(ens[0],ens[1],(step+1)*7,rate)
   init=np.zeros(norb);init[m['occupied']]=1
   for i in range(batch):
    ol,oc=pair_observables(old[0][i],old[1][i]);target=((oc-init*ol)*np.exp(2*rate*(step+1)*7)).real
    obserr=max(obserr,float(np.max(abs(target-obs[i]))))
 assert stateerr<1e-10 and obserr<1e-10,(stateerr,obserr)
 # Independently integrate transformed polynomial moments. Missing/doubled
 # Jacobians or incorrect endpoint mapping would fail these identities.
 x,w=roots_legendre(48);u=(x+1)/2;w=w/2
 times,jac=transform_schedule(np.c_[u,np.full(len(u),np.inf)],1.)
 transform_error=max(abs(np.dot(w,jac*times[:,0]**k)-1/(k+1)) for k in range(9))
 assert transform_error<2e-14
 # Check stable conditional Poisson quantiles against direct PMF sums,
 # including probabilities too small for inverse-survival subtraction.
 tailchecks=0
 for mu in [1e-6,.000548,.05,.7,2.,5.]:
  for label,low in [('positive',1),('tail',4)]:
   grid=(np.arange(4096)+.5)/4096
   un=np.full((4096,65),.25);un[:,0]=grid
   schedules=counts_and_times(label,mu,1.,un);counts=np.isfinite(schedules).sum(axis=1)
   survival=poisson.sf(low-1,mu)
   cdf=0.
   for k in range(low,int(counts.max())+1):
    prev=cdf;cdf+=np.exp(poisson.logpmf(k,mu)-np.log(survival))
    selected=grid[counts==k]
    if len(selected):assert selected.min()>=prev-2e-12 and selected.max()<=cdf+2e-12
    tailchecks+=1
   assert abs(sum(category_probability(c,mu) for c in ['1','2','3','tail'])-poisson.sf(0,mu))<2e-14
 return dict(checks=checks,max_state_difference=stateerr,max_observable_difference=obserr,
             tested_orbitals_electrons=sizes,jacobian_moment_max_error=float(transform_error),conditional_poisson_bins_checked=tailchecks,
             status='PASS; kernel, estimator and sampler checks; physical accuracy assessed separately against full QM')
if __name__=='__main__':
 result=validate();dest=Path(__file__).resolve().parents[1]/'scan_slater_stream_100au_20260909';dest.mkdir(exist_ok=True)
 (dest/'kernel_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result))
