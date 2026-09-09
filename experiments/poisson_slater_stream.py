"""Low-memory continuous-time Poisson Slater prototype and small-space QM oracle.

No reference wavefunction, phase filtering, or bath compression is used.
The prefix baseline propagates each pair once through the full interval.
The selected scaled-time estimator rebuilds paths at each observation time
using common randomized integration points; it trades recomputation for lower noise.
The centered observable n_j-n_j(0) has exactly zero contribution if either
path has no jumps. Conditioning each complete Poisson path on >=1 jump is
therefore unbiased after multiplying by the two conditioning probabilities.
Fixed jump-count strata and an open-ended >=4 stratum retain the entire law.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
import argparse, json, math, time, sys, hashlib, platform
import scipy
from pathlib import Path
from itertools import combinations
import numpy as np
import psutil
from scipy.stats import poisson, qmc
from scipy.sparse import csr_matrix, diags, kron, eye
from scipy.sparse.linalg import expm_multiply
from scipy.special import gammaln

EV=27.211386245988

def model(norb=10,nel=5,bath='original'):
    if not 0<nel<norb or norb<3: raise ValueError('require 0 < electrons < orbitals and orbitals >= 3')
    eta=3e-7 if bath=='original' else (.2/EV)**2
    wc=(2.7 if bath=='original' else 4.)/EV
    ef=-4.5/EV
    if bath=='original':
        eb=-10/EV; a=ef-wc/2-eb; b=ef+wc/2-eb
        ratio=(b/a)**(2/3); off=(norb-1-ratio)/(ratio-1)
        en=eb+a*((np.arange(1,norb)+off)/(1+off))**1.5
    else:
        from scipy.optimize import brentq
        def cdf(x): return .5+(x*np.sqrt(max(0,1-x*x))+np.arcsin(x))/np.pi
        en=ef+(wc/2)*np.array([brentq(lambda x:cdf(x)-(j+.5)/(norb-1),-1,1) for j in range(norb-1)])
    frequency=3.6749323758566211e-3; mass=14583.1067146087
    return dict(norb=norb,nel=nel,bath=bath,eta=eta,wc_ev=wc*EV,
        energies=np.r_[-.6319,en].tolist(),coupling=np.r_[0,np.full(norb-1,np.sqrt(eta/(norb-1)))].tolist(),
        frequency=frequency,mass=mass,displacement=2*np.sqrt(mass*frequency/2),
        occupied=list(range(nel)),source_model='v1.14-signed-csr original bath; semicircle is separately labelled')

class Ensemble:
    def __init__(self,m,batch):
        self.m=m;self.n=batch
        self.u=np.broadcast_to(np.eye(m['norb'],dtype=complex)[:,m['occupied']],(batch,m['norb'],m['nel'])).copy()
        self.w=np.ones(batch,complex);self.imp=np.full(batch,0 in m['occupied'],bool)
        self.alpha=np.full(batch,m['displacement'],complex);self.t=np.zeros(batch);self.count=np.zeros(batch,int)
    def advance(self,idx,endpoint):
        if len(idx)==0:return
        endpoint=np.broadcast_to(endpoint,(len(idx),));dt=endpoint-self.t[idx]
        if np.min(dt)<-1e-10:raise AssertionError('negative propagation interval')
        center=self.m['displacement']*self.imp[idx]
        old=self.alpha[idx];new=(old-center)*np.exp(-1j*self.m['frequency']*dt)+center
        self.w[idx]*=np.exp(-.5j*self.m['frequency']*dt+1j*center*(old.imag-new.imag))
        self.alpha[idx]=new
        u=self.u[idx]*np.exp(-1j*dt[:,None]*np.asarray(self.m['energies'])[None,:])[:,:,None]
        has=self.imp[idx]
        self.w[idx[has]]*=u[has,0,0];u[has,0,0]=1
        self.u[idx]=u;self.t[idx]=endpoint
    def jump(self,idx,rate):
        if len(idx)==0:return
        c=np.asarray(self.m['coupling']);nel=self.m['nel'];original_imp=self.imp[idx].copy()
        for occupied in [True,False]:
            ids=idx[original_imp==occupied]
            if len(ids)==0:continue
            u=self.u[ids];out=np.zeros_like(u)
            if occupied:
                block=np.concatenate((np.broadcast_to(c[None,1:,None],(len(ids),len(c)-1,1)),u[:,1:,1:]),axis=2)
                q,r=np.linalg.qr(block,mode='reduced');out[:,1:]=q
                factor=np.prod(np.diagonal(r,axis1=-2,axis2=-1),axis=1)
            else:
                a=np.einsum('r,brj->bj',c,u);pivot=np.argmax(abs(a),axis=1)
                ap=a[np.arange(len(ids)),pivot];safe=np.where(abs(ap)>0,ap,1.)
                keep=np.arange(nel-1)[None,:]+(np.arange(nel-1)[None,:]>=pivot[:,None])
                block=np.take_along_axis(u[:,1:,:],keep[:,None,:],axis=2)
                up=np.take_along_axis(u[:,1:,:],pivot[:,None,None],axis=2)
                ak=np.take_along_axis(a,keep,axis=1)
                block=block-up*(ak/safe[:,None])[:,None,:]
                q,r=np.linalg.qr(block,mode='reduced');out[:,0,0]=1;out[:,1:,1:]=q
                factor=(-1.)**pivot*ap*np.prod(np.diagonal(r,axis1=-2,axis2=-1),axis=1)
            self.u[ids]=out;self.w[ids]*=factor*(-1j/rate)
            self.imp[ids]=not occupied
        self.count[idx]+=1
    def through(self,times,endpoint):
        while True:
            nextt=times[np.arange(self.n),np.minimum(self.count,times.shape[1]-1)]
            idx=np.flatnonzero(nextt<=endpoint)
            if len(idx)==0:break
            self.advance(idx,nextt[idx]);self.jump(idx,self.rate)
        self.advance(np.arange(self.n),endpoint)

def centered_pairs(bra,ket,t,rate):
    result=np.zeros((bra.n,bra.m['norb']))
    ids=np.flatnonzero((bra.count>0)&(ket.count>0)&(bra.imp==ket.imp))
    if len(ids)==0:return result
    a=ket.u[ids];b=bra.u[ids]
    mat=np.einsum('bri,brj->bij',b.conj(),a)
    left,s,right=np.linalg.svd(mat)
    phase=np.linalg.det(left)*np.linalg.det(right)
    before=np.cumprod(np.c_[np.ones(len(ids)),s[:,:-1]],axis=1)
    after=np.cumprod(np.c_[np.ones(len(ids)),s[:,1:][:,::-1]],axis=1)[:,::-1]
    excluded=before*after
    adj=phase[:,None,None]*((right.conj().swapaxes(-1,-2)*excluded[:,None,:])@left.conj().swapaxes(-1,-2))
    overlap=phase*np.prod(s,axis=1)
    occ=np.einsum('bri,bij,brj->br',a,adj,b.conj())
    init=np.zeros(bra.m['norb']);init[bra.m['occupied']]=1
    occ-=overlap[:,None]*init[None,:]
    nuclear=np.exp(-.5*abs(bra.alpha[ids]-ket.alpha[ids])**2+1j*(bra.alpha[ids].conj()*ket.alpha[ids]).imag)
    weight=bra.w[ids].conj()*ket.w[ids]*nuclear*np.exp(2*rate*t)
    result[ids]=(occ*weight[:,None]).real
    return result

def category_probability(label,mu):
    if label=='all':return 1.
    if label=='positive':return -np.expm1(-mu)
    if label=='tail':return poisson.sf(3,mu)
    return poisson.pmf(int(label),mu)

def counts_and_times(label,mu,horizon,u,ordering='sort'):
    if label in ['1','2','3']:counts=np.full(len(u),int(label),int)
    else:
        low=0 if label=='all' else (1 if label=='positive' else 4)
        if low==0:
            counts=poisson.ppf(u[:,0],mu).astype(int)
            counts=np.maximum(counts,0)
        else:
            # Recurrence in the normalized tail avoids 1-survival cancellation
            # in scipy.isf at very small mu (especially >=4 at early times).
            mass=float(np.exp(poisson.logpmf(low,mu)-poisson.logsf(low-1,mu)))
            cdf=[mass];k=low
            while cdf[-1]<float(u[:,0].max(initial=0.)):
                k+=1;mass*=mu/k;updated=cdf[-1]+mass
                if updated==cdf[-1] or k>u.shape[1]-1:
                    raise ArithmeticError('conditional Poisson tail resolution/buffer exhausted')
                cdf.append(updated)
            counts=low+np.searchsorted(cdf,u[:,0],side='left')
    # Normally only a handful of jumps. No silently discarded tail.
    maximum=int(counts.max(initial=0))
    if maximum>u.shape[1]-1:raise ValueError('random dimension insufficient for a sampled Poisson tail; increase --max-events')
    if ordering=='stick':
        remaining=counts[:,None]-np.arange(u.shape[1]-1)[None,:]
        factors=(1-u[:,1:])**(1/np.maximum(remaining,1))
        times=np.where(remaining>0,horizon*(1-np.cumprod(np.where(remaining>0,factors,1.),axis=1)),np.inf)
    else:
        times=np.sort(np.where(np.arange(u.shape[1]-1)[None,:]<counts[:,None],u[:,1:]*horizon,np.inf),axis=1)
    return np.c_[times,np.full(len(u),np.inf)]

def transform_schedule(schedule,horizon):
    finite=np.isfinite(schedule);u=np.where(finite,schedule/horizon,.5)
    jac=np.prod(np.where(finite,6*u*(1-u),1.),axis=1)
    return np.where(finite,horizon*u*u*(3-2*u),np.inf),jac

def run_poisson(m,times,n,seed,mode,batch=128,rate_scale=1.,max_events=24,time_sampling='prefix',ordering='sort'):
    horizon=float(times[-1]);rate=np.linalg.norm(m['coupling'])*rate_scale;mu=rate*horizon
    if mode.startswith('stratified-'):
        cats=['1','2','3','tail'];sectors=[(a,b,1 if i==j else 2) for i,a in enumerate(cats) for j,b in enumerate(cats) if j>=i]
    else:sectors=[(('all' if mode=='plain-mc' else 'positive'),('all' if mode=='plain-mc' else 'positive'),1)]
    delta=np.zeros((len(times),m['norb']));sector_info=[]
    for number,(ca,cb,multiplicity) in enumerate(sectors):
        if time_sampling=='scaled' and ca in ['1','2','3'] and cb in ['1','2','3'] and int(ca)%2!=int(cb)%2:
            sector_info.append(dict(bra=ca,ket=cb,pairs=0,seconds=0.,analytically_zero='opposite impurity parity'))
            continue
        prob=category_probability(ca,mu)*category_probability(cb,mu)*multiplicity
        dim=2*(max_events+1)
        rng=qmc.Sobol(dim,scramble=True,seed=seed+7919*number) if 'rqmc' in mode else np.random.default_rng(seed+7919*number)
        if 'lattice' in mode:
            if time_sampling!='scaled':raise ValueError('periodized lattice requires scaled time strata')
            generator=int(.61803398875*n)//2*2+1
            while math.gcd(generator,n)!=1:generator+=2
            lattice_vector=np.array([pow(generator,i,n) for i in range(dim)],dtype=float)
            lattice_shift=rng.random(dim)
        accum=np.zeros_like(delta);start=time.perf_counter()
        for offset in range(0,n,batch):
            size=min(batch,n-offset)
            if 'lattice' in mode:u=(np.arange(offset,offset+size)[:,None]*lattice_vector[None,:]/n+lattice_shift)%1
            else:u=rng.random(size) if 'rqmc' in mode else rng.random((size,dim))
            # Put the most influential first jump times in Sobol dimensions 0,1.
            # Count uniforms go last; fixed-count strata do not use them.
            u=u[:,[dim-2]+list(range(0,dim-2,2))+[dim-1]+list(range(1,dim-2,2))]
            if time_sampling=='scaled':
                for it,t in enumerate(times):
                    if t==0:continue
                    if ca in ['1','2','3'] and cb in ['1','2','3'] and int(ca)%2!=int(cb)%2:continue
                    ua=u[:,:max_events+1].copy();ub=u[:,max_events+1:].copy()
                    if 'lattice' in mode:
                        ja=6*ua[:,1:]*(1-ua[:,1:]);jb=6*ub[:,1:]*(1-ub[:,1:])
                        ua[:,1:]=ua[:,1:]**2*(3-2*ua[:,1:]);ub[:,1:]=ub[:,1:]**2*(3-2*ub[:,1:])
                    ta=counts_and_times(ca,rate*t,t,ua,ordering);tb=counts_and_times(cb,rate*t,t,ub,ordering)
                    bra=Ensemble(m,size);ket=Ensemble(m,size);bra.rate=ket.rate=rate
                    if 'lattice' in mode:
                        ka=np.isfinite(ta).sum(axis=1);kb=np.isfinite(tb).sum(axis=1)
                        bra.w*=np.prod(np.where(np.arange(max_events)[None,:]<ka[:,None],ja,1.),axis=1)
                        ket.w*=np.prod(np.where(np.arange(max_events)[None,:]<kb[:,None],jb,1.),axis=1)
                    bra.through(ta,t);ket.through(tb,t)
                    prob_t=category_probability(ca,rate*t)*category_probability(cb,rate*t)*multiplicity
                    accum[it]+=prob_t*centered_pairs(bra,ket,t,rate).sum(axis=0)
                continue
            ta=counts_and_times(ca,mu,horizon,u[:,:max_events+1]);tb=counts_and_times(cb,mu,horizon,u[:,max_events+1:])
            bra=Ensemble(m,size);ket=Ensemble(m,size);bra.rate=ket.rate=rate
            for it,t in enumerate(times):
                bra.through(ta,t);ket.through(tb,t)
                accum[it]+=centered_pairs(bra,ket,t,rate).sum(axis=0)
        delta+=(1 if time_sampling=='scaled' else prob)*accum/n
        sector_info.append(dict(bra=ca,ket=cb,probability_with_symmetry=float(prob),pairs=n,seconds=time.perf_counter()-start))
        print(json.dumps(sector_info[-1]),flush=True)
    init=np.zeros(m['norb']);init[m['occupied']]=1
    return init[None,:]+delta,dict(rate=rate,expected_jumps_per_leg=mu,sectors=sector_info,
        sampled_pairs=sum(s['pairs'] for s in sector_info),pairs_per_sampled_stratum=n,observation_count=len(times),batch=batch,time_sampling=time_sampling,ordering=ordering,time_transform=('cubic with exact Jacobian' if 'lattice' in mode else 'uniform'),normalization='known exact unitary norm 1; centered observable; no sample-ratio normalization',
        truncation='none: >=4 jump stratum sampled from full conditional Poisson law; numerical event buffer raises on overflow')

def qm_oracle(m,times,nfock):
    if math.comb(m['norb'],m['nel'])*nfock>2_000_000:
        raise ValueError('Full QM oracle is restricted to small validation spaces; the Poisson solver has no full-basis allocation.')
    basis=list(combinations(range(m['norb']),m['nel']));index={s:i for i,s in enumerate(basis)}
    occ=np.zeros((len(basis),m['norb']));rows=[];cols=[];values=[]
    for source,state in enumerate(basis):
        occ[source,list(state)]=1
        for j in range(1,m['norb']):
            if (0 in state)==(j in state):continue
            annihilate,create=(0,j) if 0 in state else (j,0)
            target=list(state);p=target.index(annihilate);target.pop(p)
            q=sum(x<create for x in target);target.insert(q,create)
            rows.append(index[tuple(target)]);cols.append(source);values.append((-1.)**(p+q)*m['coupling'][j])
    electronic=csr_matrix((values,(rows,cols)),shape=(len(basis),len(basis)))
    en=occ@np.asarray(m['energies']);w=m['frequency'];d=m['displacement'];levels=np.arange(nfock)
    diag=en[:,None]+w*(levels[None,:]+.5)+occ[:,0,None]*w*d*d
    h=kron(electronic,eye(nfock),format='csr')+diags(diag.ravel(),format='csr')
    x=diags([np.sqrt(np.arange(1,nfock)),np.sqrt(np.arange(1,nfock))],[-1,1],shape=(nfock,nfock),format='csr')
    h+=kron(diags(-w*d*occ[:,0]),x,format='csr')
    coherent=np.exp(-d*d/2+levels*np.log(d)-.5*gammaln(levels+1))
    psi=np.zeros((len(basis),nfock),complex);psi[index[tuple(m['occupied'])]]=coherent
    psi=psi.ravel();norm0=float(np.vdot(psi,psi).real)
    shift=float(h.diagonal().mean());generator=-1j*(h-shift*eye(len(psi),format='csr'))
    answer=[];norms=[];last=0.;tail=[]
    for t in times:
        if t>last:psi=expm_multiply(generator*(t-last),psi,traceA=0.)
        density=abs(psi.reshape(len(basis),nfock))**2
        answer.append(density.sum(axis=1)@occ);norms.append(float(density.sum()));tail.append(float(density[:,-8:].sum()));last=t
    return np.asarray(answer),dict(nfock=nfock,full_electronic_states=len(basis),wavefunction_complex_entries=len(psi),
        wavefunction_bytes=psi.nbytes,initial_norm=norm0,max_norm_error=float(np.max(abs(np.array(norms)-1))),
        max_top8_fock_probability=max(tail),propagator='sparse full many-body Fock Hamiltonian expm_multiply; continuous time, no splitting')

def main():
    p=argparse.ArgumentParser();p.add_argument('--method',choices=['qm','plain-mc','conditional-mc','conditional-rqmc','stratified-rqmc','stratified-lattice'],required=True)
    p.add_argument('--bath',choices=['original','semicircle'],default='original');p.add_argument('--norb',type=int,default=10);p.add_argument('--nel',type=int,default=5)
    p.add_argument('--tmax',type=float,default=100);p.add_argument('--output-step',type=float,default=1.)
    p.add_argument('--pairs',type=int,default=1024);p.add_argument('--seed',type=int,default=20260909);p.add_argument('--batch',type=int,default=128)
    p.add_argument('--ordering',choices=['sort','stick'],default='sort')
    p.add_argument('--time-sampling',choices=['prefix','scaled'],default='prefix')
    p.add_argument('--rate-scale',type=float,default=1.);p.add_argument('--max-events',type=int,default=24);p.add_argument('--nfock',type=int,default=192)
    p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    if args.tmax<=0 or args.output_step<=0 or args.pairs<=0 or args.batch<=0 or args.rate_scale<=0 or args.max_events<4 or args.nfock<2:raise ValueError('positive numerical parameters required')
    m=model(args.norb,args.nel,args.bath);times=np.arange(0,args.tmax+args.output_step*.1,args.output_step)
    if abs(times[-1]-args.tmax)>1e-10:raise ValueError('tmax must be a multiple of output-step')
    args.out.mkdir(parents=True,exist_ok=True)
    source=Path(__file__).read_bytes()
    (args.out/'source_snapshot.py').write_bytes(source)
    command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]]
    (args.out/'run_command.json').write_text(json.dumps(command,indent=2),encoding='utf-8')
    started=time.perf_counter();cpu=time.process_time()
    if args.method=='qm':occ,info=qm_oracle(m,times,args.nfock)
    else:occ,info=run_poisson(m,times,args.pairs,args.seed,args.method,args.batch,args.rate_scale,args.max_events,args.time_sampling,args.ordering)
    info.update(wall_seconds=time.perf_counter()-started,cpu_seconds=time.process_time()-cpu,
        peak_rss_bytes=getattr(psutil.Process().memory_info(),'peak_wset',psutil.Process().memory_info().rss),
        current_rss_bytes=psutil.Process().memory_info().rss,model=m,method=args.method,seed=args.seed,
        max_particle_number_error=float(np.max(abs(occ.sum(axis=1)-m['nel']))),
        source_sha256=hashlib.sha256(source).hexdigest(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),processor=platform.processor(),min_occupation=float(occ.min()),max_occupation=float(occ.max()))
    if not np.isfinite(occ).all():raise FloatingPointError('non-finite occupation estimate')
    np.savez_compressed(args.out/'occupations.npz',time=times,occupations=occ)
    np.savetxt(args.out/'occupations.dat',np.c_[times,occ],header='time_au '+ ' '.join('orbital_'+str(i) for i in range(m['norb'])))
    (args.out/'metrics.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    print(json.dumps({k:info[k] for k in ['method','wall_seconds','cpu_seconds','peak_rss_bytes','max_particle_number_error']}),flush=True)
if __name__=='__main__':main()
