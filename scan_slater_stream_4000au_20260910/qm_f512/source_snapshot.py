"""Checkpointed time-block execution of the validated Poisson Slater estimator.

Only batching and persistence change. The Poisson law, lattice points, Jacobians,
physical Hamiltonian, and centered observable are inherited from the 100-au code.
"""
import os,json,sys,time,math,hashlib,argparse,platform
from pathlib import Path
import poisson_slater_stream as base
np=base.np

class PlannedStop(Exception): pass

class Checkpoint:
    def __init__(self,out,configuration):
        self.out=Path(out);self.out.mkdir(parents=True,exist_ok=True)
        self.configuration=configuration;self.previous={};self.arrays={}
        self.start=time.perf_counter();self.cpu=time.process_time()
        path=self.out/'checkpoint.npz'
        if path.exists():
            with np.load(path,allow_pickle=False) as z:
                self.previous=json.loads(str(z['metadata']))
                if self.previous['configuration']!=configuration:raise ValueError('Checkpoint configuration/source mismatch')
                self.arrays={k:z[k].copy() for k in z.files if k!='metadata'}
        self.prior_wall=self.previous.get('wall_seconds',0.)
        self.prior_cpu=self.previous.get('cpu_seconds',0.)
        self.prior_peak=self.previous.get('peak_rss_bytes',0)
    def resources(self):
        mi=base.psutil.Process().memory_info()
        return dict(wall_seconds=self.prior_wall+time.perf_counter()-self.start,
            cpu_seconds=self.prior_cpu+time.process_time()-self.cpu,
            peak_rss_bytes=max(self.prior_peak,getattr(mi,'peak_wset',mi.rss)),
            current_rss_bytes=mi.rss,timing_scope='solver and checkpoint I/O; interrupted work since last checkpoint is not counted')
    def save(self,progress,**arrays):
        info=dict(configuration=self.configuration,progress=progress,**self.resources())
        tmp=self.out/'checkpoint.tmp.npz'
        np.savez_compressed(tmp,metadata=json.dumps(info),**arrays)
        os.replace(tmp,self.out/'checkpoint.npz')
        q=self.out/'progress.tmp.json';q.write_text(json.dumps(info,indent=2),encoding='utf-8')
        os.replace(q,self.out/'progress.json')
        print(json.dumps(dict(progress=progress,**self.resources())),flush=True)


def poisson_blocked(m,times,n,seed,batch=128,time_block=8,rate_scale=1.,max_events=24,checkpoint=None,stop_after_batches=None):
    cats=['1','2','3','tail']
    sectors=[(a,b,1 if i==j else 2) for i,a in enumerate(cats) for j,b in enumerate(cats) if j>=i]
    rate=float(np.linalg.norm(m['coupling'])*rate_scale);dim=2*(max_events+1)
    generator=int(.61803398875*n)//2*2+1
    while math.gcd(generator,n)!=1:generator+=2
    vector=np.array([pow(generator,i,n) for i in range(dim)],float)
    permutation=[dim-2]+list(range(0,dim-2,2))+[dim-1]+list(range(1,dim-2,2))
    contributions=np.zeros((len(sectors),len(times),m['norb']))
    offsets=np.zeros(len(sectors),int)
    if checkpoint is not None and checkpoint.arrays:
        contributions=checkpoint.arrays['contributions'];offsets=checkpoint.arrays['offsets']
    blocks=0
    for number,(ca,cb,mult) in enumerate(sectors):
        if ca!='tail' and cb!='tail' and int(ca)%2!=int(cb)%2:
            offsets[number]=n;continue
        shift=np.random.default_rng(seed+7919*number).random(dim)
        for offset in range(int(offsets[number]),n,batch):
            size=min(batch,n-offset)
            u=(np.arange(offset,offset+size)[:,None]*vector[None,:]/n+shift)%1
            u=u[:,permutation]
            ua=u[:,:max_events+1].copy();ub=u[:,max_events+1:].copy()
            ja=6*ua[:,1:]*(1-ua[:,1:]);jb=6*ub[:,1:]*(1-ub[:,1:])
            ua[:,1:]=ua[:,1:]**2*(3-2*ua[:,1:]);ub[:,1:]=ub[:,1:]**2*(3-2*ub[:,1:])
            for first in range(1,len(times),time_block):
                ts=times[first:first+time_block]
                sa=[base.counts_and_times(ca,rate*t,t,ua) for t in ts]
                sb=[base.counts_and_times(cb,rate*t,t,ub) for t in ts]
                ta=np.concatenate(sa);tb=np.concatenate(sb)
                ka=np.isfinite(ta).sum(axis=1);kb=np.isfinite(tb).sum(axis=1)
                bra=base.Ensemble(m,size*len(ts));ket=base.Ensemble(m,size*len(ts));bra.rate=ket.rate=rate
                bra.w*=np.prod(np.where(np.arange(max_events)[None,:]<ka[:,None],np.tile(ja,(len(ts),1)),1.),axis=1)
                ket.w*=np.prod(np.where(np.arange(max_events)[None,:]<kb[:,None],np.tile(jb,(len(ts),1)),1.),axis=1)
                endpoints=np.repeat(ts,size)
                bra.through(ta,endpoints);ket.through(tb,endpoints)
                # Move the scalar compensation outside the common measurement
                # kernel; this supports a different final time for each row.
                values=base.centered_pairs(bra,ket,0.,rate).reshape(len(ts),size,m['norb'])
                probability=np.array([base.category_probability(ca,rate*t)*base.category_probability(cb,rate*t)*mult for t in ts])
                contributions[number,first:first+len(ts)]+=values.sum(axis=1)*(probability*np.exp(2*rate*ts))[:,None]/n
            offsets[number]=offset+size;blocks+=1
            if checkpoint is not None:
                checkpoint.save(dict(method='poisson',sector=number,bra=ca,ket=cb,samples_completed=int(offsets[number]),samples_per_sector=n),
                    contributions=contributions,offsets=offsets)
            if stop_after_batches is not None and blocks>=stop_after_batches:raise PlannedStop('Controlled checkpoint test')
    initial=np.zeros(m['norb']);initial[m['occupied']]=1
    info=dict(method='stratified-lattice-time-blocks',rate=rate,expected_jumps_per_leg=float(rate*times[-1]),
        pairs_per_sampled_stratum=n,sampled_pairs=8*n,batch=batch,time_block=time_block,max_simultaneous_pairs=batch*time_block,
        time_sampling='scaled',time_transform='cubic with exact Jacobian',normalization='known exact unitary norm 1; centered observable',
        truncation='full conditional >=4 Poisson tail; event-buffer overflow raises',max_events=max_events,
        sectors=[dict(bra=a,ket=b,multiplicity=k,analytically_zero=(a!='tail' and b!='tail' and int(a)%2!=int(b)%2)) for a,b,k in sectors])
    return initial[None,:]+contributions.sum(axis=0),info,contributions


def qm_checkpointed(m,times,nfock,checkpoint=None,stop_after_steps=None):
    if math.comb(m['norb'],m['nel'])*nfock>2_000_000:raise ValueError('QM validation space limit')
    basis=list(base.combinations(range(m['norb']),m['nel']));index={s:i for i,s in enumerate(basis)}
    occ=np.zeros((len(basis),m['norb']));rows=[];cols=[];values=[]
    for source,state in enumerate(basis):
        occ[source,list(state)]=1
        for j in range(1,m['norb']):
            if (0 in state)==(j in state):continue
            annihilate,create=(0,j) if 0 in state else (j,0)
            target=list(state);p=target.index(annihilate);target.pop(p);q=sum(x<create for x in target);target.insert(q,create)
            rows.append(index[tuple(target)]);cols.append(source);values.append((-1.)**(p+q)*m['coupling'][j])
    electronic=base.csr_matrix((values,(rows,cols)),shape=(len(basis),len(basis)))
    en=occ@np.asarray(m['energies']);w=m['frequency'];d=m['displacement'];levels=np.arange(nfock)
    diagonal=en[:,None]+w*(levels[None,:]+.5)+occ[:,0,None]*w*d*d
    h=base.kron(electronic,base.eye(nfock),format='csr')+base.diags(diagonal.ravel(),format='csr')
    x=base.diags([np.sqrt(np.arange(1,nfock)),np.sqrt(np.arange(1,nfock))],[-1,1],shape=(nfock,nfock),format='csr')
    h+=base.kron(base.diags(-w*d*occ[:,0]),x,format='csr')
    coherent=np.exp(-d*d/2+levels*np.log(d)-.5*base.gammaln(levels+1))
    psi=np.zeros((len(basis),nfock),complex);psi[index[tuple(m['occupied'])]]=coherent;psi=psi.ravel()
    initial_norm=float(np.vdot(psi,psi).real);shift=float(h.diagonal().mean())
    generator=-1j*(h-shift*base.eye(len(psi),format='csr'))
    answer=np.zeros((len(times),m['norb']));norms=np.zeros(len(times));tail=np.zeros(len(times));next_index=0
    if checkpoint is not None and checkpoint.arrays:
        z=checkpoint.arrays;psi=z['psi'];answer=z['answer'];norms=z['norms'];tail=z['tail'];next_index=int(z['next_index'])
    last=times[next_index-1] if next_index else 0.;steps=0
    for it in range(next_index,len(times)):
        t=times[it]
        if t>last:psi=base.expm_multiply(generator*(t-last),psi,traceA=0.)
        density=abs(psi.reshape(len(basis),nfock))**2
        answer[it]=density.sum(axis=1)@occ;norms[it]=density.sum();tail[it]=density[:,-8:].sum();last=t;steps+=1
        stop=stop_after_steps is not None and steps>=stop_after_steps
        if checkpoint is not None and (it%50==0 or it==len(times)-1 or stop):
            checkpoint.save(dict(method='qm',time_au=float(t),tmax=float(times[-1]),completed_observations=it+1),
                psi=psi,answer=answer,norms=norms,tail=tail,next_index=np.array(it+1))
        if stop:raise PlannedStop('Controlled QM checkpoint test')
    return answer,dict(method='qm',nfock=nfock,full_electronic_states=len(basis),wavefunction_complex_entries=len(psi),
        wavefunction_bytes=psi.nbytes,initial_norm=initial_norm,max_norm_error=float(np.max(abs(norms-1))),
        max_top8_fock_probability=float(tail.max()),propagator='full sparse many-body Fock Hamiltonian; continuous time, no splitting')


def main():
    p=argparse.ArgumentParser();p.add_argument('--method',choices=['qm','poisson'],required=True)
    p.add_argument('--norb',type=int,default=10);p.add_argument('--nel',type=int,default=5)
    p.add_argument('--tmax',type=float,default=4000);p.add_argument('--output-step',type=float,default=2.)
    p.add_argument('--pairs',type=int,default=1024);p.add_argument('--seed',type=int,default=1164001)
    p.add_argument('--batch',type=int,default=128);p.add_argument('--time-block',type=int,default=8)
    p.add_argument('--rate-scale',type=float,default=1.);p.add_argument('--max-events',type=int,default=24)
    p.add_argument('--nfock',type=int,default=384);p.add_argument('--out',type=Path,required=True)
    args=p.parse_args()
    if min(args.tmax,args.output_step,args.pairs,args.batch,args.time_block,args.rate_scale)>0 and args.max_events>=4 and args.nfock>=2:pass
    else:raise ValueError('Invalid numerical parameters')
    times=np.arange(0,args.tmax+args.output_step*.1,args.output_step)
    if abs(times[-1]-args.tmax)>1e-10:raise ValueError('tmax must be a multiple of output-step')
    m=base.model(args.norb,args.nel)
    source=Path(__file__).read_bytes();kernel=Path(base.__file__).read_bytes()
    config={k:v for k,v in vars(args).items() if k!='out'}
    config.update(model=m,source_sha256=hashlib.sha256(source).hexdigest(),kernel_sha256=hashlib.sha256(kernel).hexdigest())
    args.out.mkdir(parents=True,exist_ok=True)
    existing=args.out/'configuration.json'
    if existing.exists() and json.loads(existing.read_text(encoding='utf-8'))!=config:raise ValueError('Output directory belongs to another configuration')
    if (args.out/'metrics.json').exists():raise ValueError('Completed output exists; choose a new directory')
    existing.write_text(json.dumps(config,indent=2),encoding='utf-8')
    (args.out/'source_snapshot.py').write_bytes(source);(args.out/'kernel_snapshot.py').write_bytes(kernel)
    (args.out/'run_command.json').write_text(json.dumps([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:]],indent=2),encoding='utf-8')
    checkpoint=Checkpoint(args.out,config)
    if args.method=='qm':answer,info=qm_checkpointed(m,times,args.nfock,checkpoint)
    else:
        answer,info,contributions=poisson_blocked(m,times,args.pairs,args.seed,args.batch,args.time_block,args.rate_scale,args.max_events,checkpoint)
        np.savez_compressed(args.out/'sector_contributions.npz',time=times,contributions=contributions)
    if not np.isfinite(answer).all():raise FloatingPointError('Non-finite result')
    info.update(checkpoint.resources());info.update(model=m,seed=args.seed,source_sha256=config['source_sha256'],kernel_sha256=config['kernel_sha256'],
        output_step=args.output_step,observation_count=len(times),python=sys.version,numpy=np.__version__,scipy=base.scipy.__version__,
        platform=platform.platform(),processor=platform.processor(),min_occupation=float(answer.min()),max_occupation=float(answer.max()),
        max_particle_number_error=float(np.max(abs(answer.sum(axis=1)-m['nel']))))
    np.savez_compressed(args.out/'occupations.npz',time=times,occupations=answer)
    np.savetxt(args.out/'occupations.dat',np.c_[times,answer],header='time_au '+' '.join('orbital_'+str(i) for i in range(m['norb'])))
    (args.out/'metrics.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='COMPLETE',method=info['method'],**checkpoint.resources())),flush=True)
if __name__=='__main__':main()
