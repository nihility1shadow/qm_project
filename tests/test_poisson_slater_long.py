"""Validate temporal batching and exact restart against the pre-existing solver."""
import sys,json,time,tempfile,contextlib,io
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
import poisson_slater_long as long
base=long.base;np=long.np

def check():
    cases=[]
    for norb,nel in [(4,2),(10,5),(30,15),(50,25)]:
        m=base.model(norb,nel);times=np.array([0.,10.,100.,4000.])
        with contextlib.redirect_stdout(io.StringIO()):
            expected,_=base.run_poisson(m,times,32,1164001,'stratified-lattice',8,time_sampling='scaled')
        for block in [1,3,8]:
            actual,_,_=long.poisson_blocked(m,times,32,1164001,batch=8,time_block=block)
            error=float(np.max(abs(expected-actual)));assert error<2e-11,(norb,nel,block,error)
            cases.append(dict(norb=norb,nel=nel,time_block=block,max_difference=error))
    m=base.model();times=np.array([0.,2.,10.,100.,1710.,4000.])
    expected,_,_=long.poisson_blocked(m,times,64,1164002,batch=16,time_block=3)
    folder=Path(tempfile.mkdtemp(prefix='slater-long-check-',dir=ROOT/'tmp'))
    with contextlib.redirect_stdout(io.StringIO()):
        c=long.Checkpoint(folder/'poisson',{'test':'poisson-resume'})
        try:long.poisson_blocked(m,times,64,1164002,16,3,checkpoint=c,stop_after_batches=3)
        except long.PlannedStop:pass
        else:raise AssertionError('Expected controlled stop')
        c=long.Checkpoint(folder/'poisson',{'test':'poisson-resume'})
        actual,_,_=long.poisson_blocked(m,times,64,1164002,16,3,checkpoint=c)
    resume_error=float(np.max(abs(expected-actual)));assert resume_error==0.
    try:long.Checkpoint(folder/'poisson',{'test':'different-configuration'})
    except ValueError:pass
    else:raise AssertionError('Mismatched checkpoint accepted')
    qt=np.arange(0.,12.,2.);qexpected,_=base.qm_oracle(m,qt,192)
    with contextlib.redirect_stdout(io.StringIO()):
        c=long.Checkpoint(folder/'qm',{'test':'qm-resume'})
        try:long.qm_checkpointed(m,qt,192,c,stop_after_steps=4)
        except long.PlannedStop:pass
        else:raise AssertionError('Expected QM stop')
        qactual,_=long.qm_checkpointed(m,qt,192,long.Checkpoint(folder/'qm',{'test':'qm-resume'}))
    qm_error=float(np.max(abs(qexpected-qactual)));assert qm_error<1e-13
    bt=np.arange(101.);start=time.perf_counter()
    with contextlib.redirect_stdout(io.StringIO()):old,_=base.run_poisson(m,bt,1024,1160101,'stratified-lattice',128,time_sampling='scaled')
    original_seconds=time.perf_counter()-start;start=time.perf_counter()
    new,_,_=long.poisson_blocked(m,bt,1024,1160101,128,8)
    blocked_seconds=time.perf_counter()-start;benchmark_error=float(np.max(abs(new-old)));assert benchmark_error<1e-13
    return dict(status='PASS',cases=cases,poisson_resume_max_difference=resume_error,qm_resume_max_difference=qm_error,
        benchmark=dict(tmax=100,pairs_per_stratum=1024,original_seconds=original_seconds,blocked_seconds=blocked_seconds,
            speedup=original_seconds/blocked_seconds,max_difference=benchmark_error),
        scope='Batching/restart equivalence; not a claim of 4000-au trajectory accuracy')
if __name__=='__main__':
    result=check();dest=ROOT/'scan_slater_stream_4000au_20260910';dest.mkdir(exist_ok=True)
    (dest/'implementation_validation.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result),flush=True)
