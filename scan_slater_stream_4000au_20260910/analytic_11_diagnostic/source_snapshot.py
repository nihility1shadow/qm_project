"""Exact (one bra jump, one ket jump) sector of the centered Poisson expansion.

This is a diagnostic second-order truncation, not the full Poisson result.
It sums the coherent-state oscillator number distribution, with a tail bound.
"""
import sys,json,time,hashlib
from pathlib import Path
import poisson_slater_stream as base
np=base.np

def one_jump_sector(m,times,tolerance=1e-15):
    lam=m['displacement']**2
    cutoff=int(base.poisson.isf(tolerance,lam));levels=np.arange(cutoff+1)
    probabilities=base.poisson.pmf(levels,lam);tail=float(base.poisson.sf(cutoff,lam))
    result=np.zeros((len(times),m['norb']));bounds=[]
    for j in range(m['nel'],m['norb']):
        gap=m['energies'][j]-m['energies'][0]
        if gap<=0:raise ValueError('Tail bound assumes positive unoccupied-to-impurity gap')
        omega=gap+m['frequency']*levels
        for first in range(0,len(times),128):
            ts=times[first:first+128]
            result[first:first+len(ts),j]=4*m['coupling'][j]**2*(np.sin(ts[:,None]*omega[None,:]/2)**2/omega[None,:]**2)@probabilities
        bounds.append(4*m['coupling'][j]**2*tail/(gap+m['frequency']*(cutoff+1))**2)
    result[:,0]=-result[:,m['nel']:].sum(axis=1)
    return result,dict(sector='bra=1, ket=1 only; higher sectors omitted',oscillator_levels=cutoff+1,
        oscillator_probability_tail=tail,maximum_occupation_tail_bound=float(sum(bounds)),
        approximation='Second order in hopping. Not the full untruncated Poisson estimator; no QM input.')

def validate(m):
    from scipy.special import roots_legendre
    checks=[];r=float(np.linalg.norm(m['coupling']))
    for t,order in [(1.,24),(10.,32),(50.,64),(100.,96)]:
        x,w=roots_legendre(order);s=(x+1)*t/2;weights=w*t/2
        sa,sb=np.meshgrid(s,s,indexing='ij');bra=base.Ensemble(m,order**2);ket=base.Ensemble(m,order**2)
        bra.rate=ket.rate=r
        bra.through(np.c_[sa.ravel(),np.full(order**2,np.inf)],t)
        ket.through(np.c_[sb.ravel(),np.full(order**2,np.inf)],t)
        sampled=base.centered_pairs(bra,ket,0.,r)
        direct=(sampled*(np.outer(weights,weights).ravel()*r*r)[:,None]).sum(axis=0)
        exact,_=one_jump_sector(m,np.array([t]));error=float(np.max(abs(direct-exact[0])))
        assert error<2e-18,(t,error)
        checks.append(dict(time=t,quadrature_order=order,max_difference=error))
    return checks

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1];dest=root/'scan_slater_stream_4000au_20260910'/'analytic_11_diagnostic';dest.mkdir(exist_ok=True)
    m=base.model();checks=validate(m);times=np.arange(0.,4002.,2.)
    started=time.perf_counter();delta,info=one_jump_sector(m,times);elapsed=time.perf_counter()-started
    initial=np.zeros(m['norb']);initial[m['occupied']]=1;answer=initial[None,:]+delta
    np.savez_compressed(dest/'occupations.npz',time=times,occupations=answer,sector_contribution=delta)
    np.savetxt(dest/'occupations.dat',np.c_[times,answer],header='DIAGNOSTIC: only exact (1,1) sector, higher orders omitted; time_au orbital_0..9')
    info.update(model=m,wall_seconds=elapsed,quadrature_checks=checks,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (dest/'source_snapshot.py').write_bytes(Path(__file__).read_bytes())
    (dest/'diagnostic.json').write_text(json.dumps(info,indent=2),encoding='utf-8');print(json.dumps(info),flush=True)
