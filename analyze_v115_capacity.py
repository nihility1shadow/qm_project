from pathlib import Path
from math import comb
import json, subprocess, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('scan_v115_generality_20260906');out.mkdir(exist_ok=True);(out/'figures').mkdir(exist_ok=True)
rows=[]
for norb,nel in [(10,5),(30,15),(50,25),(50,17),(64,32)]:
 for depth in range(4):
  a,b=nel-1,norb-nel
  states=sum(comb(a,d)*(comb(b,d)+comb(b,d+1)) for d in range(depth+1))
  edges=2*b*sum(comb(a,d)*comb(b,d) for d in range(depth+1))
  rows.append(dict(orbitals=norb,electrons=nel,depth=depth,full_determinants=comb(norb,nel),states=states,directed_edges=edges,
   one_QM_wavefunction_bytes=comb(norb,nel)*1024*16,two_reference_arrays_bytes=2*384*states*16,
   signed_graph_bytes_per_rank=8*(states+1)+4*edges))
report={'notice':'Exact combinatorial storage estimates, NOT measured RSS. QM wavefunction excludes PES/work/basis. Reference arrays exclude halo, graph copies, metadata and MPI runtime.', 'rows':rows}
(out/'capacity_estimates.json').write_text(json.dumps(report,indent=2))
fig,axes=plt.subplots(1,2,figsize=(11,4.8))
labels=['10 / 5','30 / 15','50 / 25'];x=np.arange(3)
for depth,color,offset in [(1,'#9baaaa',-.23),(2,'#298378',0),(3,'#c07d4c',.23)]:
 vals=[next(r for r in rows if (r['orbitals'],r['electrons'],r['depth'])==(*p,depth))['two_reference_arrays_bytes']/1024**3 for p in [(10,5),(30,15),(50,25)]]
 axes[0].bar(x+offset,vals,width=.22,color=color,label=f'D{depth}')
axes[0].set(xticks=x,xticklabels=labels,yscale='log',ylabel='GiB, estimated',title='Two complex reference arrays, F384');axes[0].legend()
vals=[next(r for r in rows if (r['orbitals'],r['electrons'],r['depth'])==(*p,2))['one_QM_wavefunction_bytes']/1024**3 for p in [(10,5),(30,15),(50,25)]]
axes[1].bar(x,vals,color='#7e879b',width=.6)
for ix,val in enumerate(vals):axes[1].text(ix,val*1.35,f'{val:.3g}',ha='center')
axes[1].set(xticks=x,xticklabels=labels,yscale='log',ylabel='GiB, estimated',title='One full grid-QM wavefunction, 1024 points',ylim=(1e-3,max(vals)*100))
for ax in axes:ax.set_xlabel('Orbitals / electrons');ax.grid(axis='y',alpha=.15)
fig.suptitle('Storage scaling | Estimates only, not measured process memory\nReference graph copies and other buffers are excluded from the left panel',fontsize=11)
fig.tight_layout(rect=(0,0,1,.9));fig.savefig(out/'figures/storage_scaling_estimates.png',dpi=150);plt.close(fig)
# Isolate the actual Bernoulli time-discretization used by the path sampler.
# H0=0, H1=g sigma_x; analytical expectation, not a sampled physical AHM run.
g=np.sqrt(3e-7);T=4000.;toy=[]
for dt in [1.,.5,.25,.125]:
 n=round(T/dt);eig=np.array([g,-g]);v=np.array([[1,1],[1,-1]])/np.sqrt(2)
 u=(v*(1-1j*eig*dt)**n)@v.T
 exact=(v*np.exp(-1j*eig*T))@v.T
 psi=u[:,0];truth=exact[:,0];norm=float(np.vdot(psi,psi).real)
 toy.append({'dt':dt,'steps':n,'raw_norm':norm,'raw_norm_error':abs(norm-1),
             'raw_occupation_error':float(abs(abs(psi)**2-abs(truth)**2).max()),
             'normalized_occupation_error':float(abs(abs(psi)**2/norm-abs(truth)**2).max()),
             'propagator_max_error':float(abs(u-exact).max())})
(out/'bernoulli_time_step_audit.json').write_text(json.dumps({'model':'H0=0, H1=sqrt(eta)*sigma_x; analytical diagnostic only',
 'identity':'(1+r*dt)[(1-p)I+p*(-i H1/r)]=I-i H1 dt, p=r*dt/(1+r*dt)',
 'warning':'Finite-dt expected propagator differs from exp(-i H1 T). Normalization can hide leading norm drift. This toy is not a bound on physical AHM error.', 'rows':toy},indent=2))
tags=subprocess.check_output(['git','tag','--list'],text=True).splitlines();catalog=[]
for tag in tags:
 if not(tag.startswith(('v0.92','v0.93','v0.94','v0.96','v0.98','v0.99','v1.'))):continue
 src=subprocess.check_output(['git','show',f'{tag}:ahm-mb-sep.cpp'],text=True,encoding='utf-8')
 commit=subprocess.check_output(['git','rev-parse',f'{tag}^{{commit}}'],text=True).strip()
 catalog.append({'tag':tag,'commit':commit,'path_local_option':'SEP_MB_PATH_LOCAL_BASIS' in src,
  'bounded_reference_option':'SEP_MB_REFERENCE_DISTANCE' in src,'mpi_reference_option':'SEP_MB_REFERENCE_MPI' in src,
  'rqmc_option':'SEP_MB_RQMC_REPLICATES' in src,'signed_csr_option':'SEP_MB_REFERENCE_SIGNED_CSR' in src})
(out/'release_feature_inventory.json').write_text(json.dumps({'notice':'Code-feature inventory, not a matched runtime benchmark. Presence of a flag does not prove correctness or scalability.','releases':catalog},indent=2))
print(json.dumps({'capacity50_D2':next(r for r in rows if (r['orbitals'],r['electrons'],r['depth'])==(50,25,2)), 'dt_diagnostic':toy,'release_count':len(catalog)},indent=2))
