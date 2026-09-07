"""Validate downloaded 50-orbital jobs; never label a bounded reference as exact QM."""
from pathlib import Path
from math import comb
import csv,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_v113_large_4000 import validate_grid,load_reference
from analyze_midpoint_probability import midpoint_count
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'scan_v115_generality_20260906';FIG=OUT/'figures'
BINARY_SHA='4354060f4de44e5ef7a644b35fb401111385e7e934f202c7b283589d685edda2'


def resources(folder,reference_active=True):
    log=(folder/'program.out').read_text();timing=(folder/'time.txt').read_text()
    def field(name,text=log):
        match=re.search(re.escape(name)+r'=([0-9.eE+-]+)',text)
        if not match:raise ValueError(f'{folder}: missing {name}')
        return float(match[1])
    wall=field('wall_seconds',timing);ranks=int(field('ranks'))
    return {'wall_seconds':wall,'ranks':ranks,'reference_seconds':field('reference_seconds') if reference_active else None,
            'sampling_seconds':field('sampling_seconds'),'max_rank_peak_mib':field('max_rank_peak_rss_kb')/1024,
            'sum_rank_peaks_gib':field('sum_rank_peak_rss_kb')/1024**2,'allocated_rank_hours':wall*ranks/3600}


def read_run(row):
    folder=OUT/row['job'];norb,nel,steps,paths,depth,fock,ranks,seed=[int(row[k]) for k in
        ['orbitals','electrons','steps','paths','distance','fock_levels','ranks','seed']]
    config=dict(line.split('=',1) for line in (folder/'config.txt').read_text().splitlines() if '=' in line)
    expected={'AHM_WC_EV':2.7,'AHM_ETA':3e-7,'AHM_DELE_EV':-17.194874968839816,'AHM_NSTEP':steps,
              'AHM_SEED':seed,'SEP_MB_PATH_LOCAL_BASIS':1,'SEP_MB_DETERMINISTIC_FOCK':0,
              'SEP_MB_EXACT_ORBITALS':0,'SEP_MB_REFERENCE_DISTANCE':depth,
              'SEP_MB_REFERENCE_FOCK_STATES':fock,'SLURM_NTASKS':ranks}
    for key,value in expected.items():
        if float(config.get(key,'nan'))!=value:raise ValueError(f'{folder}: mismatched {key}')
    if (folder/'binary.sha256').read_text().split()[0]!=BINARY_SHA:raise ValueError('unexpected binary')
    data=np.loadtxt(folder/f'ahm-sepmb-s{norb}-n{nel}-{paths}.dat');validate_grid(data,norb,steps,.5)
    log=(folder/'program.out').read_text()
    for marker in [f'nproc={ranks}',f'seed_base={seed}',f'nstep={steps} dt=0.5','selected=stochastic-poisson','full_determinant_basis=skipped']:
        if marker not in log:raise ValueError(f'{folder}: actual solver missing {marker}')
    match=re.search(r'requested_distance=(-?\d+) selected_distance=(-?\d+) active=(\d+) states=(\d+) max_states=\d+ fock_states=(\d+)',log)
    if not match or int(match[2])!=depth:raise ValueError('reference distance silently changed')
    states=0 if depth<0 else sum(comb(nel-1,d)*(comb(norb-nel,d)+comb(norb-nel,d+1)) for d in range(depth+1))
    if int(match[3])!=(depth>=0) or int(match[4])!=states or int(match[5])!=fock:raise ValueError('reference geometry mismatch')
    reference=None;stats={}
    header='\n'.join(line for line in (folder/f'ahm-sepmb-s{norb}-n{nel}-{paths}.dat').read_text().splitlines() if line.startswith('#'))
    sampling=next(line for line in header.splitlines() if line.startswith('#sampling:'))
    if 'stratify_forward_steps=1' not in sampling:raise ValueError('unexpected sampling branch')
    probability=float(re.search(r'jump_probability=([0-9.eE+-]+)',sampling)[1])
    effective=midpoint_count(paths,probability)/paths
    stats['fixed_step_probability_audit']={'nominal':probability,'actual_midpoint_marginal':effective,'relative_event_probability_error':effective/probability-1,'notice':'Event-probability rounding error, NOT the physical occupation error. This low-path pilot cannot establish sampling accuracy.'}
    if depth>=0 and row['reference_job']=='-':
        reference,reference_stats=load_reference(folder/'reference-observables.dat',norb,nel,steps,.5)
        stats.update(reference_stats)
        stats['max_pilot_difference_to_reference']=float(abs(data[:,4:]-reference[:,4:]).max())
    initial=np.r_[np.ones(nel),np.zeros(norb-nel)]
    return data,reference,{'job':row['job'],'role':row['role'],'status':'complete data validated; accuracy NOT accepted',
        'accuracy_status':'Low-path pilot: visible nonphysical fluctuations and fixed-step event rounding; not a high-accuracy production result',
        'reference_states':states,'resources':resources(folder,depth>=0),**stats,
        'maximum_occupation_change':float(abs(data[:,4:]-initial).max()),
        'minimum_occupation':float(data[:,4:].min()),'maximum_occupation':float(data[:,4:].max()),
        'normalized_particle_sum_error':float(abs(data[:,4:].sum(axis=1)-nel).max())}


def main():
    rows=list(csv.DictReader((OUT/'long_4000_jobs.tsv').open(),delimiter='\t'))
    report={'notice':'No exact full-space QM is available at 50/25. Reference convergence and pilot output differences are not true QM fitting. Final particle conservation is constrained by normalization.',
            'resource_notice':'Sum of rank peaks is not simultaneous RSS. Allocated rank-hours are not measured CPU time.',
            'time_range_au':[0,4000],'runs':[]};references={}
    for row in rows:
        folder=OUT/row['job']
        required=[folder/name for name in ['config.txt','time.txt','program.out','binary.sha256',f"ahm-sepmb-s{row['orbitals']}-n{row['electrons']}-{row['paths']}.dat"]]
        if not all(p.exists() for p in required):
            report['runs'].append({'job':row['job'],'role':row['role'],'status':'awaiting complete downloaded output'});continue
        data,ref,stats=read_run(row);report['runs'].append(stats)
        if ref is not None:references[row['role']]=ref
        initial=np.r_[np.ones(int(row['electrons'])),np.zeros(int(row['orbitals'])-int(row['electrons']))]
        fig,axes=plt.subplots(10,5,figsize=(16,22),sharex=True)
        for orbital,ax in enumerate(axes.flat):
            ax.plot(data[:,0],data[:,4+orbital]-initial[orbital],color='#238378',lw=.7,label=f"{row['paths']}-path pilot")
            if ref is not None:ax.plot(ref[:,0],ref[:,4+orbital]-initial[orbital],color='#9a7156',lw=.6,label='Bounded reference')
            ax.set_title(f'Orbital {orbital}',fontsize=8);ax.grid(alpha=.15)
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
        axes[0,0].legend(fontsize=6)
        for ax in axes[-1]:ax.set_xlabel('Time (a.u.)')
        fig.suptitle(f"50 orbitals / {row['electrons']} electrons | 0–4000 au | job {row['job']}\n{row['role']}; low-path pilot with rounded jump probability, NOT a QM accuracy proof",fontsize=12)
        fig.tight_layout(rect=(0,0,1,.96));fig.savefig(FIG/f"job{row['job']}_all_orbitals.png",dpi=120);plt.close(fig)
    comparisons=[]
    for a,b in [('reference_D1_F384','reference_D2_F384'),('reference_D2_F384','reference_D2_F512')]:
        if a not in references or b not in references:continue
        difference=references[a][:,4:]-references[b][:,4:]
        comparisons.append({'pair':[a,b],'max_orbital_difference':float(abs(difference).max()),'rms_orbital_difference':float(np.sqrt(np.mean(difference**2)))})
        fig,ax=plt.subplots(figsize=(10,4));ax.plot(references[a][:,0],abs(difference).max(axis=1),color='#238378')
        ax.set(xlabel='Time (a.u.)',ylabel='Maximum orbital difference',title=f'{a} vs {b}\nReference convergence test; shared bias is not excluded');ax.grid(alpha=.15)
        fig.tight_layout();fig.savefig(FIG/f'{a}_vs_{b}.png',dpi=150);plt.close(fig)
    report['reference_comparisons']=comparisons
    (OUT/'long_4000_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps({'runs':report['runs'],'reference_comparisons':comparisons},indent=2))

if __name__=='__main__':main()
