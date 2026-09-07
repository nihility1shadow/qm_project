"""Audit full-time sampling comparisons and actual process resources; no full-QM claim."""
from pathlib import Path
import csv,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_v113_large_4000 import validate_grid,load_reference,windows
from analyze_versions_4000 import bath
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'scan_rate_4000_20260907';FIG=OUT/'figures'
SHA='4354060f4de44e5ef7a644b35fb401111385e7e934f202c7b283589d685edda2'


def profile(folder,ranks):
    files=list(folder.glob('process-[0-9]*.txt'))
    expected={f'process-{i}.txt' for i in range(ranks)}
    if {p.name for p in files}!=expected:raise ValueError('Incomplete per-process records')
    records=np.array([np.loadtxt(folder/f'process-{i}.txt') for i in range(ranks)])
    if records.shape!=(ranks,4) or not np.isfinite(records).all() or (records<0).any():raise ValueError('Invalid process resources')
    summary={k:float(v) for k,v in re.findall(r'(\w+)=([0-9.eE+-]+)',(folder/'process-resources.txt').read_text())}
    actual={'ranks':ranks,'user_cpu_seconds_sum':records[:,0].sum(),'system_cpu_seconds_sum':records[:,1].sum(),'max_process_wall_seconds':records[:,2].max(),'max_rank_peak_rss_kb':records[:,3].max(),'sum_rank_peak_rss_kb':records[:,3].sum()}
    for key,value in actual.items():
        if key not in summary or abs(summary[key]-value)>1e-5:raise ValueError(f'Wrong aggregate {key}')
    actual.update({'max_rank_peak_mib':actual['max_rank_peak_rss_kb']/1024,'sum_rank_peaks_gib':actual['sum_rank_peak_rss_kb']/1024**2,'measured_cpu_hours':float(records[:,:2].sum()/3600)})
    return {key:float(value) for key,value in actual.items()}


def main():
    FIG.mkdir(exist_ok=True)
    rows=list(csv.DictReader((OUT/'sampling_jobs.tsv').open(),delimiter='\t'))
    reference,_=load_reference(ROOT/'scan_v113_two_buffer_20260905/647702/reference-observables.dat',30,15,8000,.5)
    expected_bath=bath(ROOT/'scan_v113_two_buffer_20260905/647704/ahm-sepmb-s30-n15-1000000.dat')
    initial=np.r_[np.ones(15),np.zeros(15)];active=[0,*range(15,30)]
    report={'scope':'30 orbitals / 15 electrons, 0–4000 au, 100,000 paths, 128 ranks, cached D3/F384 reference.',
      'notice':'Signal/reference-difference ratios are NOT Q against full-space QM. Two seeds per mode do not establish confidence bounds. Peak sums are not simultaneous RSS.',
      'supersedes_cancelled_array':647806,'runs':[]};ready=[]
    for row in rows:
        job=row['array_job']+'_'+row['task'];folder=OUT/job;path=folder/'ahm-sepmb-s30-n15-100000.dat'
        if not all((folder/name).exists() for name in ['time.txt','program.out','config.txt','binary.sha256','process-resources.txt']) or not path.exists():
            report['runs'].append({**row,'job':job,'status':'awaiting complete downloaded output'});continue
        config=dict(line.split('=',1) for line in (folder/'config.txt').read_text().splitlines() if '=' in line)
        expected={'AHM_WC_EV':2.7,'AHM_ETA':3e-7,'AHM_DELE_EV':-17.194874968839816,'AHM_NSTEP':8000,'AHM_SEED':int(row['seed']),'SLURM_NTASKS':128,
          'SEP_MB_REFERENCE_DISTANCE':3,'SEP_MB_REFERENCE_FOCK_STATES':384,'SEP_MB_RATE_SCALE':float(row['rate_scale']),'SEP_MB_STRATIFY_FORWARD_STEPS':int(row['forward_steps']),
          'SEP_MB_STRATIFY_FORWARD_COUNT':1,'SEP_MB_PATH_LOCAL_BASIS':1,'SEP_MB_DETERMINISTIC_FOCK':0,'SEP_MB_EXACT_ORBITALS':0}
        for key,value in expected.items():
            if float(config.get(key,'nan'))!=value:raise ValueError(f'{job}: wrong {key}')
        if (folder/'binary.sha256').read_text().split()[0]!=SHA:raise ValueError('Wrong binary')
        data=np.loadtxt(path);validate_grid(data,30,8000,.5)
        if not np.array_equal(bath(path),expected_bath):raise ValueError('Actual bath differs')
        log=(folder/'program.out').read_text()
        for marker in ['nproc=128',f"seed_base={row['seed']}",'nstep=8000 dt=0.5','full_determinant_basis=skipped','selected_distance=3 active=1 states=715136','fock_states=384','parameter_key_and_checksum=verified','loaded=1 rows=8001','#SEP_MB_TIMING sampling_seconds=']:
            if marker not in log:raise ValueError(f'{job}: missing actual marker {marker}')
        header='\n'.join(line for line in path.read_text().splitlines() if line.startswith('#'))
        sampling=next(line for line in header.splitlines() if line.startswith('#sampling:'))
        if f"stratify_forward_steps={row['forward_steps']}" not in sampling:raise ValueError('Actual sampling mode differs')
        rate=float(re.search(r'rate_scale=([0-9.eE+-]+)',sampling)[1])
        if rate!=float(row['rate_scale']):raise ValueError('Actual rate differs')
        difference=data[:,4:]-reference[:,4:];ww=windows(data[:,0],reference[:,4:]-initial,difference,active)
        stats={**row,'job':job,'status':'complete data validated; screening result only','actual_sampling_header':sampling,'resources':profile(folder,128),
          'wall_seconds':float(re.search(r'wall_seconds=([0-9.]+)',(folder/'time.txt').read_text())[1]),
          'max_difference_to_reference':float(abs(difference).max()),'minimum_unclipped_occupation':float(data[:,4:].min()),'maximum_unclipped_occupation':float(data[:,4:].max()),
          'worst_active_window_ratio':min(w['Q_active'] for w in ww if w['Q_active'] is not None),
          'worst_active_orbital_window_ratio':min(w['Q_weakest_active'] for w in ww if w['Q_weakest_active'] is not None),'windows':ww}
        report['runs'].append(stats);ready.append((stats,difference))
    if ready:
        fig,axes=plt.subplots(2,2,figsize=(12,8))
        for stats,difference in ready:
            label=f"steps={stats['forward_steps']}, rate={stats['rate_scale']}, seed={stats['seed']}"
            axes[0,0].semilogy(reference[:,0],np.maximum(abs(difference).max(axis=1),1e-18),lw=.55,label=label)
            axes[0,1].semilogy(np.arange(50,4000,100),[w['Q_weakest_active'] for w in stats['windows']],lw=1,label=label)
        axes[0,0].set(title='Maximum difference over all 30 orbitals',xlabel='Time (a.u.)',ylabel='Absolute difference to bounded reference')
        axes[0,1].axhline(10,ls='--',color='#a95046');axes[0,1].set(title='Weakest active orbital in each 100-au window',xlabel='Window midpoint (a.u.)',ylabel='Signal / difference ratio (NOT full-QM Q)')
        labels=[stats['job'].split('_')[-1] for stats,_ in ready];x=np.arange(len(ready))
        for ax,key,title in [(axes[1,0],'measured_cpu_hours','Measured sum of process CPU hours'),(axes[1,1],'max_rank_peak_mib','Measured maximum process peak RSS (MiB)')]:
            ax.bar(x,[s['resources'][key] for s,_ in ready],color='#35877e');ax.set(xticks=x,xticklabels=labels,title=title,xlabel='Array task index')
        axes[0,0].legend(fontsize=7)
        for ax in axes.flat:ax.grid(alpha=.15)
        fig.suptitle('Sampling mode screen | 30/15 | Complete 4000 au | 100,000 paths per run\nOnly finished cases shown; two seeds per setting are a screen, not an accuracy certificate',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.92));fig.savefig(FIG/'sampling_comparison.png',dpi=145);plt.close(fig)
    (OUT/'sampling_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps([{k:v for k,v in row.items() if k not in ['windows','actual_sampling_header']} for row in report['runs']],indent=2))

if __name__=='__main__':main()
