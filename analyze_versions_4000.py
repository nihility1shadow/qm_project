"""Compare retained server binaries against one QM dataset after strict data checks."""
from pathlib import Path
import csv,json,re
import numpy as np
from analyze_v113_large_4000 import validate_grid
from analyze_v114_small_4000 import metrics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'scan_versions_4000_20260907';FIG=OUT/'figures'


def bath(path):
    values=[]
    for line in path.read_text().splitlines():
        if re.match(r'^#\s*\d+\s+[+\-\d.]',line):
            row=line.lstrip('#').split()
            if len(row)==3:
                try:values.append([float(x) for x in row])
                except ValueError:pass
    return np.asarray(values)


def main():
    FIG.mkdir(exist_ok=True)
    qmpath=ROOT/'scan_v113_two_buffer_20260905/647656/ahm-qm-s10-n5.dat'
    qm=np.loadtxt(qmpath);validate_grid(qm,10,8000,.5);qmbath=bath(qmpath)
    assert qmbath.shape==(10,3)
    rows=list(csv.DictReader((OUT/'jobs.tsv').open(),delimiter='\t'))
    report={'notice':'Same physical model, requested 10,000 paths, 64 ranks, seed, and time grid. Features absent from an old binary may ignore newer controls. Retained binary identities are recorded by SHA; this is not a claim that binaries were rebuilt from the current tags.',
     'QM_job':647656,'time_range_au':[0,4000],'jobs':[]};ready=[];arrays={}
    for row in rows:
        p=OUT/row['job'];datafile=p/'ahm-sepmb-s10-n5-10000.dat'
        if not all(x.exists() for x in [datafile,p/'time.txt',p/'program.out',p/'config.txt',p/'binary.sha256']):
            report['jobs'].append({**row,'status':'awaiting complete downloaded output'});continue
        data=np.loadtxt(datafile);validate_grid(data,10,8000,.5)
        if not np.array_equal(bath(datafile),qmbath):raise ValueError(f'{p}: actual bath Hamiltonian differs from QM')
        log=(p/'program.out').read_text();header='\n'.join(x for x in datafile.read_text().splitlines() if x.startswith('#'))
        config=dict(line.split('=',1) for line in (p/'config.txt').read_text().splitlines() if '=' in line)
        for key,val in {'AHM_WC_EV':2.7,'AHM_ETA':3e-7,'AHM_DELE_EV':-17.194874968839816,'AHM_NSTEP':8000,'AHM_SEED':1158001,'SLURM_NTASKS':64,'SEP_MB_DETERMINISTIC_FOCK':0,'SEP_MB_EXACT_ORBITALS':0}.items():
            if float(config.get(key,'nan'))!=val:raise ValueError(f'{p}: config {key} differs')
        for text in ['nproc=64','seed_base=1158001','nstep=8000 dt=0.5']:
            if text not in log:raise ValueError(f'{p}: binary did not confirm {text}')
        if 'deterministic-Fock core active' in header or 'old-qm-star-poisson active' in header:raise ValueError('Unexpected solver branch')
        def find(name,text=log):
            m=re.search(re.escape(name)+r'=([\d.eE+-]+)',text);return float(m[1]) if m else None
        time=(p/'time.txt').read_text();wall=find('wall_seconds',time)
        if wall is None:raise ValueError('missing finished time record')
        peak=find('max_rank_peak_rss_kb');summed=find('sum_rank_peak_rss_kb')
        result={**row,'status':'complete and validated','binary_sha256':(p/'binary.sha256').read_text().split()[0],
          'actual_output_headers':header,'wall_seconds':wall,'allocated_rank_hours':64*wall/3600,
          'max_rank_peak_mib':None if peak is None else peak/1024,
          'sum_rank_peaks_gib':None if summed is None else summed/1024**2,**metrics(qm,data[:,4:])}
        result['minimum_unclipped_occupation']=float(data[:,4:].min())
        result['maximum_unclipped_occupation']=float(data[:,4:].max())
        result['max_absolute_orbital_error']=float(abs(data[:,4:]-qm[:,4:]).max())
        report['jobs'].append(result);ready.append(result);arrays[row['version']]=data
    if ready:
        fig,axes=plt.subplots(1,3,figsize=(13,4.6));labels=[r['version'] for r in ready];x=np.arange(len(ready))
        for ax,key,title in zip(axes,['Q_against_QM_weakest_active_orbital_window','wall_seconds','max_rank_peak_mib'],['Worst active orbital/window Q against QM','Wall time (seconds), 64 ranks','Maximum process peak RSS (MiB)']):
            vals=[r[key] if r[key] is not None else np.nan for r in ready]
            bars=ax.bar(x,vals,color='#35877e');ax.set_xticks(x,labels);ax.set_title(title,fontsize=10);ax.grid(axis='y',alpha=.15)
            if key.startswith('Q_'):ax.set_yscale('log')
            else:ax.set_ylim(0,max(vals)*1.22)
            for bar,value in zip(bars,vals):
                if np.isfinite(value):ax.text(bar.get_x()+bar.get_width()/2,value,f'{value:.2g}' if key.startswith('Q_') else f'{value:.1f}',ha='center',va='bottom',fontsize=8)
        axes[0].axhline(10,color='#b44c49',ls='--',lw=1)
        qm_profile=ROOT/'scan_v115_generality_20260906/qm_resource_validation.json'
        if qm_profile.exists():
            measured=json.loads(qm_profile.read_text())['process_peak_mib']
            axes[2].axhline(measured,color='#b44c49',ls='--',lw=1,label=f'QM (1 process): {measured:.2f} MiB')
            axes[2].set_ylim(0,max(measured,max(r['max_rank_peak_mib'] or 0 for r in ready))*1.25);axes[2].legend(fontsize=7)
        fig.suptitle('10 orbitals / 5 electrons | 4000 au | 10,000 paths per retained binary\nActual supported branches and memory fields are recorded; peak sums are not simultaneous RSS',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.87));fig.savefig(FIG/'version_comparison.png',dpi=150);plt.close(fig)
    if 'v109' in arrays and 'v114' in arrays:
        a,b=arrays['v109'],arrays['v114']
        report['v109_v114_matched_comparison']={'max_occupation_difference':float(abs(a[:,4:]-b[:,4:]).max()),'max_difference_all_columns':float(abs(a-b).max()),'wall_speed_ratio':next(r['wall_seconds'] for r in ready if r['version']=='v109')/next(r['wall_seconds'] for r in ready if r['version']=='v114'),'notice':'One run of each binary, same physical inputs/path count/ranks/seed; not an equal-error optimized benchmark or a universal speedup.'}
    if ready:
        fig,axes=plt.subplots(2,3,figsize=(13,7),sharex=True,sharey=True)
        for ax,r in zip(axes.flat,ready):
            delta=arrays[r['version']][:,4:]-qm[:,4:]
            ax.semilogy(qm[:,0],np.maximum(abs(delta).max(axis=1),1e-18),lw=.55,color='#35877e')
            ax.set_title(r['version']+' | max error over all 10 orbitals',fontsize=9)
            ax.set(xlabel='Time (a.u.)',ylabel='Absolute occupation error');ax.grid(alpha=.15)
        fig.suptitle('Complete 0–4000 au | Same 10,000 paths per binary | Errors against full many-electron grid QM\nOriginal output; common logarithmic scale; floor 1e-18 is used only for plotting zero error',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.91));fig.savefig(FIG/'version_errors_over_time.png',dpi=145);plt.close(fig)
    (OUT/'validation.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps([{'version':r['version'],'job':r['job'],'status':r['status']} for r in report['jobs']],indent=2))

if __name__=='__main__':main()
