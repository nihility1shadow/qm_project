"""4000-au million-path small-system validation against existing exact grid QM."""
import csv,json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_v113_large_4000 import windows,ratio,validate_grid
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'scan_v114_signed_csr_20260906'
OLD=ROOT/'scan_v113_two_buffer_20260905';FIG=OUT/'figures';FIG.mkdir(exist_ok=True)

def metrics(qm,values):
    initial=np.r_[np.ones(5),np.zeros(5)];active=[0,5,6,7,8,9]
    signal=qm[:,4:]-initial;error=values-qm[:,4:]
    ww=windows(qm[:,0],signal,error,active)
    return {'Q_against_QM_global_active':ratio(signal[:,active],error[:,active]),
            'Q_against_QM_min_active_window':min(w['Q_active'] for w in ww),
            'Q_against_QM_weakest_active_orbital_window':min(w['Q_weakest_active'] for w in ww),
            'windows':ww,'max_orbital_error':float(abs(error).max()),
            'max_particle_error':float(abs(values.sum(axis=1)-5).max()),
            'minimum_occupation':float(values.min()),'maximum_occupation':float(values.max())}


def main():
    qm=np.loadtxt(OLD/'647656/ahm-qm-s10-n5.dat');validate_grid(qm,10,8000,.5)
    entries=list(csv.DictReader((OUT/'small_4000_jobs.tsv').open(),delimiter='\t'))
    report={'status':'awaiting million-path results','QM_job':647656,'orbitals':10,'electrons':5,'time_range_au':[0,4000],
            'definition':'Q = RMS(QM occupation - initial occupation) / RMS(Poisson occupation - QM occupation)',
            'notice':'The D3 reference includes 246 of 252 electronic determinants here; high small-system Q alone does not establish scalable 30/15 accuracy.',
            'jobs':[], 'resource_notice':'Cache consumers only; reference producer cost reported separately. Sum of rank peaks is not simultaneous RSS; allocated rank-hours are not measured CPU time.'}
    complete=[]
    for row in entries:
        if row['role']!='sample':continue
        folder=OUT/row['job'];path=folder/'ahm-sepmb-s10-n5-1000000.dat'
        if not all(p.exists() for p in [path,folder/'time.txt',folder/'program.out',folder/'config.txt']):
            report['jobs'].append({'job':row['job'],'status':'awaiting complete downloaded output'});continue
        config=dict(line.split('=',1) for line in (folder/'config.txt').read_text().splitlines() if '=' in line)
        for name,value in {'AHM_WC_EV':2.7,'AHM_ETA':3e-7,'AHM_DELE_EV':-17.194874968839816,'AHM_NSTEP':8000,
                           'AHM_SEED':int(row['seed']),'SEP_MB_REFERENCE_FOCK_STATES':384,'SEP_MB_REFERENCE_DISTANCE':3,
                           'SLURM_NTASKS':64}.items():
            if float(config.get(name,'nan'))!=value:raise ValueError(f'{folder}: {name} differs')
        log=(folder/'program.out').read_text()
        if '#SEP_MB_TIMING sampling_seconds=' not in log:raise ValueError(f'{folder}: no completion marker')
        if not re.search(r'selected_distance=3 active=1 states=246 max_states=\d+ fock_states=384',log):
            raise ValueError(f'{folder}: actual reference differs')
        values=np.loadtxt(path);validate_grid(values,10,8000,.5)
        timing=(folder/'time.txt').read_text()
        wall=float(re.search(r'wall_seconds=([\d.]+)',timing)[1])
        resource={'wall_seconds':wall,'ranks':64,'allocated_rank_hours':wall*64/3600,
                  'max_rank_peak_mib':int(re.search(r'max_rank_peak_rss_kb=(\d+)',log)[1])/1024,
                  'sum_rank_peaks_gib':int(re.search(r'sum_rank_peak_rss_kb=(\d+)',log)[1])/1024**2}
        report['jobs'].append({'job':row['job'],'seed':row['seed'],'status':'complete','resources':resource,**metrics(qm,values[:,4:])})
        complete.append((row,values))
    if len(complete)==3:
        if len({row['seed'] for row,_ in complete})!=3:raise ValueError('independent seeds required')
        values=np.stack([a[:,4:] for _,a in complete]);mean=values.mean(axis=0)
        report['status']='three complete 4000-au million-path runs validated against QM'
        report['three_seed_mean']=metrics(qm,mean)
        initial=np.r_[np.ones(5),np.zeros(5)];t=qm[:,0]
        fig,axes=plt.subplots(5,2,figsize=(11,12),sharex=True)
        for orbital,ax in enumerate(axes.flat):
            ax.plot(t,qm[:,4+orbital]-initial[orbital],color='#263740',lw=1.1,label='Grid QM')
            for (row,data),color in zip(complete,['#bb936b','#95aa82','#9295b4']):
                ax.plot(t,data[:,4+orbital]-initial[orbital],color=color,alpha=.6,lw=.55,label=f"Seed {row['seed']}")
            ax.plot(t,mean[:,orbital]-initial[orbital],color='#147e83',lw=.8,label='3-seed mean')
            ax.set_title(f'Orbital {orbital}');ax.grid(alpha=.15)
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
        axes[0,0].legend(fontsize=6)
        for ax in axes[-1]:ax.set_xlabel('Time (a.u.)')
        fig.suptitle('10 orbitals / 5 electrons | 3 independent runs, 1 million paths each\nSame physical parameters; unmodified output (residual measured every 4 steps, linearly interpolated)')
        fig.tight_layout(rect=(0,0,1,.96));fig.savefig(FIG/'small_million_4000_all_orbitals.png',dpi=145);plt.close(fig)
        fig,axes=plt.subplots(2,2,figsize=(12,7.6))
        old=np.loadtxt(OLD/'647657/ahm-sepmb-s10-n5-10000.dat')
        oldmetrics=metrics(qm,old[:,4:]);report['old_10000_paths']=oldmetrics
        colors=['#b07042','#538349','#8471aa']
        for row,color in zip(report['jobs'],colors):
            ww=row['windows'];centers=[(w['start']+w['stop'])/2 for w in ww]
            axes[0,0].plot(centers,[w['Q_weakest_active'] for w in ww],color=color,label=f"Seed {row['seed']}")
        ww=oldmetrics['windows']
        axes[0,0].plot(centers,[w['Q_weakest_active'] for w in ww],color='#a8a8a8',ls='--',label='Earlier 10,000 paths')
        axes[0,0].axhline(10,color='#bb4040',lw=1,ls=':');axes[0,0].set_yscale('log')
        axes[0,0].set(title='Worst active-orbital Q in each 100-au window',xlabel='Time (a.u.)',ylabel='Q against grid QM')
        axes[0,0].legend(fontsize=7)
        for (row,data),color in zip(complete,colors):
            axes[0,1].plot(qm[:,0],data[:,4]-qm[:,4],color=color,lw=.6,label=row['seed'])
        axes[0,1].set(title='Impurity occupation error against QM',xlabel='Time (a.u.)',ylabel='Absolute occupation difference')
        labels=[r['seed'] for r in report['jobs']];xs=np.arange(3)
        vals=[r['resources']['wall_seconds']/60 for r in report['jobs']]
        bars=axes[1,0].bar(xs,vals,color=colors)
        for bar,val in zip(bars,vals):axes[1,0].text(bar.get_x()+bar.get_width()/2,val,f'{val:.2f}',ha='center',va='bottom')
        axes[1,0].axhline(3287.21/60,color='#444',ls='--',label='QM: 54.79 min, 1 rank')
        axes[1,0].set(xticks=xs,xticklabels=labels,ylabel='Elapsed minutes',title='Each Poisson run: 64 ranks, cached reference',ylim=(0,75));axes[1,0].legend(fontsize=8)
        vals=[r['resources']['max_rank_peak_mib'] for r in report['jobs']]
        bars=axes[1,1].bar(xs,vals,color=colors)
        for bar,val in zip(bars,vals):axes[1,1].text(bar.get_x()+bar.get_width()/2,val,f'{val:.2f}',ha='center',va='bottom')
        axes[1,1].set(xticks=xs,xticklabels=labels,ylabel='MiB',title='Maximum rank peak RSS (sampling jobs)',ylim=(0,max(vals)*1.2))
        for ax in axes.flat:ax.grid(alpha=.14)
        fig.suptitle('10 orbitals / 5 electrons | 0–4000 au | 1 million paths per independent run\nQ uses occupation change, not the large initial occupation baseline; D3 reference has 246/252 states',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.93));fig.savefig(FIG/'small_million_4000_validation.png',dpi=150);plt.close(fig)
    (OUT/'small_4000_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps({'status':report['status'],'jobs':[{'job':r['job'],'status':r['status']} for r in report['jobs']]},indent=2))

if __name__=='__main__':main()
