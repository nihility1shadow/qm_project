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
            'jobs':[]}
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
        report['jobs'].append({'job':row['job'],'seed':row['seed'],'status':'complete',**metrics(qm,values[:,4:])})
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
        fig.suptitle('10 orbitals / 5 electrons | 3 independent runs, 1 million paths each\nSame physical parameters as grid QM; no clipping or smoothing')
        fig.tight_layout(rect=(0,0,1,.96));fig.savefig(FIG/'small_million_4000_all_orbitals.png',dpi=145);plt.close(fig)
    (OUT/'small_4000_validation.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps({'status':report['status'],'jobs':[{'job':r['job'],'status':r['status']} for r in report['jobs']]},indent=2))

if __name__=='__main__':main()
