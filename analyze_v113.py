"""Cloud v1.13 regressions, resource measurements and full 4000-au QM comparison."""
import json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'scan_v113_two_buffer_20260905';OLD=ROOT/'scan_v111_distributed_20260905'
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':150})
def folder(job):
    p=OUT/str(job)
    return p if p.exists() else OLD/str(job)
def data(job,pattern='ahm-sepmb-*.dat'):
    matches=list(folder(job).glob(pattern));assert len(matches)==1,(job,matches)
    a=np.loadtxt(matches[0]);assert np.isfinite(a).all();return a

def resources(job):
    p=folder(job);log=(p/'program.out').read_text();timing=(p/'time.txt').read_text()
    def get(pattern):
        match=re.search(pattern,log);return float(match[1]) if match else None
    mx=get(r'max_rank_peak_rss_kb=(\d+)');summed=get(r'sum_rank_peak_rss_kb=(\d+)')
    return {'job':job,'wall_seconds':float(re.search(r'wall_seconds=([\d.]+)',timing)[1]),
            'max_rank_peak_mib':mx/1024 if mx else None,'sum_rank_peak_gib':summed/1024**2 if summed else None,
            'reference_seconds':get(r'reference_seconds=([\d.]+)'),
            'sampling_seconds':get(r'sampling_seconds=([\d.]+)')}

regressions=[]
for a,b in [(647618,647652),(647619,647653),(647624,647654),(647654,647655),
            (647652,647658),(647652,647662),(647662,647663),(647621,647665),(647665,647667)]:
    x,y=data(a),data(b);assert x.shape==y.shape and np.array_equal(x[:,0],y[:,0])
    error=float(np.abs(x-y).max());assert error<2e-12,(a,b,error)
    regressions.append({'jobs':[a,b],'max_abs_all_columns':error,'max_abs_orbital_difference':float(np.abs(x[:,4:]-y[:,4:]).max())})
rejection=(folder(647664)/'program.err').read_text()
assert 'mismatched physical/solver parameters' in rejection and 'errorcode 3' in rejection

qm=data(647656,'ahm-qm-s10-n5.dat');poisson=data(647657);fock=data(647666)
assert qm.shape==poisson.shape==fock.shape and qm[-1,0]==4000
assert np.array_equal(qm[:,0],poisson[:,0]) and np.array_equal(qm[:,0],fock[:,0])
t=qm[:,0];initial=np.r_[np.ones(5),np.zeros(5)];active=[0,5,6,7,8,9]
signal=qm[:,4:]-initial
rms=lambda a:float(np.sqrt(np.mean(a*a)))
def q_ratio(a,b):return rms(a)/rms(b) if rms(b)>0 else None

def metrics(a):
    err=a[:,4:]-qm[:,4:];windows=[]
    for start in range(0,4000,100):
        ix=(t>=start)&(t<=start+100)
        per_orb=[q_ratio(signal[ix,i],err[ix,i]) for i in range(10)]
        windows.append({'start':start,'stop':start+100,'Q_active':q_ratio(signal[ix][:,active],err[ix][:,active]),
                        'Q_weakest_active':min(per_orb[i] for i in active),'Q_all_orbitals':per_orb,
                        'max_absolute_orbital_error':float(abs(err[ix]).max())})
    return {'Q_active_global':q_ratio(signal[:,active],err[:,active]),'windows':windows,
            'Q_min_active_window':min(w['Q_active'] for w in windows),
            'Q_min_active_orbital_window':min(w['Q_weakest_active'] for w in windows),
            'max_absolute_orbital_error':float(abs(err).max()),
            'max_particle_error':float(abs(a[:,4:].sum(axis=1)-5).max()),
            'minimum_occupation':float(a[:,4:].min()),'maximum_occupation':float(a[:,4:].max())}
report={'status':'4000-au 10/5 comparison complete; 30/15 4000-au production jobs running or waiting for references',
        'regressions':regressions,'cache_parameter_rejection':{'job':647664,'expected_exit':3,'passed':True},
        'resources':[resources(j) for j in [647619,647653,647624,647654,647655,647656,647657,647666,647665,647667]],
        'small_4000':{'qm_job':647656,'poisson_job':647657,'full_fock_job':647666,'active_orbitals':active,
                      'notice':'Q is against same-parameter grid QM. Poisson D3 uses 246/252 reference states; full Fock D4 uses all 252 states and 1024 oscillator levels. Neither establishes large-system accuracy.',
                      'poisson':metrics(poisson),'full_fock':metrics(fock)}}
(OUT/'cloud_validation.json').write_text(json.dumps(report,indent=2))

fig,axes=plt.subplots(5,2,figsize=(11,12),sharex=True)
for i,ax in enumerate(axes.flat):
    ax.plot(t,signal[:,i],color='#25323a',lw=1.15,label='Grid QM')
    ax.plot(t,poisson[:,4+i]-initial[i],color='#d67d32',lw=.75,alpha=.85,label='Poisson: 10,000 paths, D3/F384')
    ax.plot(t,fock[:,4+i]-initial[i],color='#258878',ls='--',lw=.7,label='Full electronic Fock reference')
    ax.set_title(f'Orbital {i}');ax.grid(alpha=.18)
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
axes[0,0].legend(fontsize=6.8)
for ax in axes[-1]:ax.set_xlabel('Time (a.u.)')
fig.suptitle('10 orbitals / 5 electrons | Complete 0-4000 a.u. data | SAME physical parameters\n'
             'All orbital occupation changes shown; no clipping or smoothing',y=.995)
fig.tight_layout(rect=(0,0,1,.965));fig.savefig(FIG/'small_4000_all_orbitals.png')
fig.savefig(FIG/'small_4000_all_orbitals_preview.png',dpi=100);plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(11,4.2))
centers=np.arange(50,4000,100)
for name,color,label in [('poisson','#d67d32','Poisson D3/F384, 10,000 paths'),('full_fock','#258878','Full electronic Fock reference')]:
    metrics_=report['small_4000'][name]
    axes[0].semilogy(centers,[w['Q_weakest_active'] for w in metrics_['windows']],color=color,label=label)
axes[0].axhline(10,color='#b34343',ls='--',lw=1,label='Diagnostic Q = 10')
axes[0].set(xlabel='100-a.u. window midpoint',ylabel='Weakest active orbital Q against grid QM')
axes[0].legend(fontsize=7.2);axes[0].grid(alpha=.18)
for arr,color,label in [(poisson,'#d67d32','Poisson'),(fock,'#258878','Full electronic Fock')]:
    axes[1].semilogy(t,np.maximum(np.max(abs(arr[:,4:]-qm[:,4:]),axis=1),1e-18),color=color,lw=.8,label=label)
axes[1].set(xlabel='Time (a.u.)',ylabel='Largest absolute orbital error');axes[1].grid(alpha=.18);axes[1].legend(fontsize=8)
fig.suptitle('4000-a.u. validation: late Poisson noise remains at 10,000 paths\n'
             'The full electronic-space reference is a small-system check, not a scalable large-system solution')
fig.tight_layout();fig.savefig(FIG/'small_4000_accuracy.png');plt.close(fig)

fig,axes=plt.subplots(2,3,figsize=(12,7.5))
for row,jobs,labels,title in [
    (0,[647619,647653],['v1.12','v1.13'],'D2: 52,656 states | 64 ranks'),
    (1,[647624,647654,647655],['v1.12 / 64','v1.13 / 64','v1.13 / 128'],'D3: 715,136 states | 64 or 128 ranks')]:
    rr=[resources(j) for j in jobs]
    for col,field,unit in [(0,'wall_seconds','Wall time (s)'),(1,'max_rank_peak_mib','Maximum rank peak (MiB)'),(2,'sum_rank_peak_gib','Sum of rank peaks (GiB)')]:
        ax=axes[row,col];vals=[r[field] for r in rr];bars=ax.bar(labels,vals,color=['#97a6b2','#32887b','#487db0'][:len(vals)],width=.6)
        for bar,value in zip(bars,vals):ax.text(bar.get_x()+bar.get_width()/2,value,f'{value:.2f}',ha='center',va='bottom',fontsize=8)
        ax.set_ylim(0,max(vals)*1.2);ax.set_title(unit);ax.grid(axis='y',alpha=.15)
        if col==0:ax.set_ylabel(title,fontsize=8)
fig.suptitle('30 orbitals / 15 electrons | Matched 10-a.u. tests | 100 paths | 384 oscillator levels\n'
             'Sum of individual peaks is not simultaneous RSS; more ranks can increase aggregate memory')
fig.tight_layout(rect=(0,0,1,.94));fig.savefig(FIG/'cloud_resources.png');plt.close(fig)
print(json.dumps({'poisson_Q_min_window':report['small_4000']['poisson']['Q_min_active_window'],
                  'poisson_Q_weakest_window':report['small_4000']['poisson']['Q_min_active_orbital_window'],
                  'full_Fock_Q_weakest_window':report['small_4000']['full_fock']['Q_min_active_orbital_window'],
                  'regression_pairs':len(regressions),'figures':str(FIG)},indent=2))
