"""Strict 30/15, 0-4000-au analysis. Repeatability is never labelled QM fitting."""
import argparse,csv,json,re,struct
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
DEFAULT=ROOT/'scan_v113_two_buffer_20260905'


def load_reference(path,norb,nel,steps,dt):
    lines=path.read_text().splitlines()
    key=lines[0].split()
    if key[:2]!=['#SEP_MB_REFERENCE_KEY','v1.13-csr-v1'] or [int(x) for x in key[2:5]]!=[norb,nel,steps] or float(key[5])!=dt:
        raise ValueError(f'{path}: physical dimensions/time key mismatch')
    footer=lines[-1].split()
    if len(footer)!=2 or footer[0]!='#SEP_MB_REFERENCE_COMPLETE':
        raise ValueError(f'{path}: incomplete reference')
    data=np.loadtxt(path);validate_grid(data,norb,steps,dt)
    checksum=14695981039346656037
    for byte in np.asarray(data[:,1:],dtype='<f8').tobytes():
        checksum=((checksum^byte)*1099511628211)&((1<<64)-1)
    if checksum!=int(footer[1]):raise ValueError(f'{path}: reference checksum mismatch')
    norm=data[:,1]/nel
    if np.any(norm<=0):raise ValueError(f'{path}: nonpositive reference norm')
    result=data.copy();result[:,1:]/=norm[:,None]
    return result,{'reference_states':int(key[7]),'fock_levels':int(key[8]),
                   'raw_norm_max_error':float(abs(norm-1).max())}


def validate_grid(data,norb,steps,dt):
    if data.shape!=(steps+1,norb+4) or not np.isfinite(data).all():
        raise ValueError('incomplete shape or nonfinite values')
    if not np.array_equal(data[:,0],np.arange(steps+1)*dt):
        raise ValueError('time grid is not complete and exact')


def rms(data):return float(np.sqrt(np.mean(np.square(data))))
def ratio(signal,error):
    denominator=rms(error)
    return rms(signal)/denominator if denominator>0 else None


def windows(t,signal,error,active):
    result=[]
    for start in range(0,4000,100):
        ix=(t>=start)&(t<=start+100)
        scores=[ratio(signal[ix,j],error[ix,j]) for j in range(signal.shape[1])]
        finite=[scores[j] for j in active if scores[j] is not None]
        result.append({'start':start,'stop':start+100,'Q_active':ratio(signal[ix][:,active],error[ix][:,active]),
                       'Q_weakest_active':min(finite) if finite else None,'Q_per_orbital':scores,
                       'max_absolute_error_or_SEM':float(abs(error[ix]).max())})
    return result


def physical_check(folder,entry):
    config=dict(line.split('=',1) for line in (folder/'config.txt').read_text().splitlines() if '=' in line)
    expected={'AHM_WC_EV':2.7,'AHM_ETA':.30e-6,'AHM_DELE_EV':-17.194874968839816,
              'AHM_NSTEP':8000,'AHM_SEED':int(entry['seed']),
              'SEP_MB_REFERENCE_DISTANCE':int(entry['distance']),
              'SEP_MB_REFERENCE_FOCK_STATES':int(entry['fock_levels']),'SLURM_NTASKS':int(entry['ranks'])}
    for key,value in expected.items():
        if key not in config or float(config[key])!=value:raise ValueError(f'{folder}: wrong {key}')
    log=(folder/'program.out').read_text()
    for marker in ['nstep=8000 dt=0.5','#SEP_MB_TIMING sampling_seconds=','#SEP_MB_RESOURCE max_rank_peak_rss_kb=']:
        if marker not in log:raise ValueError(f'{folder}: no completion/resource marker {marker}')
    if 'wall_seconds=' not in (folder/'time.txt').read_text():raise ValueError(f'{folder}: missing wall time')
    actual=re.search(r'#SEP_MB_REFERENCE requested_distance=\d+ selected_distance=(\d+) active=1 states=(\d+) max_states=\d+ fock_states=(\d+)',log)
    if not actual or int(actual[1])!=int(entry['distance']) or int(actual[3])!=int(entry['fock_levels']):
        raise ValueError(f'{folder}: actual reference differs from requested manifest')
    if int(actual[2])!={2:52656,3:715136}[int(entry['distance'])]:
        raise ValueError(f'{folder}: unexpected reference state count')
    return log


def plot_orbitals(t,curves,destination,title,initial):
    fig,axes=plt.subplots(6,5,figsize=(16,14),sharex=True)
    for orbital,ax in enumerate(axes.flat):
        for label,data,color,style in curves:
            ax.plot(t,data[:,orbital]-initial[orbital],color=color,ls=style,lw=.65,label=label)
        ax.set_title(f'Orbital {orbital}',fontsize=9);ax.grid(alpha=.15)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
    axes[0,0].legend(fontsize=6)
    for ax in axes[-1]:ax.set_xlabel('Time (a.u.)')
    fig.suptitle(title+'\nAll 30 occupation changes; original output (residual has linear interpolation)',y=.997)
    fig.tight_layout(rect=(0,0,1,.966));fig.savefig(destination,dpi=135);plt.close(fig)


def analyze(directory):
    figure_dir=directory/'figures';figure_dir.mkdir(exist_ok=True)
    entries=list(csv.DictReader((directory/'long_4000_jobs.tsv').open(),delimiter='\t'))
    ready={};states=[];t=np.arange(8001)*.5;initial=np.r_[np.ones(15),np.zeros(15)];active=[0,*range(15,30)]
    for entry in entries:
        folder=directory/entry['job'];pattern=f"ahm-sepmb-s30-n15-{entry['paths']}.dat"
        required=[folder/'program.out',folder/'time.txt',folder/'config.txt',folder/pattern]
        if not all(p.exists() for p in required):
            states.append({'job':entry['job'],'role':entry['role'],'local_status':'awaiting complete downloaded files'});continue
        log=physical_check(folder,entry)
        data=np.loadtxt(folder/pattern);validate_grid(data,30,8000,.5)
        record={'entry':entry,'data':data}
        if entry['role'].startswith('reference'):
            record['reference'],record['reference_info']=load_reference(folder/'reference-observables.dat',30,15,8000,.5)
        ready[entry['job']]=record
        states.append({'job':entry['job'],'role':entry['role'],'local_status':'complete and validated',
                       'wall_seconds':float(re.search(r'wall_seconds=([\d.]+)',(folder/'time.txt').read_text())[1]),
                       'reference_seconds':float(re.search(r'reference_seconds=([\d.]+)',log)[1]),
                       'sampling_seconds':float(re.search(r'sampling_seconds=([\d.]+)',log)[1]),
                       'ranks':int(entry['ranks']),
                       'max_particle_error':float(abs(data[:,4:].sum(axis=1)-15).max()),
                       'minimum_occupation':float(data[:,4:].min()),'maximum_occupation':float(data[:,4:].max()),
                       'max_rank_peak_mib':int(re.search(r'max_rank_peak_rss_kb=(\d+)',log)[1])/1024,
                       'sum_rank_peaks_gib':int(re.search(r'sum_rank_peak_rss_kb=(\d+)',log)[1])/1024**2})
    report={'status':'all seven runs downloaded and validated' if len(ready)==len(entries) else 'awaiting long-run outputs',
            'time_range_au':[0,4000],'orbitals':30,'electrons':15,'active_orbitals':active,'jobs':states,
            'Q_against_full_30_15_QM':None,
            'notice':'Reference agreement and independent-seed repeatability cannot measure a shared systematic bias. No full-space 30/15 QM truth is available.'}
    references=[(rec['entry']['role'],rec) for rec in ready.values() if 'reference' in rec]
    if references:
        curves=[(f"{role}: D{rec['entry']['distance']}/F{rec['entry']['fock_levels']}",rec['reference'][:,4:],color,style)
                for (role,rec),(color,style) in zip(references,[('#278875','-'),('#496fa8','--'),('#ba7735',':')])]
        plot_orbitals(t,curves,figure_dir/'large_4000_references_all_orbitals.png',
                      '30 orbitals / 15 electrons | BOUNDED REFERENCES ONLY | 0-4000 a.u.',initial)
    reference_by_role=dict(references)
    if references:
        selected=reference_by_role.get('reference384',references[0][1])
        values=selected['reference'][:,4:];change=values-initial
        fig=plt.figure(figsize=(11,8.5))
        grid=fig.add_gridspec(3,2,width_ratios=[1,.022])
        axes=[fig.add_subplot(grid[0,0])]
        axes += [fig.add_subplot(grid[row,0],sharex=axes[0]) for row in [1,2]]
        color_axes=[fig.add_subplot(grid[row,1]) for row in [1,2]]
        axes[0].plot(t,change[:,0],color='#283c46',lw=.85)
        axes[0].set_ylabel('Impurity: n0 - 1')
        for ax,orbitals,title in [(axes[1],list(range(1,15)),'Initially occupied bath: n_i - 1'),
                                   (axes[2],list(range(15,30)),'Initially empty bath: n_i')]:
            cmap=plt.get_cmap('viridis');normalizer=plt.Normalize(min(orbitals),max(orbitals))
            for orbital in orbitals:ax.plot(t,change[:,orbital],color=cmap(normalizer(orbital)),lw=.7)
            fig.colorbar(plt.cm.ScalarMappable(norm=normalizer,cmap=cmap),cax=color_axes.pop(0),label='Orbital index',ticks=[orbitals[0],orbitals[len(orbitals)//3],orbitals[2*len(orbitals)//3],orbitals[-1]])
            ax.set_ylabel(title)
        for ax in axes:
            ax.grid(alpha=.18);ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
        axes[-1].set_xlabel('Time (a.u.)')
        entry=selected['entry']
        fig.suptitle(f"30 orbitals / 15 electrons | Complete 0-4000 a.u.\nD{entry['distance']}/F{entry['fock_levels']} bounded reference ONLY; production accuracy not yet validated",fontsize=12)
        fig.tight_layout(rect=(0,0,1,.93));fig.savefig(figure_dir/'large_4000_reference_overview.png',dpi=145);plt.close(fig)
    for label,role in [('oscillator_agreement','reference512'),('reference_order_agreement','reference_d2')]:
        if 'reference384' in reference_by_role and role in reference_by_role:
            ref=reference_by_role['reference384']['reference'][:,4:];other=reference_by_role[role]['reference'][:,4:]
            report[label]={'definition':'RMS(reference384 - initial) / RMS(reference384 - comparison); not absolute QM fitting',
                           'max_orbital_difference':float(abs(ref-other).max()),'windows':windows(t,ref-initial,ref-other,active)}
    samples=[rec for rec in ready.values() if rec['entry']['role']=='sample384']
    if samples and 'reference384' in reference_by_role:
        reference=reference_by_role['reference384']['reference'][:,4:]
        single=[]
        for rec in samples:
            delta=rec['data'][:,4:]-reference
            single.append({'job':rec['entry']['job'],'seed':rec['entry']['seed'],
                'max_orbital_difference_to_reference':float(abs(delta).max()),
                'windows':windows(t,reference-initial,delta,active),
                'notice':'Signal/reference-difference ratio is a diagnostic, not a true QM fitting score. Physical residual and numerical noise are not separated by this statistic.'})
        report['single_run_reference_diagnostics']=single
        rec=samples[0];delta=rec['data'][:,4:]-reference;ww=single[0]['windows']
        curves=[('D3/F384 reference',reference,'#314751','-'),
                (f"1 million paths, seed {rec['entry']['seed']}",rec['data'][:,4:],'#bc7448','-')]
        plot_orbitals(t,curves,figure_dir/'large_4000_first_million_all_orbitals.png',
            '30 orbitals / 15 electrons | First million-path run vs bounded reference',initial)
        fig,axes=plt.subplots(2,2,figsize=(12,8))
        axes[0,0].plot(t,reference[:,0]-1,color='#314751',lw=.8,label='D3/F384 reference')
        axes[0,0].plot(t,rec['data'][:,4]-1,color='#bc7448',lw=.55,label='Million-path result')
        axes[0,0].set(title='Impurity occupation change',xlabel='Time (a.u.)',ylabel='n0 - 1');axes[0,0].legend(fontsize=8)
        axes[0,1].plot(t,delta[:,0],color='#bc7448',lw=.5)
        axes[0,1].set(title='Poisson minus bounded reference (impurity)',xlabel='Time (a.u.)',ylabel='Occupation difference')
        centers=np.arange(50,4000,100)
        axes[1,0].semilogy(centers,[w['Q_active'] for w in ww],label='All active orbitals')
        axes[1,0].semilogy(centers,[w['Q_weakest_active'] for w in ww],label='Weakest active orbital')
        axes[1,0].axhline(10,color='#b44c49',ls='--',lw=1)
        axes[1,0].set(title='Signal / reference-difference diagnostic',xlabel='100-au window midpoint',ylabel='Ratio (NOT Q against full QM)');axes[1,0].legend(fontsize=8)
        d=[r for r in states if r['job'] in ['647702','647703',rec['entry']['job']]]
        bars=axes[1,1].bar(np.arange(len(d)),[r['wall_seconds']/3600 for r in d],color=['#3e8b80','#8c9aaa','#bc7448'])
        for bar,r in zip(bars,d):axes[1,1].text(bar.get_x()+bar.get_width()/2,bar.get_height(),f"{bar.get_height():.2f} h\\n{r['max_rank_peak_mib']:.1f} MiB/rank".replace('\\n','\n'),ha='center',va='bottom',fontsize=8)
        axes[1,1].set(xticks=np.arange(len(d)),xticklabels=['Reference F384','Reference F512','Sampling 1M'][:len(d)],ylabel='Wall hours, 128 ranks',title='Measured stages; sampling reuses F384 reference',ylim=(0,max(r['wall_seconds']/3600 for r in d)*1.23))
        for ax in axes.flat:ax.grid(alpha=.15)
        fig.suptitle('30 orbitals / 15 electrons | Complete 0–4000 au | First million-path result\nLate stochastic fluctuations remain; bounded-reference agreement does not prove full-space accuracy',fontsize=11)
        fig.tight_layout(rect=(0,0,1,.92));fig.savefig(figure_dir/'large_4000_first_million_diagnostics.png',dpi=150);plt.close(fig)
    if len(samples)==3:
        seeds=[rec['entry']['seed'] for rec in samples]
        if len(set(seeds))!=3:raise ValueError('repeated seeds cannot establish independent repeatability')
        values=np.stack([rec['data'][:,4:] for rec in samples]);mean=values.mean(axis=0);sem=values.std(axis=0,ddof=1)/np.sqrt(3)
        report['independent_seed_repeatability']={'seeds':seeds,'replicates':3,'definition':'RMS(mean - initial) / RMS(sample SD/sqrt(3)) within each window',
            'warning':'Three seeds give a noisy SEM estimate; this is not a QM fitting score or a guaranteed confidence bound.',
            'windows':windows(t,mean-initial,sem,active),'max_SEM':float(sem.max()),
            'max_mean_particle_error':float(abs(mean.sum(axis=1)-15).max()),
            'minimum_raw_replicate_occupation':float(values.min())}
        curves=[(f"Poisson seed {rec['entry']['seed']}",rec['data'][:,4:],color,'-')
                for rec,color in zip(samples,['#ae774b','#609882','#7c89af'])]
        curves.append(('3-seed mean',mean,'#28353d','-'))
        plot_orbitals(t,curves,figure_dir/'large_4000_poisson_all_orbitals.png',
                      '30 orbitals / 15 electrons | 3 independent Poisson runs | 0-4000 a.u.',initial)
        np.savez_compressed(directory/'large_4000_replicate_statistics.npz',time=t,mean=mean,standard_error=sem,seeds=seeds)
        fig,ax=plt.subplots(figsize=(10,4.4));qq=report['independent_seed_repeatability']['windows']
        ax.semilogy(np.arange(50,4000,100),[w['Q_active'] for w in qq],label='All active orbitals')
        ax.semilogy(np.arange(50,4000,100),[w['Q_weakest_active'] for w in qq],label='Weakest active orbital')
        ax.axhline(10,color='#b44c49',ls='--',lw=1,label='Diagnostic Q = 10')
        ax.set(xlabel='100-a.u. window midpoint',ylabel='Repeatability Q (SEM of 3 seeds)',
               title='0-4000 a.u. independent-seed noise assessment\nNot a fitting score against full-space QM')
        ax.grid(alpha=.18);ax.legend();fig.tight_layout();fig.savefig(figure_dir/'large_4000_repeatability.png',dpi=150);plt.close(fig)
    if '647704' in ready and '647707' in ready:
        a=ready['647704']['data'][:,4:];b=ready['647707']['data'][:,4:]
        report['matched_seed_fock_difference']={'seed':1134001,'fock_levels':[384,512],
            'max_orbital_difference':float(abs(a-b).max()),'windows':windows(t,a-initial,a-b,active)}
    destination=directory/'large_4000_validation.json';destination.write_text(json.dumps(report,indent=2,allow_nan=False))
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--directory',type=Path,default=DEFAULT);args=parser.parse_args()
    report=analyze(args.directory)
    print(json.dumps({'status':report['status'],'jobs':report['jobs'],'Q_against_full_QM':report['Q_against_full_30_15_QM']},indent=2))
