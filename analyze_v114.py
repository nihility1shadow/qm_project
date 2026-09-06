"""Measured signed-CSR regressions. Preserve all timing results, including slow runs."""
import json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;OLD=ROOT/'scan_v113_two_buffer_20260905';OUT=ROOT/'scan_v114_signed_csr_20260906'
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)

def load(folder):
    matches=list(folder.glob('ahm-sepmb-*.dat'));assert len(matches)==1
    values=np.loadtxt(matches[0]);assert np.isfinite(values).all();return values

def resource(folder):
    log=(folder/'program.out').read_text();timing=(folder/'time.txt').read_text()
    def find(pattern):
        match=re.search(pattern,log);return float(match[1]) if match else None
    return {'wall_seconds':float(re.search(r'wall_seconds=([\d.]+)',timing)[1]),
            'reference_seconds':find(r'reference_seconds=([\d.]+)'),
            'max_rank_peak_mib':find(r'max_rank_peak_rss_kb=(\d+)')/1024,
            'sum_rank_peaks_gib':find(r'sum_rank_peak_rss_kb=(\d+)')/1024**2,
            'graph_metadata_mib':find(r'metadata_mib=([\d.]+)')}

pairs=[]
for a,b in [(647652,647709),(647652,647710),(647653,647711),(647654,647712),(647655,647713),(647662,647714)]:
    old,new=load(OLD/str(a)),load(OUT/str(b));assert old.shape==new.shape
    error=float(abs(old-new).max());assert error==0,(a,b,error)
    pairs.append({'old_job':a,'new_job':b,'max_abs_all_columns':error,'old':resource(OLD/str(a)),'new':resource(OUT/str(b))})
report={'regressions':pairs,'physical_parameters':'unchanged','sampling_rng_and_sum_order':'unchanged',
        'status':'Six cloud regressions passed; same-allocation AB timing checked separately',
        'timing_note':'The unpinned 128-rank run is slower than the earlier v1.13 run; it is retained in the report. Use the explicit 64-ranks/node same-allocation AB result to assess timing.',
        'peak_sum_note':'Sum of individual process peaks is not simultaneous aggregate RSS.'}
abroot=OUT/'647716'
if all((abroot/label/'time.txt').exists() for label in ['old_first','new_first','new_second','old_second']):
    ab={label:resource(abroot/label) for label in ['old_first','new_first','new_second','old_second']}
    expected=load(abroot/'old_first')
    for label in ab:assert np.array_equal(expected,load(abroot/label)),label
    ab['old_mean_reference_seconds']=float(np.mean([ab[k]['reference_seconds'] for k in ['old_first','old_second']]))
    ab['new_mean_reference_seconds']=float(np.mean([ab[k]['reference_seconds'] for k in ['new_first','new_second']]))
    ab['reference_time_ratio_new_over_old']=ab['new_mean_reference_seconds']/ab['old_mean_reference_seconds']
    report['same_allocation_AB']=ab
(OUT/'validation.json').write_text(json.dumps(report,indent=2))
fig,axes=plt.subplots(1,3,figsize=(12,4.6))
for ax,field,title in zip(axes,['max_rank_peak_mib','sum_rank_peaks_gib','graph_metadata_mib'],
                         ['Maximum rank peak (MiB)','Sum of rank peaks (GiB)','Graph table per rank (MiB)']):
    x=np.arange(3);rows=pairs[2:5]
    for offset,version,color in [(-.18,'old','#97a5b1'),(.18,'new','#288775')]:
        bars=ax.bar(x+offset,[row[version][field] for row in rows],width=.36,color=color,label='v1.13' if version=='old' else 'Signed CSR')
        for bar in bars:ax.text(bar.get_x()+bar.get_width()/2,bar.get_height(),f'{bar.get_height():.2f}',ha='center',va='bottom',fontsize=8)
    ax.set_xticks(x,['D2 / 64 ranks','D3 / 64 ranks','D3 / 128 ranks']);ax.tick_params(axis='x',labelsize=8)
    ax.set_title(title,fontsize=10);ax.grid(axis='y',alpha=.15);ax.margins(y=.2)
axes[0].legend(fontsize=8)
fig.suptitle('30 orbitals / 15 electrons | Measured memory reduction, identical numerical results\n10 a.u. / 100 paths / F384; sums of process peaks are not simultaneous RSS',fontsize=11)
fig.tight_layout(rect=(0,0,1,.9));fig.savefig(FIG/'signed_csr_memory.png',dpi=150);plt.close(fig)
print(json.dumps(report,indent=2))
