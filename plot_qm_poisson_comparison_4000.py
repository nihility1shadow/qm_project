"""Compare actual small-system cost with failed large-space precision, without mixing Q definitions."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'scan_versions_4000_20260907'

def main():
    old=json.loads((OUT/'validation.json').read_text())['jobs']
    small=json.loads((ROOT/'scan_v114_signed_csr_20260906/small_4000_validation.json').read_text())['jobs']
    qm=json.loads((ROOT/'scan_v115_generality_20260906/qm_resource_validation.json').read_text())
    large=json.loads((ROOT/'scan_v115_generality_20260906/long_4000_validation.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(12.6,8.4))
    ax=axes[0,0];x=np.arange(len(old));q=[r['Q_against_QM_weakest_active_orbital_window'] for r in old]
    ax.bar(x,q,color='#438d80');ax.set_yscale('log');ax.axhline(10,ls='--',lw=1,color='#ad5349')
    ax.set(xticks=x,xticklabels=[r['version'] for r in old],ylabel='Worst active-orbital/window Q',title='10/5: same 10,000 paths, 64 ranks')
    ax=axes[0,1];labels=['QM\n1 rank','v1.14 1M\n64 ranks'];means=[qm['process_wall_seconds']/60,np.mean([r['resources']['wall_seconds']/60 for r in small])]
    bars=ax.bar(np.arange(2),means,color=['#475866','#438d80'])
    lo=min(r['resources']['wall_seconds']/60 for r in small);hi=max(r['resources']['wall_seconds']/60 for r in small)
    ax.errorbar(1,means[1],yerr=[[means[1]-lo],[hi-means[1]]],fmt='none',color='#222',capsize=5)
    for bar,val in zip(bars,means):ax.text(bar.get_x()+bar.get_width()/2,val+2,f'{val:.2f} min',ha='center',fontsize=9)
    ax.set(xticks=np.arange(2),xticklabels=labels,ylabel='Wall minutes',title='10/5: measured high-accuracy run costs',ylim=(0,76))
    ax.text(.5,.04,'Poisson: Q_min=11.17–16.82 over three runs\nCache consumers; reference preparation adds 24.08 s',transform=ax.transAxes,ha='center',fontsize=8,bbox={'facecolor':'white','edgecolor':'none','alpha':.95})
    ax=axes[1,0];peak=[qm['process_peak_mib'],np.mean([r['resources']['max_rank_peak_mib'] for r in small])]
    bars=ax.bar(np.arange(2),peak,color=['#475866','#438d80'])
    for bar,val in zip(bars,peak):ax.text(bar.get_x()+bar.get_width()/2,val+.5,f'{val:.2f}',ha='center',fontsize=9)
    ax.set(xticks=np.arange(2),xticklabels=labels,ylabel='Maximum process peak RSS (MiB)',title='10/5: process peak is not whole-job memory',ylim=(0,28))
    ax.text(.5,.79,'64-rank Poisson sum of peaks: 0.858–0.866 GiB\nPeak sum is NOT simultaneous total RSS',transform=ax.transAxes,ha='center',fontsize=8)
    ax=axes[1,1];d={r['job']:r for r in large['runs']}
    for job,label,color in [('647780','D1 + 1,000 paths','#a3754e'),('647777','D2 + 1,000 paths','#438d80'),('647778','D2 + 10,000 paths','#6d81ad')]:
        ww=d[job]['reference_difference_windows']
        ax.semilogy(np.arange(50,4000,100),[w['max_absolute_error_or_SEM'] or np.nan for w in ww],label=label,color=color)
    ax.set(title='50/25: large residual fluctuations remain',xlabel='100-au window midpoint',ylabel='Maximum difference to own bounded reference')
    ax.legend(fontsize=8);ax.set_xlim(0,4000)
    for ax in axes.flat:ax.grid(axis='y',alpha=.15)
    fig.suptitle('QM and Poisson comparison | Complete 0–4000 au\nSmall-system Q is against grid QM; large-system reference differences are NOT full-QM errors',fontsize=12)
    fig.text(.5,.01,'Small-system Poisson bars: mean of 3 runs; time whisker: full range. Zero differences omitted on the large-system log plot.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.025,1,.92));fig.savefig(OUT/'figures/qm_poisson_comparison_4000.png',dpi=150);plt.close(fig)

if __name__=='__main__':main()
