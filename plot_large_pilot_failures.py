"""Compact view of failed large-space pilots, keeping the full unmodified series."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'scan_v115_generality_20260906'

def main():
    report=json.loads((OUT/'long_4000_validation.json').read_text())
    jobs=[r for r in report['runs'] if r['job'] in ['647780','647781','647782']]
    assert len(jobs)==3
    fig,axes=plt.subplots(2,2,figsize=(12,7.8))
    labels=['50/25 D1','50/17 D1','50/25 no reference'];colors=['#4e8b82','#a16e47','#98505b']
    for ax,row,label,color in zip(axes.flat,jobs,labels,colors):
        nel=17 if row['job']=='647781' else 25
        data=np.loadtxt(OUT/row['job']/f'ahm-sepmb-s50-n{nel}-1000.dat')
        violation=np.maximum(np.max(data[:,4:]-1,axis=1),np.max(-data[:,4:],axis=1))
        ax.semilogy(data[:,0],np.maximum(violation,1e-16),lw=.7,color=color)
        ax.set(title=f'{label}: deviation outside [0, 1]',xlabel='Time (a.u.)',ylabel='Maximum occupation-bound violation')
        ax.grid(alpha=.15)
    ax=axes[1,1];x=np.arange(3);bars=ax.bar(x,[r['resources']['max_rank_peak_mib'] for r in jobs],color=colors)
    for bar,row in zip(bars,jobs):
        ax.text(bar.get_x()+bar.get_width()/2,bar.get_height(),f"{bar.get_height():.2f} MiB\n{row['resources']['wall_seconds']/60:.2f} min",ha='center',va='bottom',fontsize=8)
    ax.set(xticks=x,xticklabels=labels,ylabel='Maximum process peak RSS (MiB)',title='Measured resources; 64 ranks per job',ylim=(0,37));ax.grid(axis='y',alpha=.15)
    fig.suptitle('50-orbital tests | Full 0–4000 au | 1,000-path pilots FAILED accuracy checks\nEvent probabilities were rounded by the fixed-step midpoint grid; low memory alone is insufficient',fontsize=11)
    fig.text(.5,.012,'Full series retained, no clipping or extra smoothing. Log plots floor zero violation at 1e-16. Bound violation is not error against QM.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.025,1,.92));fig.savefig(OUT/'figures/large_pilot_failures.png',dpi=145);plt.close(fig)

if __name__=='__main__':main()
