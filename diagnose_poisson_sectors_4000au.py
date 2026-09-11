"""Replace only the sampled (1,1) sector by its exact value; retain every higher sector.

The separate second-order curve omits higher sectors and is explicitly diagnostic.
No QM result enters either construction; QM is loaded only for scoring afterwards.
"""
from pathlib import Path
import json,time,sys,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'experiments'))
from poisson_one_jump_sector import one_jump_sector,base
np=base.np
DEST=ROOT/'scan_slater_stream_4000au_20260910'
ACTIVE=[0,5,6,7,8,9];INITIAL=np.r_[np.ones(5),np.zeros(5)]

def scores(t,qm,value):
    rows=[]
    for lo in range(0,4000,100):
        mask=(t>=lo)&(t<=lo+100)
        signal=np.sqrt(np.mean((qm[mask]-INITIAL)**2,axis=0))
        error=np.sqrt(np.mean((value[mask]-qm[mask])**2,axis=0));q=signal/np.maximum(error,1e-300)
        rows.append(dict(start=lo,stop=lo+100,q_per_orbital=q.tolist(),worst_active_q=float(q[ACTIVE].min())))
    return dict(windows=rows,worst_active_window_q=min(s['worst_active_q'] for s in rows),
        max_absolute_error=float(np.max(abs(value-qm))))

def main():
    built=[]
    # Estimator construction is complete before loading the QM validation data.
    for parent_name in ['pilot_n256_step10','poisson_n4096_seed1160101']:
        parent=DEST/parent_name
        if not (parent/'metrics.json').exists():continue
        started=time.perf_counter();cpu=time.process_time()
        m=json.loads((parent/'metrics.json').read_text());raw=np.load(parent/'occupations.npz')
        sectors=np.load(parent/'sector_contributions.npz')['contributions']
        t=raw['time'];original=raw['occupations'];exact,details=one_jump_sector(m['model'],t)
        assert np.max(abs(original-(INITIAL+sectors.sum(axis=0))))<1e-12
        rb=INITIAL+exact+sectors[1:].sum(axis=0)
        out=DEST/(parent_name+'_rb11');out.mkdir(exist_ok=True)
        np.savez_compressed(out/'occupations.npz',time=t,occupations=rb)
        np.savetxt(out/'occupations.dat',np.c_[t,rb],header='Full Poisson: exact (1,1) sector + all sampled higher sectors; time_au orbital_0..9')
        mi=base.psutil.Process().memory_info();peak=getattr(mi,'peak_wset',mi.rss)
        info=dict(method='Full Poisson with exact (1,1) sector replacement',parent_run=parent_name,
            parent_solver_seconds=m['wall_seconds'],postprocess_seconds=time.perf_counter()-started,
            parent_cpu_seconds=m['cpu_seconds'],postprocess_cpu_seconds=time.process_time()-cpu,
            parent_peak_rss_bytes=m['peak_rss_bytes'],postprocess_peak_rss_bytes=peak,
            separate_process_pipeline_peak_bytes=max(m['peak_rss_bytes'],peak),analytic_sector=details,
            full_expansion_retained=True,qm_input=False,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            analytic_source_sha256=hashlib.sha256((ROOT/'experiments/poisson_one_jump_sector.py').read_bytes()).hexdigest(),
            sampling_reuse='Reuses the parent samples for all higher sectors; not an independent random run')
        (out/'estimator_metadata.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
        built.append((parent_name,t,original,rb,INITIAL+exact,out,info))
    reference_name='qm_f768' if (DEST/'qm_f768/metrics.json').exists() else 'qm_f512'
    qm=np.load(DEST/reference_name/'occupations.npz')
    allrows=[]
    for name,t,raw,rb,truncated,out,info in built:
        idx=np.searchsorted(qm['time'],t);assert np.array_equal(qm['time'][idx],t)
        reference=qm['occupations'][idx]
        row=dict(parent=name,qm_reference=reference_name,full_raw=scores(t,reference,raw),full_rb11=scores(t,reference,rb),
            second_order_only=scores(t,reference,truncated),resources=info,
            interpretation='The second-order-only curve is biased/truncated and cannot be reported as full Poisson accuracy.')
        (out/'comparison.json').write_text(json.dumps(row,indent=2),encoding='utf-8');allrows.append(row)
    (DEST/'sector_diagnostics.json').write_text(json.dumps(allrows,indent=2),encoding='utf-8')
    if not built:raise ValueError('No complete parent')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font_manager.fontManager.addfont('C:/Windows/Fonts/msyh.ttc')
    plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False})
    name,t,raw,rb,truncated,out,info=built[-1];row=allrows[-1];idx=np.searchsorted(qm['time'],t);ref=qm['occupations'][idx]
    fig,axes=plt.subplots(1,2,figsize=(13,4.7))
    versions=[('完整泊松：原采样',raw,row['full_raw'],'#d55e00'),('完整泊松：解析替换(1,1)层',rb,row['full_rb11'],'#0072b2'),('仅二阶截断：省略高阶项',truncated,row['second_order_only'],'#009e73')]
    for label,value,sc,color in versions:
        axes[0].plot(t,np.maximum(np.max(abs(value-ref),axis=1),1e-18),lw=1,color=color,label=label)
        axes[1].plot([s['start']+50 for s in sc['windows']],[s['worst_active_q'] for s in sc['windows']],lw=1,color=color,label=label)
    axes[0].set(yscale='log',xlabel='时间 / a.u.',ylabel='全轨道最大绝对误差',title='最低阶解析积分可减少部分噪声')
    axes[1].set(yscale='log',xlabel='窗口中心 / a.u.',ylabel='最差活跃轨道 Q',title='绿色截断曲线不能当作完整泊松达标');axes[1].axhline(10,color='#555',ls='--')
    for ax in axes:ax.grid(alpha=.18)
    axes[0].legend(fontsize=9,frameon=False);fig.suptitle(f'泊松分层诊断｜{name}',fontsize=13);fig.tight_layout(rect=(0,0,1,.93))
    for dpi,suffix in [(160,''),(85,'_preview')]:fig.savefig(DEST/f'sector_diagnostics{suffix}.png',dpi=dpi,facecolor='white')
    plt.close(fig)
    fig,axes=plt.subplots(5,2,figsize=(14,15),sharex=True)
    for j,ax in enumerate(axes.flat):
        ax.plot(t,rb[:,j]-INITIAL[j],color='#0072b2',lw=.8,label='完整泊松：解析(1,1)＋全部高阶采样')
        ax.plot(t,ref[:,j]-INITIAL[j],color='#20252b',lw=1.1,label='完整 QM',zorder=4)
        worst=min(w['q_per_orbital'][j] for w in row['full_rb11']['windows'])
        ax.set_title(f'轨道 {j}｜最差100 a.u.窗口 Q={worst:.3g}',fontsize=10)
        ax.set_ylabel('占据数变化 Δn');ax.set_xlim(0,4000);ax.grid(alpha=.18)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
        if j not in ACTIVE:ax.set_facecolor('#fff4ed')
    for ax in axes[-1]:ax.set_xlabel('时间 / a.u.')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.945),ncol=2,frameon=False,fontsize=9)
    fig.suptitle(f'5电子 / 10轨道｜完整泊松解析分层降噪｜0–4000 a.u.\nη=3e-7 Ha²，wc=2.7 eV；父运行 {name}\n最差活跃窗口 Q={row["full_rb11"]["worst_active_window_q"]:.3g}；全部高阶项仍保留',fontsize=12,y=.995)
    fig.text(.5,.012,'原始单次高阶样本；未裁剪、平滑或平均种子。此图与原泊松共享高阶样本，不是独立重复。',ha='center',fontsize=9)
    fig.tight_layout(rect=(.02,.025,.99,.925),h_pad=1.5)
    for dpi,suffix in [(160,''),(80,'_preview')]:fig.savefig(DEST/f'rb11_all_orbitals_4000au{suffix}.png',dpi=dpi,facecolor='white')
    plt.close(fig)
    lines=['# 泊松层解析积分诊断','',
        '这里严格区分三个估计器：完整随机泊松；只把(1,1)层替换为其解析值、仍保留全部高阶采样的完整泊松；仅保留(1,1)层的二阶截断近似。三者均不读取QM进行构造，QM仅用于事后评分。','',
        '对当前初态和正能隙，最低非零层可以写为：','',
        'Δn_j^(1,1)(t) = 4 c_j² Σν Pν sin²[(ε_j−ε_d+νω)t/2] / (ε_j−ε_d+νω)²，j为初始空轨道；Δn_0 = −Σj Δn_j。', '',
        'Pν是均值d²的泊松分布，即初始位移谐振子基态在未位移振子数态中的概率。这里积分掉的是一个已明确分层的泊松项，不引入参考波函数。核级数尾部有明确的占据误差上界；这个尾界不约束被二阶截断省略的高阶电子跳跃。','',
        '与原始Slater路径核的双重Gauss积分在1、10、50、100 a.u.交叉验证一致，具体误差见analytic_11_diagnostic/diagnostic.json。','',
        '| 父运行 | 完整原采样最差窗口Q | 完整解析替换(1,1)最差窗口Q | 仅二阶截断最差活跃窗口Q |','| --- | ---: | ---: | ---: |']
    for r in allrows:lines.append(f'| {r["parent"]} | {r["full_raw"]["worst_active_window_q"]:.6g} | {r["full_rb11"]["worst_active_window_q"]:.6g} | {r["second_order_only"]["worst_active_window_q"]:.6g} |')
    lines+=['','仅二阶截断对于初始占据的浴轨道1–4给出零变化，因此这些轨道的Q为1；即使活跃轨道拟合很好，也不是全部轨道达标，更不是完整泊松算法已解决4000 a.u.问题。','',
        '解析替换复用原运行高阶样本，因此它与原曲线相关，不是新的独立随机验证。完整方法仍受高阶项相位相消和抽样噪声限制。','',
        '[完整解析降噪泊松全部轨道](rb11_all_orbitals_4000au.png) · [诊断图](sector_diagnostics.png) · [全部数值对照](sector_diagnostics.json)','',
        '复现：先运行 experiments/poisson_one_jump_sector.py，再运行 diagnose_poisson_sectors_4000au.py。']
    (DEST/'SECTOR_DIAGNOSTICS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps([{r['parent']:{k:r[k]['worst_active_window_q'] for k in ['full_raw','full_rb11','second_order_only']}} for r in allrows],indent=2))
if __name__=='__main__':main()
