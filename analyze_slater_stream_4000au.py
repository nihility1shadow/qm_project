"""Audit and plot 4000-au raw Poisson trajectories against converged full QM."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parent
DEST=ROOT/'scan_slater_stream_4000au_20260910'
ACTIVE=[0,5,6,7,8,9]
INITIAL=np.r_[np.ones(5),np.zeros(5)]

def read_run(name):
    p=DEST/name;m=json.loads((p/'metrics.json').read_text(encoding='utf-8'))
    z=np.load(p/'occupations.npz');t=z['time'];v=z['occupations']
    assert len(t)==m['observation_count'] and t[0]==0 and t[-1]==4000 and np.isfinite(v).all()
    for file,key in [('source_snapshot.py','source_sha256'),('kernel_snapshot.py','kernel_sha256')]:
        assert hashlib.sha256((p/file).read_bytes()).hexdigest()==m[key],file
    return t,v,m

def score(ref,value):
    signal=np.sqrt(np.mean((ref-INITIAL)**2,axis=0));err=value-ref
    rmse=np.sqrt(np.mean(err**2,axis=0));q=signal/np.maximum(rmse,1e-300)
    return dict(q_per_orbital=q.tolist(),rmse_per_orbital=rmse.tolist(),signal_rms_per_orbital=signal.tolist(),
        worst_active_q=float(q[ACTIVE].min()),worst_all_orbital_q=float(q.min()),max_absolute_error=float(np.max(abs(err))))

def main():
    for font in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/simhei.ttf']:
        if Path(font).exists():font_manager.fontManager.addfont(font)
    plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    t512,q512,m512=read_run('qm_f512');t384,q384,m384=read_run('qm_f384')
    tq,qm,qmm=read_run('qm_f768') if (DEST/'qm_f768/metrics.json').exists() else (t512,q512,m512)
    primary_nfock=qmm['nfock']
    assert np.array_equal(tq,t384) and np.array_equal(tq,t512) and qmm['model']==m384['model']==m512['model']
    convergence=float(np.max(abs(qm-q512))) if primary_nfock==768 else float(np.max(abs(q512-q384)))
    runs=[]
    for p in sorted(DEST.iterdir()):
        if not p.is_dir() or not (p/'metrics.json').exists():continue
        m=json.loads((p/'metrics.json').read_text(encoding='utf-8'))
        if not m['method'].startswith('stratified-lattice'):continue
        t,v,m=read_run(p.name);idx=np.searchsorted(tq,t)
        assert np.array_equal(tq[idx],t) and m['model']==qmm['model']
        ref=qm[idx];windows=[];prefix=0
        for lo in range(0,4000,100):
            mask=(t>=lo)&(t<=lo+100);s=score(ref[mask],v[mask])
            windows.append(dict(start=lo,stop=lo+100,observations=int(mask.sum()),**s))
            if prefix==lo and s['worst_active_q']>=10:prefix=lo+100
        row=dict(name=p.name,seed=m['seed'],pairs_per_stratum=m['pairs_per_sampled_stratum'],sampled_pairs=m['sampled_pairs'],
            batch=m['batch'],time_block=m['time_block'],output_step=m['output_step'],wall_seconds=m['wall_seconds'],cpu_seconds=m['cpu_seconds'],
            peak_rss_mib=m['peak_rss_bytes']/2**20,full_interval=score(ref,v),windows=windows,
            worst_active_window_q=min(s['worst_active_q'] for s in windows),q10_contiguous_end_au=prefix,
            windows_passed=sum(s['worst_active_q']>=10 for s in windows),max_particle_number_error=m['max_particle_number_error'],
            min_occupation=m['min_occupation'],max_occupation=m['max_occupation'])
        if m['seed']==1160101 and m['pairs_per_sampled_stratum']==4096:
            old=np.load(ROOT/'scan_slater_stream_100au_20260909/final_seed1/occupations.npz')
            mask=t<=100;ii=np.searchsorted(old['time'],t[mask]);assert np.array_equal(old['time'][ii],t[mask])
            row['first_100au_vs_archived_run_max_difference']=float(np.max(abs(v[mask]-old['occupations'][ii])))
            assert row['first_100au_vs_archived_run_max_difference']<1e-12
        runs.append((row,t,v,ref))
    if not runs:raise ValueError('No completed Poisson output')
    primary=next((x for x in runs if x[0]['name']=='poisson_n4096_seed1160101'),runs[0])
    row,t,v,ref=primary
    summary=dict(scope='5 electrons / 10 orbitals; continuous-time propagation to 4000 au; no reference wavefunction in Poisson',
        status='FORMAL_COMPLETE' if row['pairs_per_stratum']==4096 else 'PILOT_ONLY_FORMAL_RUNNING',
        physical_parameters=qmm['model'],active_orbitals=ACTIVE,q_definition='RMS(QM occupation minus initial) / RMS(Poisson minus QM)',
        qm_f384=m384,qm_f512=m512,qm_primary=qmm,qm_primary_nfock=primary_nfock,
        qm_max_f384_f512_difference=float(np.max(abs(q512-q384))),qm_primary_comparison_max_difference=convergence,runs=[r[0] for r in runs],
        implementation_validation=json.loads((DEST/'implementation_validation.json').read_text()),
        limitations=['Q is evaluated on the explicitly recorded observation grid. The pilot grid of 10 au is coarse.',
            'No clipping, smoothing, reference substitution, or averaging independent seeds.',
            'A complete output at 4000 au does not by itself meet the accuracy target.',
            'The earlier interrupted 20260909 launches produced no complete output and are not included as completed benchmarks.',
            'Memory is whole-process peak working set; time includes checkpoint I/O and excludes uncheckpointed interrupted work.'])
    (DEST/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    with (DEST/'window_metrics.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['run','start_au','stop_au','observations','orbital','active','signal_RMS','RMSE','Q'])
        for r,_,_,_ in runs:
            for x in r['windows']:
                for j in range(10):w.writerow([r['name'],x['start'],x['stop'],x['observations'],j,j in ACTIVE,x['signal_rms_per_orbital'][j],x['rmse_per_orbital'][j],x['q_per_orbital'][j]])
    fig,axes=plt.subplots(5,2,figsize=(14,15),sharex=True)
    for j,ax in enumerate(axes.flat):
        ax.plot(t,v[:,j]-INITIAL[j],color='#d55e00',lw=.8,label='泊松：原始单次输出')
        ax.plot(tq,qm[:,j]-INITIAL[j],color='#20252b',lw=1.05,label=f'完整 QM，F{primary_nfock}',zorder=4)
        worst=min(x['q_per_orbital'][j] for x in row['windows'])
        ax.set_title(f'轨道 {j}｜最差 100 a.u. 窗口 Q={worst:.3g}',fontsize=10)
        ax.set_ylabel('占据数变化 Δn');ax.set_xlim(0,4000);ax.grid(alpha=.18)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
        if j not in ACTIVE:ax.set_facecolor('#fff4ed')
    for ax in axes[-1]:ax.set_xlabel('时间 / a.u.')
    h,l=axes[0,0].get_legend_handles_labels();fig.legend(h,l,loc='upper center',bbox_to_anchor=(.5,.947),ncol=2,frameon=False)
    title='正式运行' if row['pairs_per_stratum']==4096 else '小样本预跑（正式运行尚未完成）'
    fig.suptitle(f'5 电子 / 10 轨道｜0–4000 a.u.｜{title}\nη=3e-7 Ha²，wc=2.7 eV，εd=-17.194875 eV；每层 {row["pairs_per_stratum"]} 点，种子 {row["seed"]}\n最差活跃窗口 Q={row["worst_active_window_q"]:.3g}；求解 {row["wall_seconds"]:.1f} 秒；峰值 {row["peak_rss_mib"]:.1f} MiB',fontsize=12,y=.995)
    fig.text(.5,.012,'纵轴独立缩放；全部原始点，未平滑或裁剪。若泊松噪声远大于 QM 信号，QM 曲线会接近横轴。',ha='center',fontsize=9)
    fig.tight_layout(rect=(.02,.025,.99,.925),h_pad=1.5)
    for dpi,suffix in [(160,''),(80,'_preview')]:fig.savefig(DEST/f'all_orbitals_4000au{suffix}.png',dpi=dpi,facecolor='white')
    plt.close(fig)
    fig,axes=plt.subplots(3,2,figsize=(12,8),sharex=True)
    early=t<=100
    for j,ax in zip(ACTIVE,axes.flat):
        ax.plot(t[early],ref[early,j]-INITIAL[j],color='#20252b',lw=1.8,label='QM')
        ax.plot(t[early],v[early,j]-INITIAL[j],color='#d55e00',lw=1,label='泊松')
        ax.set_title(f'轨道 {j}');ax.set_ylabel('Δn');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False);ax.grid(alpha=.18)
    for ax in axes[-1]:ax.set_xlabel('时间 / a.u.')
    fig.suptitle('同一次 4000 a.u. 运行的前 100 a.u.｜未拼接其他运行',fontsize=13)
    axes[0,0].legend(frameon=False);fig.tight_layout(rect=(0,0,1,.95))
    for dpi,suffix in [(160,''),(80,'_preview')]:fig.savefig(DEST/f'first_100au{suffix}.png',dpi=dpi,facecolor='white')
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(15,4.7))
    colors=['#0072b2','#d55e00','#009e73']
    for i,(r,tt,vv,rr) in enumerate(runs):
        label=f'N={r["pairs_per_stratum"]}/层，间隔{r["output_step"]:g}'
        x=[s['start']+50 for s in r['windows']];y=[s['worst_active_q'] for s in r['windows']]
        axes[0].plot(x,y,'o-',ms=3,lw=1,label=label,color=colors[i%3])
        err=np.max(abs(vv-rr),axis=1);axes[1].plot(tt,np.maximum(err,1e-18),lw=.8,label=label,color=colors[i%3])
    axes[0].axhline(10,color='#555',ls='--');axes[0].set(yscale='log',xlabel='窗口中心 / a.u.',ylabel='最差活跃轨道 Q',title='40 个窗口逐一检查');axes[0].legend(fontsize=8,frameon=False)
    axes[1].set(yscale='log',xlabel='时间 / a.u.',ylabel='全轨道最大绝对误差',title='误差随时间变化')
    names=['QM F384',f'QM F{primary_nfock}']+[f'泊松 N={r[0]["pairs_per_stratum"]}' for r in runs]
    mem=[m384['peak_rss_bytes']/2**20,qmm['peak_rss_bytes']/2**20]+[r[0]['peak_rss_mib'] for r in runs]
    axes[2].bar(names,mem,color=['#777','#333']+[colors[i%3] for i in range(len(runs))]);axes[2].set(ylabel='峰值进程内存 / MiB',title='相同模型，不同精度与预算');axes[2].tick_params(axis='x',rotation=25)
    for i,a in enumerate(mem):axes[2].text(i,a+2,f'{a:.1f}',ha='center',fontsize=9)
    for ax in axes:ax.grid(axis='y',alpha=.18)
    fig.suptitle('4000 a.u. 验证｜低内存与拟合精度分别报告',fontsize=13);fig.tight_layout(rect=(0,0,1,.94))
    for dpi,suffix in [(160,''),(80,'_preview')]:fig.savefig(DEST/f'quality_and_memory{suffix}.png',dpi=dpi,facecolor='white')
    plt.close(fig)
    verdict='全程活跃窗口达到 Q≥10' if row['worst_active_window_q']>=10 else '全程活跃窗口未达到 Q≥10，当前参数不满足4000 a.u.高拟合目标'
    lines=['# 无参考泊松 4000 a.u. 阶段验证（2026-09-10）','',verdict+'。','',
        f'状态：{summary["status"]}。正式运行沿用100 a.u.测试的每层4096点、种子1160101和原始物理参数，仅延长时间并采用经验证的时间批处理。QM只作为外部对照。','',
        '| 运行 | 每层采样点 | 观察间隔 / a.u. | 最差活跃窗口 Q | 从零连续达标到 / a.u. | 墙钟秒 | CPU秒 | 峰值MiB |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r,_,_,_ in runs:lines.append(f'| {r["name"]} | {r["pairs_per_stratum"]} | {r["output_step"]} | {r["worst_active_window_q"]:.6g} | {r["q10_contiguous_end_au"]} | {r["wall_seconds"]:.2f} | {r["cpu_seconds"]:.2f} | {r["peak_rss_mib"]:.2f} |')
    lines+=['',f'主要完整泊松运行的全轨道最大绝对误差为 {row["full_interval"]["max_absolute_error"]:.6e}；原始占据值范围为 {row["min_occupation"]:.9g} 到 {row["max_occupation"]:.9g}。', '',f'QM F384/F512的全轨道最大差异为 {np.max(abs(q512-q384)):.6e}；当前F{primary_nfock}与次大基组的差异为 {convergence:.6e}。F{primary_nfock}耗时 {qmm["wall_seconds"]:.2f} 秒、峰值 {qmm["peak_rss_bytes"]/2**20:.2f} MiB；F384耗时 {m384["wall_seconds"]:.2f} 秒、峰值 {m384["peak_rss_bytes"]/2**20:.2f} MiB。', '',
        'Q = RMS(QM占据变化)/RMS(泊松−QM)，活跃轨道固定为0、5–9，仍展示全部10轨道。每个100 a.u.窗口分别评分，不能用全程平均掩盖坏区间。各时刻复用格点，窗口误差相关，不是40次独立随机实验。两次运行的预算、种子与观察网格不同，不构成严格的样本收敛序列。观察间隔不是传播积分步长：QM和泊松均采用连续时间传播。10 a.u.预跑网格较粗，不作为通过高精度目标的证明。','',
        '## 图与原始数据','',
        '[全部轨道](all_orbitals_4000au.png) · [同次运行前100 a.u.](first_100au.png) · [逐窗口精度与内存](quality_and_memory.png) · [逐轨道CSV](window_metrics.csv) · [完整参数与统计](summary.json)','',
        '每次运行目录保存完整命令、物理和数值参数、源代码及原路径核的SHA256与快照、原始占据数、资源记录。所有图均使用原始输出，未平滑、裁剪、拼接或平均独立种子。','',
        '## 与历史版本的关系','',
        '历史六版本4000 a.u.测试使用每版10000条前向路径、64进程及历史网格QM：v0.93/v0.96/v0.98最差窗口Q约为1.93e-5、5.61e-6、4.59e-6；带参考控制的v1.03/v1.09/v1.14约为1.58–1.59。百万路径v1.14曾达到Q约11–17，但对应D3参考单独使用约2781，不能把混合输出达标当作泊松修正有益的证明。', '',
        '本轮是单进程、每个非零层4096个配对积分点及连续时间QM，物理参数相同，样本单位、核表示、硬件与资源口径不同，不能把这些数值作为相同精度下的加速排名。详见[历史比较](../QM_POISSON_COMPARISON_4000.md)。', '',
        '本轮另外实现了[解析泊松层降噪诊断](SECTOR_DIAGNOSTICS.md)。完整解析替换仍保留全部高阶采样；仅二阶截断则省略高阶项，二者不能混称。', '',
        '## 算法、低内存与限制','',
        '保留Slater集体跳跃、1/2/3/≥4跳跃分层、随机平移格点和三次时间变换的完整雅可比。每条路径存Norb×Nel矩阵与核相干态；不分配完整电子空间波函数。批量128，每次合并8个观察时刻，最多同时处理1024对路径，内存不随总采样点数线性增加。','',
        '固定跳数的两个异奇偶层严格为零；其余8层均计算。≥4尾部包含完整条件泊松分布，24个事件缓冲不足时抛出错误，不丢弃高阶事件。平均每腿跳跃数rT从100 a.u.的0.0548增至4000 a.u.的2.1909，长时间相位相消的统计困难仍然存在。','',
        '仅改善批处理速度和断点能力不会自动解决长时间方差。即使内存较低，若未达到同一精度目标，也不能据此宣称优于QM或完成了论文中用时间换内存的目标。大体系内核检查不等于大体系完整轨道拟合验证。','',
        '内存是Windows整个进程峰值工作集，含Python和数值库。墙钟/CPU记录包含求解及断点写入；多个单线程进程曾同时运行。断点可恢复，但被中断时未保存的最后一批工作不计入累计值。20260909的三次中断启动没有完整输出，不冒充已完成测量。','',
        '## 验证与复现','',
        '时间批处理在4/2、10/5、30/15、50/25（轨道/电子）与原实现一致；泊松及QM强制中断再恢复检查通过。数值验证与速度记录见implementation_validation.json。','',
        '```powershell',
        'python tests/test_poisson_slater_long.py',
        'python experiments/poisson_slater_long.py --method qm --nfock 512 --out scan_slater_stream_4000au_20260910/reproduce_qm',
        'python experiments/poisson_slater_long.py --method poisson --pairs 4096 --seed 1160101 --out scan_slater_stream_4000au_20260910/reproduce_poisson',
        'python analyze_slater_stream_4000au.py','```','',
        '确认没有同目录实例仍在运行后，以相同参数和输出目录重新启动会恢复未完成断点；参数或源代码不匹配时拒绝恢复，已完成目录拒绝覆盖。绘图读取归档正式目录，不自动改用reproduce目录。','',
        '本阶段本地回退标签：v1.17-slater-poisson-4000au。前一100 a.u.本地回退点：v1.16-slater-poisson-100au。GitHub上传仍受之前自动审批拒绝的限制，本报告不声称已经上传。']
    if 'first_100au_vs_archived_run_max_difference' in row:lines.insert(4,f'本次前100 a.u.与已归档同种子结果在共同观察点上的最大差异为 {row["first_100au_vs_archived_run_max_difference"]:.6e}；该区间没有替换或拼接历史曲线。\n')
    (DEST/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=summary['status'],qm_convergence=convergence,runs=[{k:r[0][k] for k in ['name','wall_seconds','peak_rss_mib','worst_active_window_q','q10_contiguous_end_au']} for r in runs]),indent=2))
if __name__=='__main__':main()
