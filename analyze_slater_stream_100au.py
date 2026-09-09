"""Render all orbitals and report single-run metrics; no smoothing or seed averaging."""
from pathlib import Path
import json,hashlib,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parent
DEST=ROOT/'scan_slater_stream_100au_20260909'
ACTIVE=[0,5,6,7,8,9]
WEAK=[1,2,3,4]
COLORS=['#d55e00','#0072b2','#009e73']

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read_run(name):
 p=DEST/name
 d=np.load(p/'occupations.npz')
 m=json.loads((p/'metrics.json').read_text(encoding='utf-8'))
 assert np.isfinite(d['occupations']).all()
 if (p/'source_snapshot.py').exists():assert digest(p/'source_snapshot.py')==m['source_sha256']
 return d['time'],d['occupations'],m

def score(ref,value,initial):
 signal=np.sqrt(np.mean((ref-initial)**2,axis=0));error=value-ref
 rmse=np.sqrt(np.mean(error**2,axis=0));q=signal/np.maximum(rmse,1e-300)
 return dict(q_per_orbital=q.tolist(),signal_rms_per_orbital=signal.tolist(),rmse_per_orbital=rmse.tolist(),
             worst_active_q=float(min(q[ACTIVE])),worst_all_orbital_q=float(q.min()),max_absolute_error_all_orbitals=float(np.max(abs(error))),
             active_aggregate_q=float(np.sqrt(np.mean((ref[:,ACTIVE]-initial[ACTIVE])**2))/np.sqrt(np.mean(error[:,ACTIVE]**2))))

def historical_subset(source,name,t):
 p=ROOT/source;raw=np.loadtxt(p);idx=np.searchsorted(raw[:,0],t)
 assert np.allclose(raw[idx,0],t,rtol=0,atol=1e-10)
 selected=raw[idx]
 np.savetxt(DEST/name,selected,fmt='%.17e',header='Exact common observation times; source '+source)
 return selected[:,4:],dict(source_repository_path=source,source_sha256=digest(p),subset=name,subset_sha256=digest(DEST/name))

def main():
 for font in ['C:/Windows/Fonts/msyh.ttc','C:/Windows/Fonts/simhei.ttf']:
  if Path(font).exists():font_manager.fontManager.addfont(font)
 plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
 t,qm,qmm=read_run('final_qm224');_,qm192,qm192m=read_run('final_qm192')
 initial=np.r_[np.ones(5),np.zeros(5)]
 assert t[0]==0 and t[-1]==100 and len(t)==101
 grid,gridsource=historical_subset('scan_v113_two_buffer_20260905/647656/ahm-qm-s10-n5.dat','historical_grid_qm_0_100.dat',t)
 old,oldsource=historical_subset('pure_poisson_figures_20260908/s10n5_t100_n2m_job647976/ahm-sepmb-s10-n5-2000000.dat','historical_v093_0_100.dat',t)
 rows=[];runs=[]
 for i in range(1,4):
  tt,v,m=read_run('final_seed'+str(i));assert np.array_equal(tt,t);assert m['model']==qmm['model']
  row=dict(case='final_seed'+str(i),seed=m['seed'],continuous_qm=score(qm,v,initial),historical_grid_qm=score(grid,v,initial),
           wall_seconds=m['wall_seconds'],cpu_seconds=m['cpu_seconds'],peak_rss_mib=m['peak_rss_bytes']/2**20,
           sampled_pairs=m['sampled_pairs'],pairs_per_stratum=m['pairs_per_sampled_stratum'],batch=m['batch'],
           max_particle_number_error=m['max_particle_number_error'],source_sha256=m['source_sha256'])
  rows.append(row);runs.append(v)
 _,sobol,sobolm=read_run('final_sobol4096');sobolscore=score(qm,sobol,initial)
 assert all(r['continuous_qm']['worst_active_q']>10 for r in rows)
 summary=dict(scope='5 electrons / 10 orbitals, 0-100 au, 101 exact observation times; original fixed bath',
  physical_parameters=qmm['model'],active_orbitals=ACTIVE,near_stationary_orbitals=WEAK,
  q_definition='RMS(QM occupation minus initial occupation) / RMS(estimator minus QM), per orbital, one full 100-au window',
  active_selection='Inherited unchanged from the prior 4000-au comparison; all 10 orbitals remain visible and have metrics',
  independent_runs=rows,qm=dict(f224=qmm,f192=qm192m,max_f192_f224_difference=float(np.max(abs(qm-qm192))),
                             grid_vs_continuous_max_difference=float(np.max(abs(grid-qm)))),
  same_budget_sobol=dict(metrics=sobolscore,resources=sobolm),
  historical_v093=dict(continuous_qm=score(qm,old,initial),historical_grid_qm=score(grid,old,initial),
                      forward_paths=2000000,mpi_ranks=384,max_rank_peak_rss_kib=11452,wall_seconds=None,
                      note='Historical v093 uses 512 backward replicas and dt=0.5; path counts are not equivalent units to new paired cubature samples. Wall time unavailable.'),
  source_data=[gridsource,oldsource],kernel_validation=json.loads((DEST/'kernel_validation.json').read_text()),
  limitations=['The near-stationary orbitals do not pass relative Q>10.','No 4000-au or large-system trajectory accuracy is established.',
               'No small-space reference wavefunction is used in the Poisson estimator. QM is an external validation comparator.',
               'The selected estimator rebuilds scaled paths at each observation time. It does not retain the original full-prefix reuse strategy.',
               'The lattice samples within each randomized run are correlated; three seed results are reported individually, never averaged.',
               'Resource values include Python runtime. Full QM is faster and more accurate at this small size.',
               'No BSD/AAA, HEOM, Filinov filter, or Inchworm cache is implemented in this stage.'])
 (DEST/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
 with (DEST/'per_orbital_metrics.csv').open('w',newline='',encoding='utf-8-sig') as f:
  writer=csv.writer(f);writer.writerow(['run','seed','orbital','active','QM_signal_RMS','RMSE','Q_continuous_QM','Q_grid_QM'])
  for row in rows:
   for j in range(10):writer.writerow([row['case'],row['seed'],j,j in ACTIVE,row['continuous_qm']['signal_rms_per_orbital'][j],row['continuous_qm']['rmse_per_orbital'][j],row['continuous_qm']['q_per_orbital'][j],row['historical_grid_qm']['q_per_orbital'][j]])

 fig,axes=plt.subplots(5,2,figsize=(14,15),sharex=True)
 for j,ax in enumerate(axes.flat):
  ax.plot(t,qm[:,j]-initial[j],color='#20252b',lw=1.8,label='完整 QM（连续时间）',zorder=5)
  for i,v in enumerate(runs):ax.plot(t,v[:,j]-initial[j],color=COLORS[i],lw=1.05,alpha=.85,label=f'泊松：种子 {1160101+i}')
  qs=[r['continuous_qm']['q_per_orbital'][j] for r in rows]
  tag='近静止轨道，相对拟合未通过' if j in WEAK else '活跃轨道'
  ax.set_title(f'轨道 {j}｜{tag}\n三次独立 Q：'+', '.join(f'{v:.3g}' for v in qs),fontsize=10)
  if j in WEAK:ax.set_facecolor('#fff4ed')
  ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
  ax.set_ylabel('占据数变化 Δn');ax.grid(alpha=.18);ax.set_xlim(0,100)
 for ax in axes[-1]:ax.set_xlabel('时间 / a.u.')
 handles,labels=axes[0,0].get_legend_handles_labels()
 fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.935),ncol=4,frameon=False,fontsize=9)
 fig.suptitle('5 电子 / 10 轨道｜低内存泊松：0–100 a.u.\nη=3e-7 Ha²，wc=2.7 eV，εd=−17.194875 eV，r=√η；每个非零跳跃层 4096 个采样点',fontsize=13,y=.995)
 fig.text(.5,.012,'原始独立运行；无种子平均、平滑或裁剪。橙色底图显示变化极小轨道的不足，不能用初始占据数 1 掩盖误差。',ha='center',fontsize=9)
 fig.tight_layout(rect=(.02,.025,.99,.915),h_pad=1.6)
 fig.savefig(DEST/'all_orbitals_100au.png',dpi=160);fig.savefig(DEST/'all_orbitals_100au_preview.png',dpi=85);plt.close(fig)

 fig,axes=plt.subplots(3,2,figsize=(13,9),sharex=True)
 for j,ax in zip(ACTIVE,axes.flat):
  ax.plot(t,qm[:,j]-initial[j],color='#20252b',lw=2,label='完整 QM')
  for i,v in enumerate(runs):ax.plot(t,v[:,j]-initial[j],color=COLORS[i],lw=1,alpha=.85,label=f'种子 {1160101+i}')
  ax.set_title(f'轨道 {j}');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False);ax.grid(alpha=.18);ax.set_ylabel('Δn');ax.set_xlim(0,100)
 for ax in axes[-1]:ax.set_xlabel('时间 / a.u.')
 fig.suptitle('活跃轨道：三次独立 100 a.u. 结果与完整 QM\n泊松集体跳跃＋分层随机格点＋带雅可比的时间变量变换',fontsize=13,y=.995)
 handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.93),ncol=4,frameon=False)
 fig.tight_layout(rect=(.015,.02,.99,.90));fig.savefig(DEST/'active_orbitals_100au.png',dpi=160);fig.savefig(DEST/'active_orbitals_100au_preview.png',dpi=85);plt.close(fig)

 fig,axes=plt.subplots(1,3,figsize=(15,4.5))
 for i,row in enumerate(rows):axes[0].plot(range(10),row['continuous_qm']['q_per_orbital'],'o-',color=COLORS[i],label=str(row['seed']))
 axes[0].axhline(10,color='#555',ls='--',lw=1);axes[0].set(yscale='log',xlabel='轨道编号',ylabel='Q（越高越好）',title='全部轨道的相对精度',xticks=range(10));axes[0].grid(alpha=.18)
 for i,v in enumerate(runs):
  err=np.max(abs(v-qm),axis=1);mask=err>0;axes[1].plot(t[mask],err[mask],color=COLORS[i],label=str(1160101+i))
 axes[1].set(yscale='log',xlabel='时间 / a.u.',ylabel='最大绝对占据误差',title='每个时刻覆盖全部 10 轨道');axes[1].grid(alpha=.18)
 vals=[sobolscore['worst_active_q']]+[r['continuous_qm']['worst_active_q'] for r in rows]
 labels=['Sobol\n同采样预算','格点\n种子1','格点\n种子2','格点\n种子3']
 axes[2].bar(labels,vals,color=['#888']+COLORS);axes[2].set(yscale='log',ylabel='最差活跃轨道 Q',title='采样优化的对照');axes[2].axhline(10,color='#555',ls='--',lw=1)
 for i,v in enumerate(vals):axes[2].text(i,v*1.08,f'{v:.2f}',ha='center',fontsize=9)
 axes[2].set_ylim(min(vals)*.6,max(vals)*2);axes[2].grid(axis='y',alpha=.18)
 fig.suptitle('100 a.u. 验证｜近静止轨道 1–4 尚未达到相对精度目标',fontsize=13)
 fig.tight_layout(rect=(0,0,1,.92));fig.savefig(DEST/'accuracy_comparison.png',dpi=160);fig.savefig(DEST/'accuracy_comparison_preview.png',dpi=85);plt.close(fig)

 # Discrete correlation functions are retained exactly, following the new
 # paper's requirement to control bath representation separately from dynamics.
 m=qmm['model'];energies=np.asarray(m['energies'])[1:];weights=np.asarray(m['coupling'])[1:]**2
 f=np.array([j in m['occupied'] for j in range(1,m['norb'])],float)
 cp=np.exp(1j*t[:,None]*energies)@(weights*f);cm=np.exp(-1j*t[:,None]*energies)@(weights*(1-f))
 np.savetxt(DEST/'exact_finite_bath_correlations.dat',np.c_[t,cp.real,cp.imag,cm.real,cm.imag],header='time_au Re_C_plus Im_C_plus Re_C_minus Im_C_minus; exact discrete sum, no BSD approximation')
 report=['# 100 a.u. 低内存泊松阶段验证（2026-09-09）','',
 '已完成 5 电子 / 10 轨道、0–100 a.u. 的三次独立运行。活跃轨道 0、5–9 均通过 Q>10；近静止轨道 1–4 的相对拟合尚未通过。泊松求解过程没有读取 QM 或小空间参考波函数。','',
 '## 结果','',
 '| 独立种子 | 最差活跃轨道 Q（连续 QM） | 最差活跃轨道 Q（历史网格 QM） | 全轨道最大绝对误差 | 墙钟秒 | CPU 秒 | 峰值 RSS / MiB |',
 '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
 for row in rows:report.append(f"| {row['seed']} | {row['continuous_qm']['worst_active_q']:.4f} | {row['historical_grid_qm']['worst_active_q']:.4f} | {row['continuous_qm']['max_absolute_error_all_orbitals']:.5e} | {row['wall_seconds']:.2f} | {row['cpu_seconds']:.2f} | {row['peak_rss_mib']:.2f} |")
 report+=['',f"同采样预算 Sobol 对照的最差活跃 Q={sobolscore['worst_active_q']:.4f}，耗时 {sobolm['wall_seconds']:.2f} 秒，峰值 RSS={sobolm['peak_rss_bytes']/2**20:.2f} MiB。两个求解器都使用相同的集体 Slater、跳跃次数分层、中心化观测量和连续时间传播；格点版本另外采用随机平移格点与带雅可比的三次时间变换。",'',
 f"完整 QM（F224）耗时 {qmm['wall_seconds']:.2f} 秒、峰值 RSS={qmm['peak_rss_bytes']/2**20:.2f} MiB。本小体系上 QM 更快、更准确；低内存泊松的价值是避免在电子数和轨道数增大时存储组合数规模的波函数，本阶段尚未证明大空间或 4000 a.u. 的实际优势。",'',
 f"F192/F224 完整 QM 的最大差异为 {np.max(abs(qm-qm192)):.5e}。历史 dt=0.5 网格 QM 与连续时间 QM 的最大差异为 {np.max(abs(grid-qm)):.5e}，故并列报告两个基准；不能把算法时间离散误差和抽样误差混为一谈。",'',
 'Q=RMS(QM占据变化)/RMS(泊松−QM)，在完整的 100 a.u. 窗口逐轨道计算。活跃集合沿用历史标准；所有 10 轨道都画出、给出 Q，并进入最大绝对误差。没有用大初始占据数掩盖误差，没有平均随机种子、平滑、裁剪或删除异常点。', '',
 '## 物理与数值参数','',
 '- 原 4000 a.u. 模型：η=3e-7 Ha²，wc=2.7 eV，EF=−4.5 eV，Eb=−10 eV，εd=−0.6319 Ha，均匀实耦合 sqrt(η/9)。能级列表见 summary.json。',
 '- 质量 14583.1067146087，核频率 0.0036749323758566211 a.u.，两谐势位移 2 a.u.；初态占据轨道 0–4，核相干态位于位移势最低点。',
 '- 连续时间 Poisson 速率 r=sqrt(η)=0.0005477225575 a.u.⁻¹；每 1 a.u. 输出一次，1 a.u. 是观察间隔而非传播步长。',
 '- 每腿跳跃次数分为 1、2、3、≥4；包含完整 ≥4 条件泊松尾部。bra/ket 交换对称性合并后有 10 类，其中 2 类因奇偶性严格为零，剩余 8 类各 4096 点，共 32768 个积分点/观察时刻。随机数在时刻间复用；不是把时刻或同一格点内样本当作独立重复。',
 '- 数值事件缓冲为每腿 24；超过缓冲会报错，不丢弃事件、不把高阶贡献设为零。批量 128，只保存当前批量状态。内存含 Python、NumPy、SciPy 运行时。',
 '- 原型限制：固定粒子数、星形单体电子耦合、按杂质占据切换的位移谐振子、初始电子占据基态。未验证一般电子双体作用或任意初态。', '',
 '## 当前算法及正确性','',
 '每条电子路径用 Norb×Nel 的 Slater 矩阵表示，泊松事件作用整个星形跳跃算子，解析合并浴轨道标签。核路径用一个复相干态参数及相位表示。QR 和非正交行列式矩阵元保留费米子符号。', '',
 '测量中心化算符 n_j−n_j(0)。由于初态是占据基矢，任意一腿零跳跃时，该中心化贡献严格为零，可排除这些样本并补偿相应泊松概率。使用幺正演化的已知范数 1，不使用带噪声的分母。这个恒等式不需要参考解。', '',
 '选中的 scaled-time 方案在每个观察时刻 t 构造 τ=tφ(u)，φ(u)=3u²−2u³，对每个事件乘 φ′(u)=6u(1−u)。随机平移格点中的每一点边缘上仍均匀分布；雅可比补偿使积分不变。这是积分变量变换，不是给量子路径附加随机相位或衰减因子。', '',
 '最初试验的整段路径前缀复用噪声较大，因此最终选择逐观察时刻重新计算缩放路径，以更多计算换更低噪声和低内存。程序保留 prefix 方式用于对照，不宣称选中方案每条轨迹只传播一次。', '',
 '本阶段检查：280 次批量路径比较（覆盖 3/1、4/2、10/5、30/15、50/25、50/1、50/49）；最大状态差 1.61e-15、矩阵元差 1.18e-14；时间变换矩积分与条件泊松尾概率检查通过。这些大尺寸内核检查不是完整大体系动力学验证。', '',
 '## 新论文与前两篇的实际关联','',
 '- Dan 等，Efficient low temperature simulations for fermionic reservoirs with the hierarchical equations of motion method: Application to the Anderson impurity model，[arXiv:2211.04089v2](https://arxiv.org/abs/2211.04089v2)，Phys. Rev. B 107, 195429 (2023)。第 II C 节用 AAA/BSD 压缩谱函数，同时保持粒子/空穴关联的配对结构；第 IV 节说明近零温精度对应有限时间范围。本轮保留原离散浴，并输出其精确关联函数，未把谱函数换成更容易拟合的模型，也未实现 HEOM/AAA。',
 '- Church 等，[arXiv:1705.06617](https://arxiv.org/abs/1705.06617)：启发前后轨迹结构与振荡积分处理。未使用会改变近似程度的 Filinov 过滤，采用带雅可比的积分变量变换。',
 '- Cohen 等，[arXiv:1510.03534](https://arxiv.org/abs/1510.03534)：传播信息复用是长期候选。本轮没有 Inchworm 图重求和或两时间传播缓存。', '',
 '原论文 PDF 不上传；仅上传本项目的原创分析、链接、代码和数值结果。', '',
 '## 复现','',
 '```powershell',
 'python tests/test_poisson_slater_stream.py',
 'python experiments/poisson_slater_stream.py --method qm --nfock 224 --out scan_slater_stream_100au_20260909/reproduce_qm',
 'python experiments/poisson_slater_stream.py --method stratified-lattice --time-sampling scaled --pairs 4096 --seed 1160101 --out scan_slater_stream_100au_20260909/reproduce_seed1',
 'python analyze_slater_stream_100au.py','```','',
 '每次 final_* 目录保存完整运行命令、源码快照、代码 SHA256、原始占据数及资源记录。CSV 保存所有轨道、各个种子的指标。', '',
 '下一步优先解决近静止轨道的高阶相消噪声，再测试更长时间及更大电子空间。本阶段没有新的 4000 a.u. 数据。']
 (DEST/'README.md').write_text('\n'.join(report)+'\n',encoding='utf-8')
 print(json.dumps(dict(independent_runs=rows,sobol_worst_active_q=sobolscore['worst_active_q'],qm_convergence=np.max(abs(qm-qm192))),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
