"""Directly compare the reference alone and existing Poisson-corrected outputs."""
from pathlib import Path
import json,hashlib,re
import numpy as np
from plot_best_orbitals_20260908 import plt
from analyze_v113_large_4000 import load_reference,validate_grid
from analyze_v114_small_4000 import metrics
from analyze_versions_4000 import bath
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'scan_reference_only_ablation_20260908';OUT.mkdir(exist_ok=True)

def resources(folder):
 log=(folder/'program.out').read_text();time=(folder/'time.txt').read_text()
 def number(key,text):return float(re.search(r'\b'+key+r'=([0-9.eE+-]+)',text)[1])
 return {'wall_seconds':number('wall_seconds',time),'max_rank_peak_mib':number('max_rank_peak_rss_kb',log)/1024,'sum_rank_peaks_gib':number('sum_rank_peak_rss_kb',log)/1024**2}

def main():
 small=ROOT/'scan_v114_signed_csr_20260906';ref_folder=small/'647717'
 qm_path=ROOT/'scan_v113_two_buffer_20260905/647656/ahm-qm-s10-n5.dat'
 qm=np.loadtxt(qm_path);validate_grid(qm,10,8000,.5)
 ref_path=ref_folder/'reference-observables.dat';ref,info=load_reference(ref_path,10,5,8000,.5)
 assert info['reference_states']==246 and info['fock_levels']==384
 sources=[qm_path,ref_path,ref_folder/'config.txt',ref_folder/'program.out',ref_folder/'time.txt']
 rows=[{'case':'reference_only','job':'647717','paths_added':0,'resources':resources(ref_folder),**metrics(qm,ref[:,4:])}]
 values=[ref[:,4:]]
 for job in ['647718','647719','647720']:
  folder=small/job;path=folder/'ahm-sepmb-s10-n5-1000000.dat';data=np.loadtxt(path);validate_grid(data,10,8000,.5)
  assert np.array_equal(bath(path),bath(qm_path))
  cfg=dict(line.split('=',1) for line in (folder/'config.txt').read_text().splitlines() if '=' in line)
  assert cfg['SEP_MB_REFERENCE_CACHE'].split('/')[-2]=='647717'
  assert int(cfg['SEP_MB_REFERENCE_DISTANCE'])==3 and int(cfg['SEP_MB_REFERENCE_FOCK_STATES'])==384
  rr=resources(folder);rr['first_total_wall_seconds']=rr['wall_seconds']+rows[0]['resources']['wall_seconds']
  rows.append({'case':'reference_plus_poisson','job':job,'seed':cfg['AHM_SEED'],'paths_added':1000000,'resources':rr,**metrics(qm,data[:,4:])})
  values.append(data[:,4:]);sources += [path,folder/'config.txt',folder/'program.out',folder/'time.txt']
 result={'scope':'5 electrons / 10 orbitals, 0-4000 au, dt=0.5, D3/F384, same reference and grid-QM comparator',
  'conclusion':'For this tested case, the reference alone is substantially more accurate against grid QM than any of the three million-path corrected runs. No practical accuracy benefit from the Poisson correction is demonstrated.',
  'limits':'Grid-QM agreement is not an exact continuum error bound. The result does not prove large-space reference accuracy, universal uselessness of Poisson, or the source of all residual errors. Reference producer wall time includes its small pilot; consumer time is additional. Sum of rank peaks is not simultaneous RSS.',
  'reference_info':info,'rows':rows,'source_hashes':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
 (OUT/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False),encoding='utf-8')
 colors=['#173f4b','#c6774c','#909849','#926ca4'];labels=['参考单独使用','参考 + 泊松，种子1144001','参考 + 泊松，种子1144002','参考 + 泊松，种子1144003']
 fig,axes=plt.subplots(1,2,figsize=(14,5.6));fig.subplots_adjust(left=.075,right=.98,bottom=.2,top=.77,wspace=.29)
 for r,v,color,label in zip(rows,values,colors,labels):
  ww=r['windows'];axes[0].plot([(w['start']+w['stop'])/2 for w in ww],[w['Q_weakest_active'] for w in ww],color=color,lw=1.5 if r['case']=='reference_only' else 1,label=label)
  error=abs(v-qm[:,4:]).max(axis=1)
  axes[1].plot(qm[:,0],np.where(error>0,error,np.nan),color=color,lw=1.5 if r['case']=='reference_only' else .7,label=label)
 axes[0].axhline(10,color='#777',ls=':',lw=.8);axes[0].set(title='最差活跃轨道窗口 Q：越高越好',xlabel='100 a.u. 窗口中心 / a.u.',ylabel='Q（对照相同网格 QM）',yscale='log',xlim=(0,4000))
 axes[1].set(title='每个时刻的全轨道最大绝对误差：越低越好',xlabel='时间 / a.u.',ylabel='占据数绝对误差',yscale='log',xlim=(0,4000))
 for ax in axes:ax.grid(alpha=.16);ax.set_xticks([0,1000,2000,3000,4000])
 fig.suptitle('5电子 / 10轨道，4000 a.u.：泊松修正是否改善参考结果？',fontsize=17,y=.975,weight='bold')
 fig.text(.5,.885,'同一 D3/F384 参考（246/252态）；wc=2.7 eV，η=3e-7，dt=0.5；每次修正100万路径、64进程',ha='center',fontsize=10)
 fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.055),ncol=2,frameon=False,fontsize=9)
 fig.text(.075,.015,'Q覆盖活跃轨道0、5–9；误差图覆盖全部10轨道。零误差点不显示在对数轴上；未修改原始输出。',fontsize=8.5,color='#59676d')
 fig.savefig(OUT/'reference_vs_corrected.png',dpi=160,facecolor='white')
 fig.savefig(OUT/'reference_vs_corrected_preview.png',dpi=85,facecolor='white');plt.close(fig)
 lines=['# 参考单独使用与泊松修正的直接对照（2026-09-08）','',
 '**在现有 5电子/10轨道、4000 a.u. 测试中，直接使用 D3/F384 参考更准确，也省去了约一小时的百万路径修正。此次对照没有证明泊松修正有实用收益。**','',
 '比较使用同一个参考生产作业647717、同一个固定网格QM基准647656，以及原有三个独立百万路径结果。重新核对8001点时间网格、缓存参数与校验和及浴参数；未平均种子、平滑或裁剪。参考原始观测量先按缓存规范归一化，再与原有完整输出进行相同指标比较。','',
 '| 输出 | 最差活跃轨道窗口 Q | 全轨道最大绝对误差 | 本阶段墙钟秒 |','| --- | ---: | ---: | ---: |']
 for r in rows:lines.append(f"| {'参考单独使用' if r['case']=='reference_only' else '参考+泊松 '+r['seed']} | {r['Q_against_QM_weakest_active_orbital_window']:.6f} | {r['max_orbital_error']:.9e} | {r['resources']['wall_seconds']:.2f} |")
 lines += ['',
 'Q = RMS(QM占据数变化) / RMS(输出与QM的差异)，每100 a.u.窗口取活跃轨道0、5–9的最小值。越大越好。近静止轨道1–4包含在全轨道绝对误差中。参考生产24.08秒包括当时的小样本试跑，是现有参考生产作业的实测成本；三个修正耗时是复用参考后额外支付的时间。全程均64进程，不与单进程QM简单比较速度。','',
 '参考空间246/252态接近完整，现有参数下参考与QM已高度一致。百万路径估计带来的随机波动和可能的数值偏差，超过当前参考偏差；不能把“混合输出Q>10”解释为“泊松使参考更准确”。此比较不能单独分解随机方差、采样偏差与时间离散误差。','',
 '参考近似的优势是在可控内存内快速计算；泊松修正的研究价值，是在参考受内存限制、截断误差已不可忽略时，以增加采样时间降低这一误差。只有在固定内存预算下，校正后的误差实际小于参考单独误差，且收益可重复，才能证明这一价值。当前尚无这样的完整大体系实证。','',
 '15/30与25/50没有完整QM真值。相邻D/F参考一致只构成截断收敛证据，不能证明大体系参考已准确；现有修正结果也不能证明它纠正了参考。不能将本次小体系结论直接推广为所有情况下都无需随机修正。','',
 '后续评价必须包含：参考单独、参考加修正及可核验基准三者；比较固定内存下的实际误差与总时间。应先在可用完整QM的小体系设置明显的参考截断误差，再检验修正能否稳定恢复遗漏贡献，同时独立检查有限时间步和抽样偏差。未证明前，将“符合目标的研究设计”与“已验证有效的算法收益”分开陈述。','',
 '![参考单独与修正后对比](reference_vs_corrected.png)','',
 '[完整指标与来源校验值](validation.json)。重现：在仓库根目录执行 `python analyze_reference_only_ablation.py`。这里只重算已有数据，不修改求解器或启动新采样。','',
 '阶段回退点：[v1.14-reference-only-ablation](https://github.com/nihility1shadow/qm_project/tree/v1.14-reference-only-ablation)。']
 (OUT/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
 print(json.dumps([{k:r[k] for k in ['case','job','Q_against_QM_weakest_active_orbital_window','max_orbital_error']} for r in rows],indent=2))
if __name__=='__main__':main()
