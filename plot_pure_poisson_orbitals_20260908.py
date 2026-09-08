"""Annotate existing no-reference Poisson runs; never create missing simulations."""
from pathlib import Path
import json,re,math,hashlib
import numpy as np
from plot_best_orbitals_20260908 import plt,Line2D,hh
from analyze_v113_large_4000 import validate_grid
from analyze_versions_4000 import bath
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'pure_poisson_figures_20260908';OUT.mkdir(exist_ok=True)
CASES=[
 dict(name='5e10o_v093_2500au',job='646440',version='v0.93',nel=5,norb=10,tmax=2500,paths=50000000,ranks=384,seed=2500824,
      folder='scan_v093_poisson_t2500_20260822/final_50m_job_646440',config='run-info.txt',
      qm='scan_v093_poisson_t2500_20260822/pilots/qm_646419_eta03/ahm-qm-s10-n5.dat',role='历史单次主结果',wall_seconds=53120,
      wall_source='scan_v093_poisson_t2500_20260822/final_report.md',measurement='自适应直接采样点 414；其余输出点由原程序插值',
      sampling='后向副本 B = 256→512，随时间线性增加；前向轨道条件求和；精确后向低阶跳跃 = 4'),
 dict(name='5e10o_v093_4000au',job='647800',version='v0.93',nel=5,norb=10,tmax=4000,paths=10000,ranks=64,seed=1158001,
      folder='scan_versions_4000_20260907/647800',config='config.txt',qm='scan_v113_two_buffer_20260905/647656/ahm-qm-s10-n5.dat',
      role='同条件低样本版本基线，未达标',measurement='每 4 步直接采样，中间由原程序插值',
      sampling='后向副本 B = 16；前向轨道采样；精确后向低阶跳跃 = 0'),
 dict(name='25e50o_v114_noref_4000au',job='647782',version='v1.14（关闭参考）',nel=25,norb=50,tmax=4000,paths=1000,ranks=64,seed=1155099,
      folder='scan_v115_generality_20260906/647782',config='config.txt',qm=None,role='关闭参考的低样本试跑，严重失效',
      measurement='每 4 步直接采样，中间由原程序插值',sampling='后向副本 B = 16；速率倍率 = 1.5；RQMC = 4；按步分层开启')]

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def val(name,txt):return float(re.search(r'\b'+re.escape(name)+r'=([0-9.eE+-]+)',txt)[1])
def main():
 results=[]
 for c0 in CASES:
  c=dict(c0);folder=ROOT/c['folder'];raw=folder/f"ahm-sepmb-s{c['norb']}-n{c['nel']}-{c['paths']}.dat"
  data=np.loadtxt(raw);validate_grid(data,c['norb'],2*c['tmax'],.5)
  cfg=dict(line.split('=',1) for line in (folder/c['config']).read_text().splitlines() if '=' in line)
  log=(folder/'program.out').read_text();headers='\n'.join(line for line in raw.read_text().splitlines() if line.startswith('#'))
  assert 'selected=stochastic-poisson' in log
  assert val('nproc',log)==c['ranks'] and val('seed_base',log)==c['seed']
  for k,v in {'AHM_WC_EV':2.7,'AHM_ETA':3e-7,'AHM_DELE_EV':-17.194874968839816,'AHM_NSTEP':2*c['tmax'],'SEP_MB_DETERMINISTIC_FOCK':0}.items():assert float(cfg[k])==v
  if c['job']=='647782':assert 'selected_distance=-1 active=0 states=0' in log and 'full_determinant_basis=skipped' in log
  else:assert '#SEP_MB_REFERENCE ' not in log # These audited v0.93 binaries predate the reference-control implementation.
  peaks=val('max_rank_peak_rss_kb',log)/1024
  if 'wall_seconds' not in c:c['wall_seconds']=val('wall_seconds',(folder/'time.txt').read_text());c['wall_source']=str((folder/'time.txt').relative_to(ROOT)).replace('\\','/')
  if 'sum_rank_peak_rss_kb=' in log:
   summed=val('sum_rank_peak_rss_kb',log)/1024**2;mem=f'进程峰值和 {summed:.3f} GiB';mem_kind='sum_of_individual_peaks_not_simultaneous_RSS'
  else:
   summed=peaks*c['ranks']/1024;mem=f'最大进程峰值 × 进程数的上界 {summed:.3f} GiB';mem_kind='maximum_rank_peak_times_rank_count_upper_bound_not_measured_sum'
  initial=np.r_[np.ones(c['nel']),np.zeros(c['norb']-c['nel'])]
  assert np.max(abs(data[0,4:]-initial))<1e-12
  q=None;qm=None;window=[];pass_until=0
  sources=[raw,folder/c['config'],folder/'program.out',ROOT/c['wall_source']]
  if c['qm']:
   qmpath=ROOT/c['qm'];qm=np.loadtxt(qmpath);validate_grid(qm,c['norb'],2*c['tmax'],.5)
   assert np.array_equal(bath(raw),bath(qmpath)),'Bath mismatch'
   signal=qm[:,4:]-initial;error=data[:,4:]-qm[:,4:];active=[0,*range(c['nel'],c['norb'])]
   for start in range(0,c['tmax'],100):
    mask=(data[:,0]>=start)&(data[:,0]<=start+100)
    qq=[float(np.sqrt(np.mean(signal[mask,i]**2))/np.sqrt(np.mean(error[mask,i]**2))) for i in active]
    window.append({'start':start,'stop':start+100,'minimum_active_Q':min(qq)})
    if pass_until==start and min(qq)>=10:pass_until=start+100
   q=min(w['minimum_active_Q'] for w in window);sources.append(qmpath)
  if c['job']=='646440':
   previous=json.loads((folder/'analysis/q_metrics.json').read_text());assert np.isclose(q,previous['minimum_strict_Q'],rtol=1e-10)
   assert pass_until==700
   status=f"前 0–700 a.u. 各 100 a.u. 窗口 Q ≥ 10；全程最差 Q = {q:.4f}，未全程达标"
  elif qm is not None:status=f"4000 a.u. 低样本对照；全程最差活跃轨道窗口 Q = {q:.3g}，未达标"
  else:status=f"无完整 QM 对照；占据数范围 {data[:,4:].min():.1f} 至 {data[:,4:].max():.1f}，非物理波动"
  sampling_line=next(line for line in headers.splitlines() if line.startswith('#sampling:'))
  probability=val('jump_probability',sampling_line)
  c.update(config_values=cfg,actual_sampling_header=sampling_line,max_rank_peak_mib=peaks,memory_aggregate_gib=summed,memory_aggregate_kind=mem_kind,
           minimum_occupation=float(data[:,4:].min()),maximum_occupation=float(data[:,4:].max()),strict_Q=q,Q10_contiguous_end_au=pass_until if qm is not None else None,
           windows=window,data_file=str(raw.relative_to(ROOT)).replace('\\','/'),sources_sha256={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in sources})
  rows=math.ceil(c['norb']/5);height=4.2+2.18*rows
  fig,axes=plt.subplots(rows,5,figsize=(16.6,height),squeeze=False)
  fig.subplots_adjust(left=.057,right=.99,bottom=1.05/height,top=1-3.55/height,hspace=.62,wspace=.36)
  fig.text(.055,1-.34/height,f"{c['nel']} 电子 / {c['norb']} 轨道 · 纯泊松 · 0–{c['tmax']} a.u.",fontsize=19,weight='bold',va='top',color='#173641')
  fig.text(.055,1-.83/height,status,fontsize=12.2,color='#a84c24',va='top',weight='bold')
  lines=[f"{c['version']}；{c['role']}；无独立确定性参考；单次运行，未平均种子",
   f"电子参数：wc = 2.7 eV；η = 3e-7（代码单位）；ΔE = −17.194875 eV；初态占据轨道 0–{c['nel']-1}",
   f"时间与路径：dt = 0.5 a.u.，{2*c['tmax']} 步；前向路径 N = {c['paths']:,}；种子 {c['seed']}；{c['ranks']} 进程",
   c['sampling'],f"测量：{c['measurement']}；名义单步事件概率 p = {probability:.10g}",
   f"资源：完整作业耗时 {hh(c['wall_seconds'])}；最大进程峰值 {peaks:.3f} MiB；{mem}",
   f"作业 {c['job']}；无参考生成成本；展示原始程序输出，未追加归一化、平滑或越界裁剪"]
  for i,line in enumerate(lines):fig.text(.055,1-(1.28+.245*i)/height,line,fontsize=9.8,va='top',color='#34464e')
  handles=[Line2D([0],[0],color='#c6633b',lw=1.5,label='无参考泊松：单次完整输出')]
  if qm is not None:handles.append(Line2D([0],[0],color='#1f4054',lw=1.2,ls='--',label='匹配物理参数的网格 QM'))
  fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.05,1-3.04/height),ncol=2,frameon=False,fontsize=10)
  for i,ax in enumerate(axes.flat):
   if i>=c['norb']:ax.set_visible(False);continue
   ax.plot(data[:,0],data[:,i+4]-initial[i],color='#c6633b',lw=.72,alpha=.9,zorder=2)
   if qm is not None:ax.plot(qm[:,0],qm[:,i+4]-initial[i],color='#1f4054',lw=.85,ls='--',zorder=3)
   ax.axhline(0,color='#c6ced1',lw=.4,zorder=0)
   ax.set_title(f'轨道 {i}   初始占据 {int(initial[i])}',fontsize=10,pad=6)
   ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
   ax.tick_params(labelsize=8.2);ax.yaxis.get_offset_text().set_size(8.2)
   ax.set_xlim(0,c['tmax']);ax.set_xticks(np.arange(0,c['tmax']+1,500 if c['tmax']==2500 else 1000));ax.grid(alpha=.16,lw=.55)
   if i%5==0:ax.set_ylabel(r'$\Delta n_i=n_i(t)-n_i(0)$',fontsize=9)
   if i//5==rows-1:ax.set_xlabel('时间 / a.u.',fontsize=9)
  fig.text(.055,.20/height,'每个子图纵轴独立缩放。纯泊松指无确定性参考；仍保留条件求和/分层。不同样本数与时间范围不可直接排名。',fontsize=9,color='#59676d',va='bottom')
  c['figure']=f"pure_poisson_figures_20260908/{c['name']}_all_orbitals.png"
  fig.savefig(ROOT/c['figure'],dpi=160,facecolor='white')
  fig.savefig(OUT/f"{c['name']}_preview.png",dpi=80,facecolor='white');plt.close(fig)
  results.append(c)
 (OUT/'parameters_and_validation.json').write_text(json.dumps({'notice':'Existing no-reference Poisson runs only. No verified 15e30o no-reference long-time output found. These are different budgets, not a fair best-method ranking. No new simulations were launched.',
  'hash_notice':'SHA256 values refer to the local source bytes at generation time; Git checkout line-ending conversion can change byte hashes without changing numbers.', 'runs':results},indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
 print(json.dumps([{k:c[k] for k in ['name','job','strict_Q','Q10_contiguous_end_au','wall_seconds','max_rank_peak_mib','minimum_occupation','maximum_occupation']} for c in results],indent=2))
if __name__=='__main__':main()
