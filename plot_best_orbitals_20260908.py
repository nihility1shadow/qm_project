"""Build annotated all-orbital figures from audited, completed 4000-au runs.

Selection is retrospective. Large-space reference agreement is NOT QM accuracy.
The script never modifies raw observations, averages seeds, clips, or smooths.
"""
from pathlib import Path
from functools import lru_cache
import hashlib,json,math,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from analyze_v113_large_4000 import load_reference,validate_grid,windows

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'best_orbital_figures_20260908'
OUT.mkdir(exist_ok=True)
SMALL=ROOT/'scan_v114_signed_csr_20260906'
MID=ROOT/'scan_v113_two_buffer_20260905'
LARGE=ROOT/'scan_v115_generality_20260906'
SCREEN=ROOT/'scan_rate_4000_20260907'
FONT=Path('C:/Windows/Fonts/msyh.ttc')
if FONT.exists():
    font_manager.fontManager.addfont(str(FONT))
    plt.rcParams['font.family']=font_manager.FontProperties(fname=str(FONT)).get_name()
plt.rcParams.update({'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,
                     'font.size':10,'axes.titleweight':'medium','svg.fonttype':'none'})

def read_json(path):return json.loads(path.read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def relative(path):return path.relative_to(ROOT).as_posix()
@lru_cache(None)
def reference(path,norb,nel):return load_reference(Path(path),norb,nel,8000,.5)
def config(folder):return dict(line.split('=',1) for line in (folder/'config.txt').read_text().splitlines() if '=' in line)
def worst_q(ww):return min(w['Q_weakest_active'] for w in ww if w['Q_weakest_active'] is not None)
def hh(seconds):
    seconds=round(seconds);h,seconds=divmod(seconds,3600);m,s=divmod(seconds,60)
    return f'{h}小时{m:02d}分{s:02d}秒' if h else f'{m}分{s:02d}秒'

def all_candidates():
    records=[]
    for row in read_json(SMALL/'small_4000_validation.json')['jobs']:
        records.append({'group':'5e10o','job':row['job'],'folder':relative(SMALL/row['job']),
            'electrons':5,'orbitals':10,'paths':1000000,'reference_job':'647717','reference_folder':relative(SMALL/'647717'),
            'score':row['Q_against_QM_weakest_active_orbital_window'],'score_kind':'QM_Q','version':'v1.14'})
    mid_report=read_json(MID/'large_4000_validation.json')
    assert mid_report['status']=='all seven runs downloaded and validated'
    for job,refjob in [('647704','647702'),('647705','647702'),('647706','647702'),('647707','647703')]:
        data=np.loadtxt(MID/job/'ahm-sepmb-s30-n15-1000000.dat');validate_grid(data,30,8000,.5)
        ref,_=reference(str(MID/refjob/'reference-observables.dat'),30,15)
        initial=np.r_[np.ones(15),np.zeros(15)]
        score=worst_q(windows(data[:,0],ref[:,4:]-initial,data[:,4:]-ref[:,4:],[0,*range(15,30)]))
        records.append({'group':'15e30o','job':job,'folder':relative(MID/job),'electrons':15,'orbitals':30,'paths':1000000,
            'reference_job':refjob,'reference_folder':relative(MID/refjob),'score':score,'score_kind':'reference_diagnostic','version':'v1.13'})
    for row in read_json(SCREEN/'sampling_validation.json')['runs']:
        records.append({'group':'15e30o','job':row['job'],'folder':relative(SCREEN/row['job']),'electrons':15,'orbitals':30,
            'paths':int(row['paths']),'reference_job':'647702','reference_folder':relative(MID/'647702'),
            'score':row['worst_active_orbital_window_ratio'],'score_kind':'reference_diagnostic','version':'v1.14'})
    for row in read_json(LARGE/'long_4000_validation.json')['runs']:
        if row['job'] not in ['647777','647778','647779','647780']:continue
        refjob=row['reference_source_job']
        paths=10000 if row['job']=='647778' else 1000
        records.append({'group':'25e50o','job':row['job'],'folder':relative(LARGE/row['job']),'electrons':25,'orbitals':50,
            'paths':paths,'reference_job':refjob,'reference_folder':relative(LARGE/refjob),
            'score':row['worst_active_orbital_reference_difference_ratio'],'score_kind':'reference_diagnostic','version':'v1.14'})
    return records

def resources(folder):
    log=(folder/'program.out').read_text()
    timing=(folder/'time.txt').read_text()
    return {'wall_seconds':float(re.search(r'wall_seconds=([0-9.]+)',timing)[1]),
        'max_rank_peak_mib':int(re.search(r'max_rank_peak_rss_kb=(\d+)',log)[1])/1024,
        'sum_rank_peaks_gib':int(re.search(r'sum_rank_peak_rss_kb=(\d+)',log)[1])/1024**2}

def make_figure(case):
    nel,norb=case['electrons'],case['orbitals'];folder=ROOT/case['folder'];ref_folder=ROOT/case['reference_folder']
    raw_path=folder/f"ahm-sepmb-s{norb}-n{nel}-{case['paths']}.dat"
    data=np.loadtxt(raw_path);validate_grid(data,norb,8000,.5)
    settings=config(folder)
    for key,value in {'AHM_WC_EV':2.7,'AHM_ETA':3e-7,'AHM_DELE_EV':-17.194874968839816,'AHM_NSTEP':8000,
        'SEP_MB_TMAX':4000,'SEP_MB_PATH_LOCAL_BASIS':1,'SEP_MB_DETERMINISTIC_FOCK':0,
        'SEP_MB_EXACT_ORBITALS':0,'SEP_MB_MEASURE_STRIDE':4,'SEP_MB_RATE_SCALE':1.5,
        'SEP_MB_BACK_REPLICAS':16,'SEP_MB_RQMC_REPLICATES':4}.items():
        assert float(settings[key])==value,(case['job'],key)
    ref_path=ref_folder/'reference-observables.dat';ref,info=reference(str(ref_path),norb,nel)
    assert int(settings['SEP_MB_REFERENCE_FOCK_STATES'])==info['fock_levels']
    assert settings['SEP_MB_REFERENCE_CACHE'].split('/')[-2]==case['reference_job']
    key=ref_path.read_text().splitlines()[0].split()
    mass,omega,displacement,alpha_real,alpha_imag=map(float,key[13:18])
    x_init=alpha_real*math.sqrt(2/(mass*omega));p_init=alpha_imag*math.sqrt(2*mass*omega)
    assert abs(x_init-2)<1e-12 and p_init==0
    initial=np.r_[np.ones(nel),np.zeros(norb-nel)]
    # Both curves use the identical known initial occupation; original output is retained.
    initial_in_key=np.array([int(key[22+3*i]) for i in range(norb)])
    assert np.array_equal(initial,initial_in_key)
    if case['score_kind']=='QM_Q':
        comparator_path=MID/'647656/ahm-qm-s10-n5.dat'
        comparator=np.loadtxt(comparator_path);validate_grid(comparator,norb,8000,.5)
        comparator_label='完整组态网格 QM（数值基准）'
        status=f"已通过既定 QM 对照指标 | 最差活跃轨道窗口 Q = {case['score']:.3f}"
        status_color='#216b55'
    else:
        comparator_path=ref_path;comparator=ref
        comparator_label='有界确定性参考（不是完整 QM）'
        status=f"当前候选，尚未达标 | 相对参考诊断 = {case['score']:.3g}（不是 QM 拟合 Q）"
        status_color='#a84c24'
    verified_score=worst_q(windows(data[:,0],comparator[:,4:]-initial,data[:,4:]-comparator[:,4:],[0,*range(nel,norb)]))
    assert np.isclose(verified_score,case['score'],rtol=1e-10)
    rr,pr=resources(folder),resources(ref_folder)
    case.update({'settings':settings,'reference_info':info,'resources':rr,'reference_producer_resources':pr,
        'first_total_wall_seconds':rr['wall_seconds']+pr['wall_seconds'],
        'mass_au':mass,'omega_au':omega,'nuclear_displacement_au':displacement,'x_init_au':x_init,'p_init_au':p_init,
        'data_file':relative(raw_path),'comparison_file':relative(comparator_path),
        'minimum_occupation':float(data[:,4:].min()),'maximum_occupation':float(data[:,4:].max()),
        'maximum_difference_to_comparator':float(abs(data[:,4:]-comparator[:,4:]).max()),
        'source_hashes':{relative(p):sha(p) for p in [raw_path,comparator_path,ref_path,folder/'config.txt',folder/'binary.sha256']}})
    rows=math.ceil(norb/5);fig_height=4.2+2.18*rows
    fig,axes=plt.subplots(rows,5,figsize=(16.6,fig_height),squeeze=False)
    chart_top=1-3.55/fig_height;chart_bottom=1.05/fig_height
    fig.subplots_adjust(left=.057,right=.99,bottom=chart_bottom,top=chart_top,hspace=.62,wspace=.36)
    fig.text(.055,1-.34/fig_height,f'{nel} 电子 / {norb} 轨道  ·  0–4000 a.u.  ·  {case["version"]}',fontsize=19,weight='bold',va='top',color='#173641')
    fig.text(.055,1-.83/fig_height,status,fontsize=12.2,color=status_color,va='top',weight='bold')
    text_lines=[
      'Anderson–Holstein 星形模型；有界参考解 + 泊松残差；展示单次运行，未平均随机种子',
      f"电子参数：wc = 2.7 eV；η = 3e-7（代码单位）；ΔE = −17.194875 eV；初态占据轨道 0–{nel-1}",
      f"核参数：m = {mass:.8f}；ω = {omega:.12f}；Δx = 2；x0 = 2，p0 = 0（以上均为原子单位）",
      f"数值参数：dt = 0.5 a.u.，8000 步；D{settings['SEP_MB_REFERENCE_DISTANCE']} / F{info['fock_levels']}，{info['reference_states']:,} 参考态；前向路径 N = {case['paths']:,}",
      f"采样参数：B = 16；速率倍率 = 1.5；RQMC = 4；残差每 4 步直接采样，中间线性插值；种子 {settings['AHM_SEED']}",
      f"资源：{settings['SLURM_NTASKS']} 进程；本次采样 {hh(rr['wall_seconds'])}；峰值/进程 {rr['max_rank_peak_mib']:.2f} MiB，峰值和 {rr['sum_rank_peaks_gib']:.3f} GiB",
      f"首次参考另计：{hh(pr['wall_seconds'])}；峰值和 {pr['sum_rank_peaks_gib']:.3f} GiB；作业 {case['job']}，参考生产 {case['reference_job']}"
    ]
    for index,line in enumerate(text_lines):fig.text(.055,1-(1.28+.245*index)/fig_height,line,fontsize=9.8,va='top',color='#34464e')
    colors=['#c6633b','#1f4054']
    fig.legend(handles=[Line2D([0],[0],color=colors[0],lw=1.5,label='参考解 + 泊松修正（完整输出）'),
                        Line2D([0],[0],color=colors[1],lw=1.25,ls='--',label=comparator_label)],
               loc='upper left',bbox_to_anchor=(.05,1-3.04/fig_height),ncol=2,frameon=False,fontsize=10)
    t=data[:,0]
    for orbital,ax in enumerate(axes.flat):
        if orbital>=norb:ax.set_visible(False);continue
        ax.plot(t,data[:,orbital+4]-initial[orbital],color=colors[0],lw=.72,alpha=.9,zorder=2)
        ax.plot(t,comparator[:,orbital+4]-initial[orbital],color=colors[1],lw=.85,ls='--',zorder=3)
        ax.axhline(0,color='#c6ced1',lw=.4,zorder=0)
        ax.set_title(f'轨道 {orbital}   初始占据 {int(initial[orbital])}',fontsize=10,pad=6)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
        ax.tick_params(labelsize=8.2);ax.yaxis.get_offset_text().set_size(8.2)
        ax.set_xlim(0,4000);ax.set_xticks([0,1000,2000,3000,4000]);ax.grid(alpha=.16,lw=.55)
        if orbital%5==0:ax.set_ylabel(r'$\Delta n_i=n_i(t)-n_i(0)$',fontsize=9)
        if orbital//5==rows-1:ax.set_xlabel('时间 / a.u.',fontsize=9)
    fig.text(.055,.20/fig_height,'每个子图纵轴独立缩放；展示全部轨道和完整时域，未裁剪越界值、未平滑。峰值和并非同时总内存。',fontsize=9,color='#59676d',va='bottom')
    filename=f"{case['group']}_all_orbitals_4000au.png"
    fig.savefig(OUT/filename,dpi=160,facecolor='white')
    fig.savefig(OUT/filename.replace('.png','_preview.png'),dpi=80,facecolor='white');plt.close(fig)
    case['figure']=relative(OUT/filename)
    return case

def main():
    pool=all_candidates()
    selected=[]
    for group in ['5e10o','15e30o','25e50o']:
        candidates=[c for c in pool if c['group']==group]
        selected.append(make_figure(dict(max(candidates,key=lambda c:c['score']))))
    output={'time_range_au':[0,4000],'selection_notice':'Retrospective best observed weakest active-orbital 100-au window ratio in the stated candidate pool. Only 5e10o has full QM. Large-space ranking is a bounded-reference diagnostic, not actual accuracy. No seed averaging.',
        'excluded_notice':'30-orbital 100-path reference-production runs excluded: undersampled rare events can falsely resemble the reference. 50-orbital 17-electron run is a different filling and excluded. No-reference failed pilot is not the requested hybrid candidate.',
        'representation':'All orbitals, original normalized output minus known initial occupation; existing every-4-step residual interpolation retained; no added smoothing, clipping or time truncation.',
        'candidates':pool,'selected':selected}
    (OUT/'selection_and_parameters.json').write_text(json.dumps(output,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps([{'case':c['group'],'job':c['job'],'score':c['score'],'figure':c['figure'],'total_hours':c['first_total_wall_seconds']/3600} for c in selected],indent=2))

if __name__=='__main__':main()
