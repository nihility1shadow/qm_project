"""Long-time validation: replicate noise is explicitly distinguished from QM error."""
import json,re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/"scan_v111_distributed_20260905"
FIG=OUT/"figures";FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"savefig.dpi":160})
def load(j):
    return np.loadtxt(next((OUT/str(j)).glob("ahm-sepmb-*.dat")))
ids=[647613,647614,647622]
data=[load(j) for j in ids]
t=data[0][:,0]
assert all(np.array_equal(a[:,0],t) for a in data)
values=np.stack([a[:,4:] for a in data])
initial=np.r_[np.ones(15),np.zeros(15)]
active=[0]+list(range(15,30))
mean=values.mean(axis=0);single_sd=values.std(axis=0,ddof=1)
rms=lambda a:float(np.sqrt(np.mean(a*a)))
def ratio(a,b):
    return rms(a)/rms(b) if rms(b)>0 else None
rows=[]
for start in range(0,2500,100):
    ix=(t>=start)&(t<=start+100)
    orbital_q=[ratio((mean-initial)[ix,i],single_sd[ix,i]) for i in range(30)]
    rows.append({"start":start,"stop":start+100,
        "Q_repeat_active":ratio((mean-initial)[ix][:,active],single_sd[ix][:,active]),
        "min_Q_repeat_active":min(orbital_q[i] for i in active if orbital_q[i] is not None),
        "Q_repeat_per_orbital":orbital_q,
        "signal_rms_active":rms((mean-initial)[ix][:,active]),
        "single_run_noise_rms_active":rms(single_sd[ix][:,active])})
f512=load(647623)
fock_difference=np.abs(data[0][:,4:]-f512[:,4:])
resources=[]
for j in ids+[647623,647624]:
    folder=OUT/str(j);log=(folder/"program.out").read_text()
    wall=(folder/"time.txt").read_text()
    def extract(pattern):
        m=re.search(pattern,log);return float(m[1]) if m else None
    resources.append({"job":j,
        "wall_seconds":float(re.search(r"wall_seconds=([\d.]+)",wall)[1]),
        "max_rank_peak_mib":extract(r"max_rank_peak_rss_kb=(\d+)")/1024,
        "sum_rank_peak_gib":extract(r"sum_rank_peak_rss_kb=(\d+)")/1024**2,
        "reference_seconds":extract(r"reference_seconds=([\d.]+)")})
report={"metric_notice":"Q_repeat uses the RMS standard deviation of SINGLE runs across 3 seeds. It is not absolute error against QM and cannot detect common deterministic bias.",
        "jobs":ids,"paths_per_run":100000,"active_orbitals":active,
        "windows":rows,"min_Q_repeat_active":min(r["Q_repeat_active"] for r in rows),
        "min_active_orbital_Q_repeat":min(r["min_Q_repeat_active"] for r in rows),
        "fock384_vs512_matched_seed_max_orbital_difference":float(fock_difference.max()),
        "max_particle_error":float(np.max(np.abs(values.sum(axis=2)-15))),
        "min_occupation":float(values.min()),"max_occupation":float(values.max()),
        "resources":resources,"meets_Q10_every_active_orbital":all(r["min_Q_repeat_active"] >= 10 for r in rows)}
(OUT/"long_validation.json").write_text(json.dumps(report,indent=2))
fig,axes=plt.subplots(6,5,figsize=(18,16),sharex=True)
colors=["#278b98","#d68133","#875aaa"]
for orbital,ax in enumerate(axes.flat):
    for k in range(3):
        ax.plot(t,values[k,:,orbital]-initial[orbital],lw=.65,color=colors[k],alpha=.8,label=f"Seed {k+1}")
    ax.set_title(f"Orbital {orbital}",fontsize=9)
    ax.ticklabel_format(axis="y",style="sci",scilimits=(0,0),useOffset=False)
    ax.grid(alpha=.15)
axes[0,0].legend(fontsize=7)
for ax in axes[-1]:ax.set_xlabel("Time (a.u.)")
fig.suptitle("30 orbitals / 15 electrons | 2,500 a.u. | 100,000 paths PER RUN | 3 independent seeds\n"
             "52,656 bounded reference states; 384 oscillator levels; 64 MPI ranks per run",y=.995)
fig.tight_layout(rect=(0,0,1,.967));fig.savefig(FIG/"long_30_orbitals_three_seeds.png");plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(11,4.1))
centers=[(r["start"]+r["stop"])/2 for r in rows]
axes[0].semilogy(centers,[r["Q_repeat_active"] for r in rows],"-o",ms=3,label="Combined active orbitals")
axes[0].semilogy(centers,[r["min_Q_repeat_active"] for r in rows],"-s",ms=3,label="Weakest active orbital")
axes[0].axhline(10,color="#b64545",ls="--",label="Target Q = 10")
axes[0].set(xlabel="Window midpoint (a.u.)",ylabel="Single-run repeatability Q")
axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
axes[1].semilogy(t,np.maximum(fock_difference.max(axis=1),1e-18),color="#338b80")
axes[1].set(xlabel="Time (a.u.)",ylabel="Max orbital difference: 384 vs 512 levels")
axes[1].grid(alpha=.2)
fig.suptitle("Noise grows late; oscillator truncation is much smaller than sampling noise\n"
             "Repeatability Q is not a QM fitting score")
fig.tight_layout();fig.savefig(FIG/"long_q_and_fock_convergence.png");plt.close(fig)

old=np.loadtxt(OUT/"baseline_v110_30_15_t2500.dat")
fig,axes=plt.subplots(3,1,figsize=(11,8),sharex=True)
for ax,i in zip(axes,[0,15,29]):
    ax.plot(old[:,0],old[:,4+i]-initial[i],color="#a5a5a5",lw=.7,label="v1.10: 1,696 reference states")
    ax.plot(t,values[2,:,i]-initial[i],color="#258879",lw=.8,label="v1.12: 52,656 reference states")
    ax.set_title(f"Orbital {i}");ax.grid(alpha=.15)
    ax.ticklabel_format(axis="y",style="sci",scilimits=(0,0),useOffset=False)
axes[0].legend(fontsize=8);axes[-1].set_xlabel("Time (a.u.)")
fig.suptitle("Same physical parameters | 100,000 paths per curve | 30 orbitals / 15 electrons\n"
             "Different seeds; visual comparison, not an absolute error measurement",y=.995)
fig.tight_layout(rect=(0,0,1,.95));fig.savefig(FIG/"long_old_vs_new.png");plt.close(fig)
print(json.dumps({k:v for k,v in report.items() if k!="windows"},indent=2))
