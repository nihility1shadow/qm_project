"""Reproducible figures and accuracy/resource checks for the v1.11 MPI reference."""
import json, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "scan_v111_distributed_20260905"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 160})

def read_run(job):
    folder = OUT / str(job)
    data = np.loadtxt(next(folder.glob("ahm-sepmb-*.dat")))
    log = (folder/"program.out").read_text()
    wall = (folder/"time.txt").read_text()
    def value(pattern):
        found = re.search(pattern, log)
        return float(found.group(1)) if found else None
    return data, {"job": job, "wall_seconds": float(re.search(r"wall_seconds=([\d.]+)",wall)[1]),
        "max_rank_peak_mib": value(r"max_rank_peak_rss_kb=(\d+)") / 1024,
        "sum_rank_peak_mib": (value(r"sum_rank_peak_rss_kb=(\d+)") or 0) / 1024 or None,
        "reference_seconds": value(r"reference_seconds=([\d.]+)"),
        "max_particle_error": float(np.max(np.abs(data[:,4:].sum(axis=1)-data[0,4:].sum()))),
        "finite": bool(np.isfinite(data).all())}

report = {"regressions": [], "resources": []}
for a,b in [(647609,647606),(647606,647607),(647610,647611)]:
    da,_=read_run(a); db,_=read_run(b)
    error=float(np.max(np.abs(da-db)))
    assert da.shape == db.shape and error < 2e-12, (a,b,error)
    report["regressions"].append({"jobs":[a,b], "max_abs_difference":error})
for job in [647609,647606,647607,647608]:
    _, resource = read_run(job)
    report["resources"].append(resource)

sim, _ = read_run(647607)
qm_path = ROOT/"scan_v093_poisson_t2500_20260822/pilots/qm_646419_eta03/ahm-qm-s10-n5.dat"
if not qm_path.exists():
    qm_path = OUT / "qm_10_5_t100.dat"
qm = np.loadtxt(qm_path)
qm = qm[qm[:,0] <= sim[-1,0]+1e-10]
assert np.array_equal(qm[:,0],sim[:,0])
np.savetxt(OUT/"qm_10_5_t100.dat", qm, fmt="%.16e",
           header="Same-parameter grid QM excerpt, original job 646419; 10 orbitals / 5 electrons.")
initial = np.r_[np.ones(5),np.zeros(5)]
active = [0,5,6,7,8,9]
signal = qm[:,4:]-initial
error = sim[:,4:]-qm[:,4:]
rms=lambda x:float(np.sqrt(np.mean(x*x)))
q = rms(signal[:,active])/rms(error[:,active])
report["small_system_qm"] = {"reference":str(qm_path),"active_orbitals":active,
    "tmax":float(sim[-1,0]),"Ntraj":1000,
    "Q_active":q,"min_active_orbital_Q":min(rms(signal[:,i])/rms(error[:,i]) for i in active),
    "error_rms_active":rms(error[:,active]),"max_abs_error_all":float(np.abs(error).max())}
fig,axes=plt.subplots(5,2,figsize=(11,12),sharex=True)
for i,ax in enumerate(axes.flat):
    ax.plot(qm[:,0],signal[:,i],color="#263238",lw=1.3,label="Grid QM")
    ax.plot(sim[:,0],sim[:,4+i]-initial[i],color="#e87b36",ls="--",lw=1.1,label="Poisson + MPI reference")
    ax.set_title(f"Orbital {i}")
    ax.ticklabel_format(axis="y",style="sci",scilimits=(0,0),useOffset=False)
    ax.grid(alpha=.18)
axes[0,0].legend(fontsize=8)
for ax in axes[-1]:ax.set_xlabel("Time (a.u.)")
fig.suptitle(f"10 orbitals / 5 electrons | 1,000 paths | 4 ranks | 100 a.u.\n"
             f"Same Hamiltonian; active-orbital Q = {q:.3g}; bounded reference = 246 / 252 states",y=.995)
fig.tight_layout(rect=(0,0,1,.965))
fig.savefig(FIG/"small_system_qm.png");plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(10,3.7))
labels=["v1.10","v1.11 MPI off","v1.11 MPI on"]
rr=report["resources"][:3]
for ax,field,title in zip(axes,["wall_seconds","max_rank_peak_mib"],
                          ["Total wall time (seconds)","Maximum rank peak memory (MiB)"]):
    vals=[r[field] for r in rr]
    bars=ax.bar(labels,vals,color=["#96a6b4","#749cb5","#338b80"],width=.6)
    for bar,val in zip(bars,vals):ax.text(bar.get_x()+bar.get_width()/2,val,f"{val:.2f}",ha="center",va="bottom")
    ax.set_ylim(0,max(vals)*1.22);ax.set_title(title);ax.grid(axis="y",alpha=.15)
fig.suptitle("Identical 10/5 test: 100 a.u., 1,000 paths, 4 ranks; output differences = 0")
fig.tight_layout();fig.savefig(FIG/"regression_resources.png");plt.close(fig)
(OUT/"stage1_validation.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
print(json.dumps(report,indent=2))
