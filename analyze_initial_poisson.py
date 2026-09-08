"""Source and archival-data audit of initial Poisson variants; NOT a matched runtime benchmark."""
from pathlib import Path
from math import comb
import hashlib,json,subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_versions_4000 import bath
ROOT=Path(__file__).resolve().parent;OUT=ROOT/'scan_initial_poisson_audit_20260908'


def main():
    specs={'sources/ahm-cs-mb-v0.1.0.cpp':'v0.1.0:ahm-cs-mb.cpp','sources/ahm-mb-sep-before-v030.cpp':'v0.10.0:archive/ahm-mb-sep-before-v030.cpp','sources/yyy-inlines-v0.1.0.h':'v0.1.0:yyy_inlines.h','sources/Johnson-path-sampler-v0.1.0.cpp':'v0.1.0:Johnson-path-sampler.cpp','initial-import-output.dat':'v0.1.0:text 1000 10 5 (0).dat','initial-import-qm.dat':'v0.1.0:ahm-qm-s10-n5.dat'}
    provenance=[]
    for name,spec in specs.items():
        original=subprocess.check_output(['git','show',spec]);local=(OUT/name).read_bytes()
        # Git text checkout can change LF to CRLF; content must otherwise remain exact.
        if local.replace(b'\r\n',b'\n')!=original.replace(b'\r\n',b'\n'):raise ValueError('Snapshot contents differ')
        provenance.append({'path':name,'git_spec':spec,'git_blob':subprocess.check_output(['git','rev-parse',spec]).decode().strip(),'original_sha256':hashlib.sha256(original).hexdigest()})
    source=(OUT/'sources/ahm-cs-mb-v0.1.0.cpp').read_text();sep=(OUT/'sources/ahm-mb-sep-before-v030.cpp').read_text()
    for marker in ['#define _YYY_REALSPACE_AV_','nwf   = nstep/5','array3d<dcomplex>(nwf+1, Nhs, npt)','array3d<dcomplex>(nwf+1, Nhs, nmem)']:
        if marker not in source:raise ValueError('Initial MB storage formula changed')
    for marker in ['array2d<dcomplex>(nstep+1, Norb)','nbloc = nstep/nwf','prb[j/nbloc]','t*dt, prb[t][0]']:
        if marker not in sep:raise ValueError('Original SepMB audit markers changed')
    a=np.loadtxt(OUT/'initial-import-output.dat');q=np.loadtxt(OUT/'initial-import-qm.dat')
    assert a.shape==(201,14) and q.shape==(8001,14) and np.isfinite(a).all() and np.isfinite(q).all()
    assert np.array_equal(bath(OUT/'initial-import-output.dat'),bath(OUT/'initial-import-qm.dat'))
    current=ROOT/'scan_v113_two_buffer_20260905/647656/ahm-qm-s10-n5.dat'
    header_old=bath(OUT/'initial-import-output.dat');header_current=bath(current)
    current_matches=np.array_equal(header_old,header_current)
    sums=a[:,4:].sum(axis=1);qsum=q[:,4:].sum(axis=1)
    rows=[]
    for norb,nel in [(10,5),(30,15),(50,25)]:
        one=16*1024*comb(norb,nel);mb=1601*one
        rows.append({'orbitals':norb,'electrons':nel,'determinants':comb(norb,nel),'grid_QM_one_wavefunction_bytes':one,'initial_MB_wf_array_bytes':mb,'initial_MB_MPI_wf_plus_avg_arrays_bytes_per_rank':2*mb,'initial_MB_wf_array_GiB':mb/1024**3})
    report={'scope':'Initial-source and stored-file audit; no original executable was timed or memory-profiled.',
      'initial_tag':'v0.1.0','initial_commit':subprocess.check_output(['git','rev-parse','v0.1.0^{commit}']).decode().strip(),
      'initial_import_date':'2026-06-14','provenance':provenance,
      'variant_distinction':'Initial v0.1.0 implements MBpoisson storing time-resolved full-basis wavefunctions. The separately preserved original SepMBpoisson directly estimates observables and does not allocate that wavefunction array. Neither should be identified as v0.93.',
      'initial_snapshot_incomplete_for_standalone_build':'Initial tree does not contain a main program, build recipe or implementations for all declared model/nuclear routines; Johnson .cpp is a class declaration placeholder and yyy_inlines.h has an invalid two-argument malloc call.',
      'archived_output':{'shape':list(a.shape),'printed_time_range':[float(a[0,0]),float(a[-1,0])],'initial_particle_sum':float(sums[0]),'first_nonzero_row_particle_sum':float(sums[1]),'last_row_particle_sum':float(sums[-1]),'minimum_particle_sum':float(sums.min()),'negative_particle_sum_rows':int((sums<0).sum()),'maximum_particle_sum_deviation_from_5':float(abs(sums-5).max()),'QM_max_particle_sum_deviation_from_5':float(abs(qsum-5).max()),
      'matches_archived_QM_bath_header':True,'matches_current_QM_bath_header':bool(current_matches),'initial_impurity_energy_Ha':float(header_old[0,2]),'current_impurity_energy_Ha':float(header_current[0,2]),'initial_bath_coupling_Ha':float(header_old[1,1]),'current_bath_coupling_Ha':float(header_current[1,1]),
      'warning':'Output filename/build provenance and actual physical-time mapping are not established. Original SepMB stores prb[j/nbloc] but prints t*dt. Do not relabel this file to 4000 au or compute a trustworthy matched-Q score from row alignment. Archived QM is also only historical data, not revalidated truth.'},
      'capacity_assumptions':{'nstep':8000,'dt':.5,'time_range_au':[0,4000],'nuclear_grid_points':1024,'complex_bytes':16,'stored_time_slices':1601,'notice':'Component-size estimates under unified dimensions, not measured RSS or successfully allocated memory. Full-basis integer/index limits can fail before allocation. These arrays belong to MBpoisson, NOT original SepMBpoisson.'},'capacity_estimates':rows,
      'original_runtime_seconds':None,'original_process_peak_rss_bytes':None,'original_matched_4000_Q':None}
    (OUT/'audit.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    fig,axes=plt.subplots(1,2,figsize=(12.6,4.9));row=np.arange(len(a))
    axes[0].plot(row,a[:,4],label='Imported Poisson output, orbital 0',color='#aa6e48',lw=1)
    axes[0].plot(row,q[:len(a),4],label='Imported QM file, orbital 0',color='#415c6a',ls='--',lw=1)
    axes[0].set(xlabel='Stored row index',ylabel='Raw stored occupation',title='Stored occupations show a normalization/estimation problem')
    axes[0].legend(fontsize=8)
    axes[1].semilogy(row,np.maximum(abs(sums),1e-300),color='#aa6e48',label='Absolute sum of stored occupations')
    axes[1].axhline(5,color='#415c6a',ls='--',label='Expected particle count: 5')
    axes[1].set(xlabel='Stored row index',ylabel='Absolute occupation sum',title='The imported output collapses in magnitude');axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.15)
    fig.suptitle('Initial v0.1.0 import: archived data audit, NOT a matched accuracy/speed benchmark\nRow index is shown because shared physical-time mapping and generating binary are unverified',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.88));fig.savefig(OUT/'figures/initial_output_audit.png',dpi=145);plt.close(fig)
    print(json.dumps({'initial_commit':report['initial_commit'],'archived_output':report['archived_output'],'capacity_estimates':rows},indent=2))

if __name__=='__main__':main()
