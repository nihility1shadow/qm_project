"""Finite-grid jump probability audit, matching the binary's strict u < p comparison."""
from pathlib import Path
from math import sqrt
import json


def midpoint_count(paths,probability):
    # Binary search counts the actual floating-point grid, including equality cases.
    lo,hi=0,paths
    while lo<hi:
        mid=(lo+hi)//2
        if (mid+.5)/paths<probability:lo=mid+1
        else:hi=mid
    return lo


def main():
    cases=[(10,5,10000),(10,5,1000000),(30,15,100000),(30,15,1000000),(50,25,1000),(50,17,1000),(50,25,10000),(50,25,1000000)]
    rows=[]
    for norb,nel,paths in cases:
        strength=1.5*.5*sqrt(nel*(norb-nel))*sqrt(3e-7/(norb-1))
        p=strength/(1+strength);count=midpoint_count(paths,p);effective=count/paths
        rows.append({'orbitals':norb,'electrons':nel,'paths':paths,'nominal_jump_probability':p,'actual_jumps_per_step':count,'midpoint_actual_probability':effective,'relative_probability_error':effective/p-1,'expected_jumps_in_4000_au':8000*p,'midpoint_actual_mean_jumps':8000*effective})
    # Check strict-boundary counting independently by enumerating small grids.
    checks=0
    for paths in [1,2,7,64,101]:
        for p in [0.,1.,.137,*[(k+.5)/paths for k in range(paths)]]:
            if midpoint_count(paths,p)!=sum((k+.5)/paths<p for k in range(paths)):raise ValueError('Count audit boundary failure')
            checks+=1
    report={'scope':'Exact marginal event-count audit of the current fixed-midpoint per-step stratification. It is NOT the physical occupation error, and hybrid control can cancel some biases.',
      'cause':'u=(permuted_index+0.5)/N; K=# {k: (k+0.5)/N < p}. Away from midpoint equality K=floor(N*p+0.5); the strict boundary is handled by binary search using the actual floating-point comparisons.',
      'available_alternative':'SEP_MB_STRATIFY_FORWARD_STEPS=0 with SEP_MB_STRATIFY_FORWARD_COUNT=1 uses randomized shifted CDF sampling already present in v1.14; benchmark before recommending.',
      'boundary_checks':checks,'rows':rows}
    Path('scan_rate_4000_20260907/finite_strata_probability_audit.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
