"""CAR norm identities for the collective Slater jump at large fillings.
This validates an electronic operator, not a full AHM Monte Carlo simulation.
"""
from pathlib import Path
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from experiments.slater_path_kernel import SlaterPath
rng=np.random.default_rng(115925)
start=time.perf_counter();cases=[]
for norb,nel in [(10,5),(30,15),(50,25),(50,17),(64,32),(50,1),(50,49)]:
    error=0.;largest=0.
    for impurity in [False,True]:
        nb=nel-int(impurity)
        for trial in range(32):
            z=rng.normal(size=(norb-1,nb))+1j*rng.normal(size=(norb-1,nb))
            q=np.linalg.qr(z,mode='reduced')[0]
            u=np.zeros((norb,nel),dtype=complex)
            if impurity:u[0,0]=1;u[1:,1:]=q
            else:u[1:]=q
            c=np.r_[0.,rng.normal(size=norb-1)];c/=np.linalg.norm(c)
            projected=float(np.linalg.norm(q.conj().T@c[1:])**2)
            expected=1-projected if impurity else projected
            source=SlaterPath(u,1+0j,impurity);jumped=source.jump(c)
            measured=float(abs(jumped.weight)**2)
            error=max(error,abs(measured-expected));largest=max(largest,measured)
            assert error<2e-13,(norb,nel,impurity,error)
            assert jumped.impurity != impurity
            if measured>1e-25:
                assert np.max(abs(jumped.orbitals.conj().T@jumped.orbitals-np.eye(nel)))<2e-13
                assert np.max(abs(jumped.orbitals[0]-(np.eye(1,nel)[0] if jumped.impurity else 0)))<2e-13
    cases.append({'orbitals':norb,'electrons':nel,'random_states':64,'max_norm_identity_error':error,'max_normalized_jump_norm_squared':largest})
report={'scope':'Electronic collective-jump operator only; includes non-half-filled and one-hole systems. No physical 4000-au variance claim.',
 'identity':'impurity occupied: ||H1 psi||^2/g^2 = 1 - ||Q_bath^dagger c_hat||^2; impurity empty: same norm = ||Q_bath^dagger c_hat||^2',
 'cases':cases,'wall_seconds':time.perf_counter()-start}
out=Path(__file__).resolve().parents[1]/'scan_rate_4000_20260907/collective_jump_norm_validation.json';out.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
