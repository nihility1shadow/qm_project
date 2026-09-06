"""Independent coherent-state/Poisson weight prototype, not a production solver.

The continuous-time propagator tested here differs from the production grid
split propagator. Kernel agreement does not demonstrate Monte Carlo accuracy.
"""
from dataclasses import dataclass
from pathlib import Path
import json,math,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import gammaln
from slater_path_kernel import SlaterPath,overlap_and_occupations

@dataclass
class CoherentSlaterPath:
    electronic: SlaterPath
    alpha: complex

    def advance(self,energies,frequency,displacement,interval):
        center=displacement if self.electronic.impurity else 0.0
        alpha=(self.alpha-center)*np.exp(-1j*frequency*interval)+center
        phase=np.exp(-.5j*frequency*interval+1j*center*(self.alpha.imag-alpha.imag))
        electronic=self.electronic.advance(energies,interval)
        electronic.weight*=phase
        return CoherentSlaterPath(electronic,alpha)

    def jump(self,coupling,multiplier=1.0):
        electronic=self.electronic.jump(coupling)
        electronic.weight*=multiplier
        return CoherentSlaterPath(electronic,self.alpha)


def path_at_time(norb,occupied,energies,coupling,frequency,displacement,alpha,
                 horizon,jump_times,rate):
    if rate<=0 or horizon<0:raise ValueError('positive rate and nonnegative horizon required')
    times=np.asarray(jump_times,dtype=float)
    if not np.isfinite(times).all() or np.any(np.diff(times)<0) or np.any(times<0) or np.any(times>horizon):
        raise ValueError('jump times must be ordered inside the horizon')
    path=CoherentSlaterPath(SlaterPath.basis(norb,occupied),complex(alpha))
    path.electronic.weight*=np.exp(rate*horizon)
    previous=0.0
    for jump in times:
        path=path.advance(energies,frequency,displacement,float(jump)-previous)
        path=path.jump(coupling,-1j/rate);previous=float(jump)
    return path.advance(energies,frequency,displacement,horizon-previous)


def pair_observables(bra,ket):
    overlap,occupations=overlap_and_occupations(bra.electronic,ket.electronic)
    # This form computes overlap modulus without cancellation of large norms.
    nuclear=np.exp(-.5*abs(bra.alpha-ket.alpha)**2+1j*(bra.alpha.conjugate()*ket.alpha).imag)
    return nuclear*overlap,nuclear*occupations


def coherent_coefficients(alpha,nfock):
    if alpha==0:
        values=np.zeros(nfock,complex);values[0]=1;return values
    n=np.arange(nfock)
    return np.exp(-.5*abs(alpha)**2+n*np.log(abs(alpha))-.5*gammaln(n+1)+1j*n*np.angle(alpha))


def validate():
    frequency=3.6749323758566211e-3;mass=14583.1067146087
    displacement=2*np.sqrt(mass*frequency/2)
    rng=np.random.default_rng(11421)
    histories=[[],[500.],[500.,1000.,2000.,3000.],list(np.sort(rng.uniform(0,4000,8)))]
    # Unshifted oscillator basis; this dense oracle is used ONLY in this test.
    maximum_alpha=abs(displacement)
    for history in histories:
        alpha=complex(displacement);previous=0.;occupied=True
        for endpoint in history+[4000.]:
            center=displacement if occupied else 0.
            maximum_alpha=max(maximum_alpha,abs(center)+abs(alpha-center))
            alpha=(alpha-center)*np.exp(-1j*frequency*(endpoint-previous))+center
            maximum_alpha=max(maximum_alpha,abs(alpha));previous=endpoint;occupied=not occupied
    nfock=int(np.ceil(maximum_alpha**2+12*maximum_alpha+100));n=np.arange(nfock)
    eig=[]
    for occupied in [False,True]:
        center=displacement if occupied else 0
        eig.append(eigh_tridiagonal(frequency*(n+.5+center**2),
                                   -frequency*center*np.sqrt(np.arange(1,nfock))))
    errors=[];norm_errors=[];saved=[]
    for history in histories:
        path=CoherentSlaterPath(SlaterPath.basis(3,[0]),complex(displacement))
        exact=coherent_coefficients(displacement,nfock);previous=0.0
        for endpoint in history+[4000.]:
            elapsed=endpoint-previous
            values,vectors=eig[int(path.electronic.impurity)]
            exact=vectors@((vectors.T@exact)*np.exp(-1j*values*elapsed))
            path=path.advance(np.zeros(3),frequency,displacement,elapsed)
            reconstructed=path.electronic.weight*coherent_coefficients(path.alpha,nfock)
            error=float(np.linalg.norm(reconstructed-exact));errors.append(error)
            norm_errors.append(float(abs(np.vdot(exact,exact)-1)))
            assert error<3e-10,(history,endpoint,error)
            previous=endpoint
            if endpoint!=4000.:
                # Toggle only the nuclear Hamiltonian for this nuclear oracle.
                path.electronic.impurity=not path.electronic.impurity
        saved.append((path,exact.copy()))
    nuclear_overlap_error=0.0
    for bra,bv in saved:
        for ket,kv in saved:
            actual=(bra.electronic.weight.conjugate()*ket.electronic.weight*
                    np.exp(-.5*abs(bra.alpha-ket.alpha)**2+1j*(bra.alpha.conjugate()*ket.alpha).imag))
            nuclear_overlap_error=max(nuclear_overlap_error,float(abs(actual-np.vdot(bv,kv))))
    assert nuclear_overlap_error<3e-10

    # Exact Poisson-count quadrature of H1 alone: there is no time sampling
    # noise because H0=0. This is a weighting/fermion test, NOT the AH model.
    results=[]
    for norb,nel in [(4,1),(6,3),(10,5),(30,15)]:
        horizon=4000.;rate=np.sqrt(3e-7);coupling=np.r_[0.,np.full(norb-1,np.sqrt(3e-7/(norb-1)))]
        h=np.zeros((norb,norb));h[0]=h[:,0]=coupling
        val,vec=np.linalg.eigh(h);u=(vec*np.exp(-1j*val*horizon))@vec.conj().T
        expected=np.sum(abs(u[:,:nel])**2,axis=1)
        weighted=[]
        for k in range(25):
            times=np.linspace(0,horizon,k+2)[1:-1]
            path=path_at_time(norb,range(nel),np.zeros(norb),coupling,0.,0.,0.,horizon,times,rate)
            probability=np.exp(-rate*horizon)*(rate*horizon)**k/math.factorial(k)
            path.electronic.weight*=probability;weighted.append(path)
        overlap=0j;occupations=np.zeros(norb,complex)
        for bra in weighted:
            for ket in weighted:
                ol,no=pair_observables(bra,ket);overlap+=ol;occupations+=no
        error=float(np.max(abs(occupations-expected)));norm_error=float(abs(overlap-1))
        assert error<2e-12 and norm_error<2e-12,(norb,nel,error,norm_error)
        results.append({'orbitals':norb,'electrons':nel,'horizon_au':horizon,'maximum_count':24,
                        'max_occupation_error':error,'norm_error':norm_error,
                        'full_electronic_basis_allocated':False})
    return {'scope':'Exact conditional nuclear propagation and Poisson weighting only; no sampled AH accuracy claim',
            'nuclear_horizon_au':4000,'nuclear_histories':len(histories),'nuclear_interval_checks':len(errors),
            'nuclear_fock_oracle_levels':nfock,'maximum_conditional_alpha':maximum_alpha,'max_nuclear_state_error':max(errors),
            'max_nuclear_oracle_norm_error':max(norm_errors),'max_nuclear_overlap_error':nuclear_overlap_error,
            'poisson_count_checks_with_H0_disabled':results,
            'limits':['The count sum is deterministic; it does not measure stochastic variance.',
                      'The production method splits nuclear T/V, which is not identical to continuous H0.',
                      'Reference projection generally destroys a single Slater representation.',
                      'No claim of improved 4000-au AH fitting follows from these kernel checks.']}

if __name__=='__main__':
    result=validate();dest=Path(__file__).resolve().parents[1]/'scan_v114_slater_kernel_20260906'
    dest.mkdir(exist_ok=True);(dest/'coherent_poisson_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
