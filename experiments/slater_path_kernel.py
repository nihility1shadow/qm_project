"""Prototype ONLY: exact conditional electronic Slater paths, no nuclear sampler.

General nonorthogonal matrix element context: Burton, arXiv:2101.10944.
A fixed impurity-occupation star jump maps a determinant to a determinant.
This module does not claim a completed 4000-au Monte Carlo solver.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import json,time
import numpy as np

@dataclass
class SlaterPath:
    orbitals: np.ndarray
    weight: complex
    impurity: bool

    @classmethod
    def basis(cls,norb,occupied):
        occupied=sorted(occupied)
        return cls(np.eye(norb,dtype=complex)[:,occupied],1+0j,0 in occupied)

    def advance(self,energies,interval):
        u=self.orbitals*np.exp(-1j*np.asarray(energies)*interval)[:,None]
        weight=self.weight
        if self.impurity:
            weight*=u[0,0];u[0,0]=1.0
        return SlaterPath(u,weight,self.impurity)

    def jump(self,coupling):
        c=np.array(coupling,dtype=float,copy=True);c[0]=0
        u=self.orbitals;rows,nel=u.shape
        out=np.zeros_like(u)
        if self.impurity:
            # c_0 annihilation has no sign in the canonical first column.
            block=np.column_stack((c[1:],u[1:,1:]))
            q,r=np.linalg.qr(block,mode='reduced');out[1:]=q
            weight=self.weight*np.prod(np.diag(r))
        else:
            alpha=c@u;pivot=int(np.argmax(abs(alpha)));a=alpha[pivot]
            if a==0:return SlaterPath(out,0j,True)
            keep=[j for j in range(nel) if j!=pivot]
            block=u[1:,keep]-np.outer(u[1:,pivot],alpha[keep]/a)
            q,r=np.linalg.qr(block,mode='reduced');out[0,0]=1;out[1:,1:]=q
            weight=self.weight*((-1)**pivot)*a*np.prod(np.diag(r))
        return SlaterPath(out,weight,not self.impurity)

    def determinant_coefficients(self,basis):
        # Exhaustive expansion is for small-system validation ONLY.
        return self.weight*np.linalg.det(self.orbitals[np.asarray(basis)])

def overlap_and_occupations(bra,ket):
    b,a=bra.orbitals,ket.orbitals
    overlap_matrix=b.conj().T@a
    left,singular,right_h=np.linalg.svd(overlap_matrix)
    phase=np.linalg.det(left)*np.linalg.det(right_h)
    excluded=np.array([np.prod(np.delete(singular,k)) for k in range(len(singular))])
    # Adjugate formulation remains defined when an overlap singular value is zero.
    adjugate=phase*(right_h.conj().T*excluded)@left.conj().T
    weight=bra.weight.conjugate()*ket.weight
    overlap=weight*phase*np.prod(singular)
    occ=weight*np.einsum('aj,jk,ak->a',a,adjugate,b.conj())
    return overlap,occ

def exact_hopping_branches(path,coupling,interval):
    """P_a exp(-i H1 dt) P_b branches; both remain Slater determinants.

    Prototype uses the small-angle regime relevant to this project. An occupied
    stay branch needs a separate limiting formula when cos(g*dt) is near zero.
    """
    c=np.asarray(coupling,dtype=float).copy();c[0]=0;g=np.linalg.norm(c)
    if g==0:return path,SlaterPath(path.orbitals.copy(),0j,not path.impurity)
    cosine=np.cos(g*interval)
    if abs(cosine)<.1:raise ValueError('prototype stay-branch factorization requires |cos(g dt)| >= 0.1')
    bright=c[1:]/g;u=path.orbitals.copy();weight=path.weight
    bath=u[1:,1:] if path.impurity else u[1:,:]
    factor=1/cosine-1 if path.impurity else cosine-1
    bath=bath+factor*np.outer(bright,bright@bath)
    q,r=np.linalg.qr(bath,mode='reduced')
    if path.impurity:u[1:,1:]=q;weight*=cosine
    else:u[1:,:]=q
    weight*=np.prod(np.diag(r))
    stay=SlaterPath(u,weight,path.impurity)
    leave=path.jump(c/g);leave.weight*=-1j*np.sin(g*interval)
    return stay,leave

def oracle_hopping(norb,basis,coupling):
    index={tuple(s):i for i,s in enumerate(basis)}
    h=np.zeros((len(basis),len(basis)),dtype=complex)
    for source,occupied in enumerate(basis):
        for orbital in range(1,norb):
            if (0 in occupied)==(orbital in occupied):continue
            annihilate,create=(0,orbital) if 0 in occupied else (orbital,0)
            target=list(occupied);p=target.index(annihilate);target.pop(p)
            q=sum(x<create for x in target);target.insert(q,create)
            h[index[tuple(target)],source]+=((-1)**(p+q))*coupling[orbital]
    return h

def validate():
    rng=np.random.default_rng(114)
    max_state_relative=max_overlap_absolute=max_density_absolute=0.0
    state_checks=matrix_checks=graphs=kick_checks=0
    max_kick_error=max_spectral_error=0.0
    for norb in range(3,11):
        for nel in range(1,norb):
            basis=list(combinations(range(norb),nel));occ_matrix=np.zeros((len(basis),norb))
            for k,state in enumerate(basis):occ_matrix[k,list(state)]=1
            energies=rng.uniform(-.7,.2,norb)
            coupling=rng.uniform(-1,1,norb);coupling[0]=0;coupling/=np.linalg.norm(coupling)
            h=oracle_hopping(norb,basis,coupling)
            eigenvalues,eigenvectors=np.linalg.eigh(h)
            max_spectral_error=max(max_spectral_error,float(abs(np.max(abs(eigenvalues))-np.linalg.norm(coupling))))
            kick=(eigenvectors*np.exp(-1j*.23*eigenvalues))@eigenvectors.conj().T
            eocc=np.array([sum(energies[list(state)]) for state in basis])
            graphs+=1
            for initial in [tuple(range(nel)),tuple(range(norb-nel,norb))]:
                path=SlaterPath.basis(norb,initial)
                exact=np.zeros(len(basis),dtype=complex);exact[basis.index(initial)]=1
                saved=[]
                for hop in range(10):
                    interval=rng.uniform(0,6)
                    path=path.advance(energies,interval);exact*=np.exp(-1j*eocc*interval)
                    reconstructed=path.determinant_coefficients(basis)
                    relative=np.linalg.norm(reconstructed-exact)/max(np.linalg.norm(exact),1e-300)
                    max_state_relative=max(max_state_relative,float(relative));state_checks+=1
                    assert relative<3e-12,(norb,nel,hop,relative)
                    saved.append((path,exact.copy()))
                    path=path.jump(coupling);exact=h@exact
                    reconstructed=path.determinant_coefficients(basis)
                    relative=np.linalg.norm(reconstructed-exact)/max(np.linalg.norm(exact),1e-300)
                    max_state_relative=max(max_state_relative,float(relative));state_checks+=1
                    assert relative<3e-12,(norb,nel,hop,relative)
                for i in range(len(saved)):
                    state,coefficients=saved[i]
                    stay,leave=exact_hopping_branches(state,coupling,.23)
                    reconstructed=stay.determinant_coefficients(basis)+leave.determinant_coefficients(basis)
                    error=np.linalg.norm(reconstructed-kick@coefficients)/max(np.linalg.norm(coefficients),1e-300)
                    max_kick_error=max(max_kick_error,float(error));kick_checks+=1
                    assert error<3e-12
                    bra,bv=saved[i];ket,kv=saved[(i*3+1)%len(saved)]
                    overlap,occupations=overlap_and_occupations(bra,ket)
                    expected=np.conj(bv)*kv
                    eo=complex(expected.sum());en=expected@occ_matrix
                    max_overlap_absolute=max(max_overlap_absolute,float(abs(overlap-eo)))
                    max_density_absolute=max(max_density_absolute,float(abs(occupations-en).max()))
                    assert abs(overlap-eo)<3e-12 and np.max(abs(occupations-en))<3e-12
                    matrix_checks+=1
    # An exactly zero determinant overlap may have NONZERO one-body elements.
    b=np.zeros((4,2),complex);a=b.copy()
    b[:,0]=np.array([1,1,0,0])/np.sqrt(2);a[:,0]=np.array([1,-1,0,0])/np.sqrt(2)
    b[2,1]=a[2,1]=1
    overlap,occ=overlap_and_occupations(SlaterPath(b,1+0j,False),SlaterPath(a,1+0j,False))
    assert abs(overlap)<1e-14 and np.max(abs(occ-[.5,-.5,0,0]))<1e-14
    matrix_checks+=1
    # Large conditional states never allocate the combinatorial basis.
    large=SlaterPath.basis(30,range(15));c=np.r_[0.,np.ones(29)/np.sqrt(29)]
    energies=np.linspace(-.63,-.12,30)
    started=time.perf_counter()
    for _ in range(100):
        large=large.advance(energies,.71).jump(c)
        ol,no=overlap_and_occupations(large,large)
        assert abs(no.sum()-15*ol)<1e-12
    elapsed=time.perf_counter()-started
    return {'status':'conditional-electronic kernel only; nuclear Poisson estimator not integrated',
            'small_hilbert_graphs':graphs,'state_checks':state_checks,'matrix_element_checks':matrix_checks,
            'exact_hopping_branch_checks':kick_checks,'max_exact_kick_relative_error':max_kick_error,
            'max_hopping_spectral_norm_error':max_spectral_error,
            'max_state_relative_error':max_state_relative,'max_overlap_absolute_error':max_overlap_absolute,
            'max_occupation_absolute_error':max_density_absolute,'singular_overlap_nonzero_density_passed':True,
            'large_case':{'orbitals':30,'electrons':15,'conditional_orbital_matrix_bytes':large.orbitals.nbytes,
                          '100_jump_and_observable_checks_seconds_python':elapsed,
                          'full_basis_allocated':False},
            'limits':'No demonstrated full trajectory variance improvement or 4000-au fitting score. Reference-subspace control and Poisson time/count weights must be rederived and validated before integration.',
            'source_context':'https://arxiv.org/abs/2101.10944'}

if __name__=='__main__':
    result=validate()
    destination=Path(__file__).resolve().parents[1]/'scan_v114_slater_kernel_20260906'
    destination.mkdir(exist_ok=True)
    (destination/'kernel_validation.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
