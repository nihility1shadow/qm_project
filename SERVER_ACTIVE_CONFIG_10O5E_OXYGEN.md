# Active server model configuration

The server model has been restored to the generic many-body executable with the
following active physical configuration:

- `Norb=10`
- `Nel=5`
- `delE=-0.6319 Ha=-17.1948749688 eV`
- oxygen orbital: atomic `2p` HF-Roothaan energy
- default `wc=10 eV`
- default `eta=0.01`
- Poisson implementation: full many-body Kondo/Johnson-state path sampler

The executable and source were not replaced because they already contain this
configuration as their defaults. The previous four-orbital experiments were
selected only by runtime arguments and environment variables.

Explicit cloud submission template:

```bash
sbatch --ntasks=64 \
  --export=ALL,AHM_WC_EV=10,AHM_ETA=0.01,AHM_DELE_EV=-17.194874968839816 \
  run_mpi_cloud.slurm NTRAJ 10 5 1
```

Current server hashes:

- `na_mpi_cloud.out`: `0feef1a464cf9bd22e76af7a64c46dc12c17f249d3c32a9d0840949fb18d6bbf`
- `ahm-mb-sep.cpp`: `843694ce9fbcca45f6cfb170163953b78a08a177b270353b329b52d08fef2849`

The last four-orbital result remains recoverable from Git commit `bf557a4`.

Rollback verification:

- Slurm smoke job: `636250`
- state: `COMPLETED`
- exit code: `0:0`
- generated file: `ahm-sepmb-s10-n5-100.dat`
- runtime header: `wc=10 eV`, `eta=0.01`,
  `delE=-17.194874968839816 eV`, `nstep=4`, `dt=0.5`
- smoke-test output was removed after verification; only
  `latest_ahm_config.txt` remains on the server as the active configuration
  marker
