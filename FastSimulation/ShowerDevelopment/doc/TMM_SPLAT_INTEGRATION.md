# TMM splatting in FastSim: table-driven integration design

The Track-B learned-mixture shower model ("TMM splatting"), integrated
without any ML runtime: every learned object is exported to plain
tables; the C++ side is a sampler plus closed-form rendering.

## The model, as the C++ sees it

A shower at incident energy E (x = ln E) is:

1. K ~ P(K | x)          -- multinomial over component counts,
                            tabulated on an x grid.
2. theta ~ copula_K(x)   -- the (7K+1)-dim deck
                            [ln pi | mu_x | mu_y | mu_z | ln sigma_r |
                             ln sigma_z | ln nu] per component + ln E_w:
                            mean_d(x)  = quadratic polynomial (3 coeffs/dim),
                            sigma_d(x), corr(x) = linear interpolation
                            between per-bin tables (16 bins),
                            drawn via Cholesky of the interpolated
                            correlation; the last dim maps through an
                            empirical quantile table (65 points) instead
                            of the Gaussian marginal.
3. render: for each layer l, component k:
     zmass = Phi((zhi_l - mu_z)/sigma_z) - Phi((zlo_l - mu_z)/sigma_z)
     kernel(r) = bivariate Student-t(nu_k, sigma_r,k)   (Gaussian when
                 nu > 1e5); energy density = E_w * sum_k pi_k * zmass *
                 kernel, evaluated at cell centres.
4. noise: cell energy = e_spot(l) * Poisson(density * A_cell / e_spot(l)).
5. calibrations (all precomputed tables vs x): per-layer e_spot,
   per-ring factors, per-layer centroid offsets. Exact rescale of the
   weighted total to E_w closes energy conservation.

Positions are in the incidence frame (nominal shower axis); the caller
transplants to the track's impact point and direction, as GFlash does.

## Files

- `interface/HGCalTMMShower.h`, `src/HGCalTMMShower.cc` -- pure-C++
  sampler + renderer (stdlib only; RNG injected as uniform/gauss
  functors so the CMSSW RandomEngineAndDistribution adapter is a
  10-line wrapper).
- `test/export_tmm_params.py` -- exports the fitted archive
  (theta populations per exact K, P(K|x), calibration grids) to the
  flat text format `data/tmm_splat_photon.params` (documented in the
  exporter header).
- Wiring point: CalorimetryManager's HGCal EM branch, behind a config
  flag `simulateTMMSplat` (default False) -- to be added once the
  standalone class is build-validated.

## Status

Skeleton branch: class + exporter + this note. Parameter file and the
CalorimetryManager wiring land after build validation on lxplus.
Source analysis lives in the fastsim_eval workspace (RESULTS.md).
