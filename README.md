# levy_packages

R packages for degradation models based on Lévy and related processes.

| Package | Content | Status |
|---|---|---|
| `levysim` | Exact simulation of one or several trajectories (gamma, Wiener, inverse Gaussian, variance-gamma, Meixner and inverse gamma processes, homogeneous or not), plots and data, inspection schedules | In progress |
| `levyfit` | Maximum likelihood estimation for the same processes, from one or several trajectories, and inspection schedules | In progress |
| `levyBridge` | Exact simulation of the bridges of the gamma, Wiener, inverse Gaussian, variance-gamma and Meixner processes, and resampling of trajectories observed at irregular times onto a regular grid (Grall-Maës, 2012); reconstructed increments for the tests | In progress |
| `levygof` | Goodness-of-fit tests for the gamma, Wiener, inverse Gaussian, variance-gamma and Meixner processes, from one or several trajectories inspected at regular or irregular times (estimation with `levyfit`, bridges with `levyBridge`, methods "Bridge" and "ExtBridge") | In progress |

## Installation

```r
remotes::install_github("rytracetito/levy_packages", subdir = "levysim")
remotes::install_github("rytracetito/levy_packages", subdir = "levyfit")
remotes::install_github("rytracetito/levy_packages", subdir = "levyBridge")
remotes::install_github("rytracetito/levy_packages", subdir = "levygof")
```
