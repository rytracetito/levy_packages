# levy_packages

R packages for degradation models based on Lévy and related processes.

| Package | Content | Status |
|---|---|---|
| `levysim` | Exact simulation of one or several trajectories (gamma, Wiener, inverse Gaussian, variance-gamma, Meixner and inverse gamma processes, homogeneous or not), plots and data, inspection schedules | In progress |
| `levyfit` | Maximum likelihood estimation for the same processes, from one or several trajectories, and inspection schedules | In progress |
| Goodness of fit | Goodness-of-fit tests for Lévy degradation processes | Planned |

## Installation

```r
remotes::install_github("rytracetito/levy_packages", subdir = "levysim")
remotes::install_github("rytracetito/levy_packages", subdir = "levyfit")
```
