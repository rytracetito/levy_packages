# levy_packages

R packages for degradation models based on Lévy and non-homogeneous Lévy processes.

| Package | Content | Status |
|---|---|---|
| `DegradSim` | Exact simulation of one or several trajectories of Lévy processes (gamma, Wiener, inverse Gaussian, variance-gamma, Meixner, normal inverse Gaussian, compound Poisson) and of non-homogeneous Lévy processes X_t = Y_{A(t)} (gamma, Wiener, inverse Gaussian), plots and data, inspection schedules | In progress |
| `levyfit` | Maximum likelihood estimation for the five Lévy processes, from one or several trajectories | In progress |
| `levyBridge` | Exact simulation of the bridges of the gamma, Wiener, inverse Gaussian, variance-gamma and Meixner processes, and resampling of trajectories observed at irregular times onto a regular grid (Grall-Maës, 2012); reconstructed increments for the tests | In progress |
| `levygof` | Goodness-of-fit tests for the gamma, Wiener, inverse Gaussian, variance-gamma and Meixner processes, from one or several trajectories inspected at regular or irregular times (estimation with `levyfit`, bridges with `levyBridge`, methods "Bridge", "ExtBridge" and "MultiBridge") | In progress |

## Installation

```r
remotes::install_github("rytracetito/levy_packages", subdir = "DegradSim")
remotes::install_github("rytracetito/levy_packages", subdir = "levyfit")
remotes::install_github("rytracetito/levy_packages", subdir = "levyBridge")
remotes::install_github("rytracetito/levy_packages", subdir = "levygof")
```
