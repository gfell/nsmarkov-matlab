#  nsmarkov-matlab

This repository provides **Matlab** functions to construct Markov chain approximations of (non-stationary) AR(1) processes as described in the paper "Markov-Chain Approximations for Life-Cycle Models"  by Giulio Fella, Giovanni Gallipoli and Jutong Pan, _Review of Economic Dynamics_ 34, 2019 ([https://doi.org/10.1016/j.red.2019.03.013](https://doi.org/10.1016/j.red.2019.03.013)). 

## Contains

- *lcrouwenhorst.m* subroutine to discretise a non-stationary AR(1) using our extension of Rouwenhorst [1995. "Asset pricing implications of equilibrium business cycle models," in `Frontiers of business cycle research', T. F. Cooley ed., Princeton University Press, Chapter 10.]
- *lctauchen.m* subroutine to discretise a non-stationary AR(1) using our extension of Tauchen [1986. "Finite State Markov-Chain Approximations to Univariate and
                  Vector Autoregressions," _Economics Letters_ 20].

## Installation

```
Pkg.clone("https://github.com/gfell/nsmarkov-matlab")
```
