# Two-body cuts

![](https://img.shields.io/badge/C%2B%2B-11-green)
[![](https://img.shields.io/badge/arXiv-2106.08329-blue)](https://arxiv.org/abs/2106.08329)


Code that implements a variety of cuts for selecting two-body kinematics
at hadron colliders, including the cuts from
[arXiv:2106.08329](https://arxiv.org/abs/2106.08329).

Various examples are provided:

- [code/example1.cc](code/example1.cc): shows simple usage of the cuts,
  given the kinematics for each of two decay products. The cuts
  themselves are available through the
  [code/TwoBodyCuts.hh](code/TwoBodyCuts.hh) file, and are in the `tbc`
  namespace.

- [code/example2.cc](code/example2.cc): illustrates a Monte Carlo
  calculation of the acceptance for specific Higgs kinematics, with the
  full setup of the cuts whose acceptance shows the least dependence on
  the Higgs boson transverse momentum, compensating boost-invariant
  (CBI<sub>HR</sub>) cuts from section 5 of the paper (including an additional
  raised hardness threshold at high rapidities).

- [code/example3.cc](code/example3.cc): a similar program that evaluates
  harmonic acceptances for Z production for a variety of cuts, 
  illustrating also the use of the CBI<sub>H,DY</sub> cuts shown in section 6 of the paper.

- [code/example4.cc](code/example4.cc): a program to illustrate the use of 
  (4-phi) defiducialisation for Higgs decays with symmetric cuts.

See the start of each file for instructions on building and running the
executables.



## Dependencies

Usage of the code requires

- a C++11 compiler
- an installation of [FastJet](http://fastjet.fr)
