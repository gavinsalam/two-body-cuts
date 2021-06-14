# Two-body cuts

![](https://img.shields.io/badge/C%2B%2B-11-green)
<!-- [![](https://img.shields.io/badge/arXiv-2006.nnnnn-blue)](https://arxiv.org/abs/2006.nnnnn) -->


Code that implements a variety of cuts for selecting two-body kinematics
at hadron colliders, including the cuts from
[arXiv:2006.nnnnn](https://arxiv.org/abs/2006.nnnnn).

Various examples are provided:

- [code/example1.cc](code/example1.cc): shows simple usage, given the kinematics
  for each of two decay products

- [code/example2.cc](code/example2.cc): illustrates a Monte Carlo calculation of
  the acceptance for specific Higgs kinematics, with the full setup of
  the cuts whose acceptance shows the least dependence on the Higgs
  boson transverse momentum, CBI_HR cuts (including an additional raised
  hardness threshold at high rapidities).

See the start of each file for instructions on building and running the
executable.

## Dependencies

Usage of the code requires

- a C++11 compiler
- an installation of [FastJet](http://fastjet.fr)
