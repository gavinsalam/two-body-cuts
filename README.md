# Two-body cuts

![](https://img.shields.io/badge/C%2B%2B-11-green)
[![](https://img.shields.io/badge/arXiv:-2006.nnnnn-blue)](https://arxiv.org/abs/2006.nnnnn)


Code that implements a variety of cuts for selecting two-body kinematics
at hadron colliders, including those from [arXiv:2006.nnnnn](https://arxiv.org/abs/2006.nnnnn).

Various examples are provided:

- [example1.cc](example1.cc): simple usage given kinematics for each of
  two decay products

- [example2.cc](example2.cc): illustrates a Monte Carlo calculation of
  the acceptance for given Higgs kinematics, with the full setup of the
  CBI_HR cuts (including an additional raised threshold at high
  rapidities).

## Dependencies

To compile, the code requires

- a C++11 compiler
- an installation of fastjet
