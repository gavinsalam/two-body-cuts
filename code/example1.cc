/// example1.cc:
/// 
/// Example program to illustrate usage of cuts, given the 
/// momenta from a two-body decay
/// 
/// Build and run it with:
///
///     make example1
///     ./example1
///

// This file is part of the two-body-cuts project.
// Copyright 2021, Gavin Salam and Emma Slade (University of Oxford)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <memory>
#include "TwoBodyCuts.hh"

using namespace std;
using namespace fastjet;
// for two-body cuts
using namespace tbc;


int main(int argc, char** argv){
  // start off with some massless photon momenta
  double px1 = 40.0, py1 = 12.0, pz1 = 45.0; 
  double E1 = sqrt(px1*px1+py1*py1+pz1*pz1);
  PseudoJet gamma1(px1,py1,pz1,E1);

  double px2 = -55.0, py2 = -15.0, pz2 = -31.8; 
  double E2 = sqrt(px2*px2+py2*py2+pz2*pz2);
  PseudoJet gamma2(px2,py2,pz2,E2);

  // Create a two-body system (here generically referred to as a "Boson", 
  // though it need not actually correspond to a Boson)
  Boson boson(gamma1,gamma2);
  cout << "Boson is\n" << boson << endl;

  // Then create a range of cuts, and apply them (here we pass the cut
  // thresholds as a fraction of the boson mass and below, when we apply
  // the cut, we pass the boson divided by its mass)
  vector<CutsBase*> cuts;
  cuts.push_back(new CutsNormal  (0.25, 0.35));
  cuts.push_back(new CutsProduct (0.25, 0.35));
  cuts.push_back(new CutsCBIPt   (0.25, 0.35));
  cuts.push_back(new CutsCBIPtRap(0.25, 0.35));

  // for each of the hardness cuts, add in the info about the rapidity cuts
  for (auto & cut: cuts) {
    cut->add_upper_absrap_cut(2.37);
    cut->add_absrap_exclusion(1.37, 1.52);
  }

  for (const auto & cut: cuts) {
    cout << "Applying " << cut->description() << endl;
    // remember to divide out the boson mass if the cut thresholds passed
    // above were dimensionless
    cout << "Boson passes cuts? " << cut->pass(boson/boson.m()) << endl << endl;
  }

  // we could save this with a 
  for (auto & cut: cuts) {
    delete cut;
  } 
}
