/// example3.cc:
///
/// Example program to illustrate a Monte Carlo evaluation of cut
/// efficiencies for a given choice of Higgs kinematics (i.e. Higgs pt,
/// rapidity and mass), using events with defiducialisation.
///
/// Build and run it with:
///
///     make example4 ./example4
///
/// It takes a few seconds to run and gives a progress indicator as it
/// does.

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
#include <random>
#include "TwoBodyCuts.hh"

using namespace std;
using namespace fastjet;
// for two-body cuts
using namespace tbc;

int main(int argc, char** argv) {

  // first prepare the cuts that we will examine
  double yband_lo = 1.37;
  double yband_hi = 1.52;
  double ymax = 2.37;

  // the same transverse momentum cut applied will be applied to both
  // decay photons (or equivalently it's a cut on the soft photon pt)
  CutsDefid cuts_defid(0.25);
  // add in the info about the rapidity cuts
  cuts_defid.add_upper_absrap_cut(ymax);
  cuts_defid.add_absrap_exclusion(yband_lo, yband_hi);

  // choose some Higgs kinematics at which we wish to examine the
  // weighted acceptance
  double ptH = 25.0;
  double yH  = 0.4;
  double mH  = 125.0;
  cout << "Boson kinematics: ptH = " << ptH 
       << ", yH = " << yH
       << ", mH = " << mH << endl;

  // choose a random number generator and a number of events.
  //
  // NB: if you want to loop over different Higgs kinematic
  //     points in order to examine the ptH dependence, it is
  //     a good idea to use the same starting seed for each
  //     ptH, so that statistical errors are correlated 
  //     across the pTH points.
  unsigned nev = 1000000;
  int seed = 5;
  mt19937 generator(seed);
  auto phi_dist      = uniform_real_distribution<double>(-pi,pi);
  auto costheta_dist = uniform_real_distribution<double>(-1.0,1.0);
  
  // do a simple MC integration to evaluate the acceptance
  // for each of the cuts (using the a common set of boson
  // decays across all the cuts)
  double defiducialised_weight_no_Born = 0.0;
  double defiducialised_weight_Born = 0.0;
  for (unsigned iev = 0; iev < nev; iev++) {
    if (iev > 0 && iev % 100000 == 0) cout << "completed " << iev << " events out of " << nev << endl;
    // choose uniform random orientation for the Higgs decay
    double phi = phi_dist(generator);
    double theta = acos(costheta_dist(generator));
    // create boson (which by default is boosted along the x axis) and
    // evaluate whether it passes each of the cuts
    Boson boson(mH,theta,phi,ptH,yH);
    bool accept = cuts_defid.pass(boson/boson.m());

    if (accept) {
      // decide whether we want to defiducialise such that the weighted
      // sum over accepted events is one (false), or is equal to the 
      // Born acceptance at that rapidity (we will try both options here)
      bool normalise_to_Born = false;
      double acceptance = cuts_defid.scalar_acceptance_four_phi(boson/boson.m(), normalise_to_Born);
      defiducialised_weight_no_Born += 1.0 / acceptance;

      normalise_to_Born = true;
      acceptance = cuts_defid.scalar_acceptance_four_phi(boson/boson.m(), normalise_to_Born);
      defiducialised_weight_Born += 1.0 / acceptance;

    }
  }

  // evaluate a reference result for the expected Born acceptance
  Boson born_boson(mH,0.0,0.0,0.0,yH);
  double expected_born_acceptance = cuts_defid.scalar_acceptance_Born(born_boson/born_boson.m());

  cout << "\nResults for \n" << cuts_defid.description() << "\n"
        << "average defiducialised weight (should be equal to 1) = " << defiducialised_weight_no_Born /nev << "\n"
        << "average Born-defiducialised weight (should be equal to Born acceptance of "
        << expected_born_acceptance << ") = " << defiducialised_weight_Born /nev << "\n";

}
