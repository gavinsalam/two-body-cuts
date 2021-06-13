/// example2.cc:
///
/// Example program to illustrate a Monte Carlo evaluation of cut
/// efficiencies for a given choice of Higgs kinematics (i.e. Higgs pt,
/// rapidity and mass), as plotted throughout the associated paper.
///
/// Build and run it with:
///
///     make example2 ./example2
///
/// It takes a few seconds to run and gives a progress indicator as it
/// does.
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

  vector<CutsBase*> cuts;
  cuts.push_back(new CutsNormal  (0.25, 0.35));
  cuts.push_back(new CutsProduct (0.25, 0.35));

  // set up a default CBI_HR
  auto cbiHR = new CutsCBIPtRap(0.25, 0.35);
  // add in the higher ptCS cut for high rapidities
  double ytransition = (yband_lo + ymax)/2.0;
  double pthard_higher = 1.0 / 2.0 / cosh(yband_hi - ytransition);
  cbiHR->add_higher_pthard(ytransition, pthard_higher);
  cuts.push_back(cbiHR);

  // for each of the hardness cuts, add in the info about the rapidity cuts
  for (auto & cut: cuts) {
    cut->add_upper_absrap_cut(ymax);
    cut->add_absrap_exclusion(yband_lo, yband_hi);
  }
    
  // choose some Higgs kinematics at which we wish to examine the acceptace
  double ptH = 25.0;
  double yH  = 0.4;
  double mH  = 125.0;
  cout << "Boson kinematics: ptH = " << ptH 
       << ", yH = " << yH
       << ", mH = " << mH << endl;

  // choose a random number generator and a number of events
  unsigned nev = 1000000;
  int seed = 5;
  mt19937 generator(seed);
  auto phi_dist      = uniform_real_distribution<double>(-pi,pi);
  auto costheta_dist = uniform_real_distribution<double>(-1.0,1.0);
  
  // do a simple MC integration to evaluate the acceptance
  // for each of the cuts (using the a common set of boson
  // decays across all the cuts)
  vector<double> npass(cuts.size(),0);
  for (unsigned iev = 0; iev < nev; iev++) {
    if (iev > 0 && iev % 100000 == 0) cout << "completed " << iev << " events out of " << nev << endl;
    // choose uniform random orientation for the Higgs decay
    double phi = phi_dist(generator);
    double theta = acos(costheta_dist(generator));
    // create boson and evaluate whether it passes each of the cuts
    Boson boson(mH,theta,phi,ptH,yH);
    for (unsigned icut = 0; icut < cuts.size(); icut++) {
      if (cuts[icut]->pass(boson/boson.m())) npass[icut] += 1.0;
    }
  }

  // get the results
  for (unsigned icut = 0; icut < cuts.size(); icut++) {
    cout << "\nResults for \n" << cuts[icut]->description() << "\n"
         << "efficiency = " << npass[icut]/nev << "\n";
  }

}