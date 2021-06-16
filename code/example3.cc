/// example3.cc:
///
/// Example program to illustrate a Monte Carlo evaluation of cut
/// efficiencies for a given choice of Z boson kinematics (i.e. Z pt,
/// rapidity and mass), as plotted in section 6 of the accompanying
/// paper.
///
/// Build and run it with:
///
///     make example3 ./example3
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

// forward definition
void set_CS_weights(double theta, double phi, vector<double> & weights);

int main(int argc, char** argv) {

  double mZ  = 91.1876;

  // first prepare the cuts that we will examine
  double ymax = 2.4;

  vector<CutsBase*> cuts;
  cuts.push_back(new CutsNormal  (25.0, 25.0));
  cuts.push_back(new CutsProduct (25.0, 30.0));

  // set up a default CBI_H,DY
  auto cbiHDY = new CutsCBIPt(25.0, 30.0);
  // this sets up the cos(theta)>cbar cut (Eq. 6.7 of the paper)
  cbiHDY->set_DY_mode();
  cuts.push_back(cbiHDY);

  // for each of the hardness cuts, add in the info about the rapidity cuts
  for (auto & cut: cuts) {
    cut->add_upper_absrap_cut(ymax);
  }
    
  // choose some Higgs kinematics at which we wish to examine the acceptace
  double ptZ = 25.0;
  double yZ  = 0.4;
  cout << "Boson kinematics: ptZ = " << ptZ 
       << ", yZ = " << yZ
       << ", mZ = " << mZ << endl;

  // choose a random number generator and a number of events.
  //
  // NB: if you want to loop over different Z-boson kinematic
  //     points in order to examine the ptZ dependence, it is
  //     a good idea to use the same starting seed for each
  //     ptZ, so that statistical errors are correlated 
  //     across the pTZ points.
  unsigned nev = 10000000;
  int seed = 5;
  mt19937 generator(seed);
  auto phi_dist      = uniform_real_distribution<double>(-pi,pi);
  auto costheta_dist = uniform_real_distribution<double>(-1.0,1.0);
  
  // do a simple MC integration to evaluate the acceptance
  // for each of the cuts (using the a common set of boson
  // decays across all the cuts)
  constexpr unsigned nweights = 9;
  vector<double> weights(nweights);
  vector<vector<double>> npass(cuts.size(),vector<double>(nweights,0.0));
  for (unsigned iev = 0; iev < nev; iev++) {
    if (iev > 0 && iev % 1000000 == 0) cout << "completed " << iev << " events out of " << nev << endl;

    // choose uniform random orientation for the Z-boson decay
    double phi = phi_dist(generator);
    double theta = acos(costheta_dist(generator));

    // create boson (which by default is boosted along the x axis) and
    // determine the Collins-Soper weights
    Boson boson(mZ,theta,phi,ptZ,yZ);
    set_CS_weights(theta,phi, weights);

    // loop over cuts and weights and fill the results
    for (unsigned icut = 0; icut < cuts.size(); icut++) {
      if (cuts[icut]->pass(boson)) {
        for (unsigned iw = 0; iw < nweights; iw++) {
          npass[icut][iw] += weights[iw];
        }
      }
    }
  }

  // get the results
  for (unsigned icut = 0; icut < cuts.size(); icut++) {
    cout << "\nResults for \n" << cuts[icut]->description() << "\n";
    for (unsigned iw = 0; iw < nweights; iw++) {
      cout << "Harmonic acceptance ";
      if (iw == 8) cout << "u: ";
      else         cout << iw << ": ";
      // in the output, we include the factor of 3/4 such that
      // these accepances multiply dsigma/dq(unpol) * Ai(q)
      // (or just dsigma/dq(unpol) for i==u)
      cout << 0.75 * npass[icut][iw]/nev  << endl;
    }
  }

}

inline double pow2(double x) {return x*x;}

/// The DY cross section is given by 
///
///    dsigma/dcostheta dphi dq = 1/4pi dsigma^{unpol}/dq *
///                        3/4(weights[8] + sum_{i=0}^n weights[i] * A_i(q))
///
/// This routine takes theta and phi and fills the weights[..] vector
void set_CS_weights(double theta, double phi, vector<double> & weights){
  assert(weights.size() == 9);
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  double cos_phi   = cos(phi);

  double sin_phi = sin(phi);
  double sin_2phi = 2*sin_phi*cos_phi;
  double cos_2phi  = 2*pow(cos_phi,2)-1;
  double sin_2theta = 2*sin_theta*cos_theta;
  weights[0] = 0.5*(1 - 3*pow2(cos_theta));
  weights[1] = sin_2theta*cos_phi;
  weights[2] = 0.5*pow2(sin_theta)*cos_2phi;
  weights[3] = sin_theta * cos_phi;
  weights[4] = cos_theta;
  weights[5] = pow2(sin_theta) * sin_2phi;
  weights[6] = sin_2theta * sin_phi;
  weights[7] = sin_theta * sin_phi;
  weights[8] = 1 + pow2(cos_theta);

}
