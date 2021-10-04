// This file is part of the two-body-cuts project.
// It implements the two-body cuts from https://arxiv.org/abs/2106.08329
//
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

#ifndef __TWOBODYCUTS_HH__
#define __TWOBODYCUTS_HH__
#include "fastjet/Selector.hh"
#include "fastjet/Error.hh"
#include "Boson.hh"
#include "RealRange.hh"
#include <limits>
#include <array>
#include <limits>
#include <stdexcept>
#include <cmath>

namespace tbc {

/// Base class for cuts
class CutsBase {
public:
  CutsBase(double ptsoft, double pthard, fastjet::Selector selector = fastjet::SelectorIdentity()) : 
      _ptsoft(ptsoft), _pthard(pthard), _selector(selector) {}

  virtual bool pass(const Boson & b) const = 0;
  virtual std::string description() const = 0;

  bool pass_other(const Boson & b) const {
    return _selector.pass(b.p1) && _selector.pass(b.p2);
  }

  void get_ptminmax(const Boson & b, double & ptmin, double & ptmax) const {
      ptmin = b.p1.pt();
      ptmax = b.p2.pt();

      if (ptmin > ptmax) std::swap(ptmin, ptmax);
  }

  // add an upper cut on rapidity
  virtual void add_upper_rap_cut(double rap) {_selector &= fastjet::SelectorRapMax(rap);}
  // add a lower cut on rapidity
  virtual void add_lower_rap_cut(double rap) {_selector &= fastjet::SelectorRapMin(rap);}
  // add a rapidity exclusion band
  virtual void add_rap_exclusion(double rap1, double rap2) {
    assert(rap1 < rap2);
    _selector &= !(fastjet::SelectorRapRange(rap1,rap2));
  }

  // add an upper cut on |rapidity|
  virtual void add_upper_absrap_cut(double absrap) {
    add_upper_rap_cut( std::abs(absrap));
    add_lower_rap_cut(-std::abs(absrap));
  }
  // add an absolute rapidity exclusion band between |absrap1| and |absrap2|
  virtual void add_absrap_exclusion(double absrap1, double absrap2) {
    add_rap_exclusion( std::abs(absrap1),  std::abs(absrap2));
    add_rap_exclusion(-std::abs(absrap2), -std::abs(absrap1));
  }

  virtual ~CutsBase() {}

  const fastjet::Selector & selector() const {return _selector;}
protected:

  /// the cuts on the harder and softer decay product
  double _ptsoft, _pthard;
  /// any cuts that are universal
  fastjet::Selector _selector;

  constexpr static double pi = 3.141592653589793238462643383279502884197;
};


/// normal cuts
class CutsNormal : public CutsBase {
public:
  CutsNormal(double ptsoft, double pthard, fastjet::Selector selector = fastjet::SelectorIdentity()) : 
    CutsBase(ptsoft, pthard, selector) {}

  bool pass(const Boson & b) const override {
    // three common lines
    if (!pass_other(b)) return false;
    double ptmin, ptmax;
    get_ptminmax(b, ptmin, ptmax);

    // then the specifics of what we do
    return ptmin > _ptsoft && ptmax > _pthard;
  }

  std::string description() const override {
    std::ostringstream ostr; 
    ostr << "normal cuts: ptmin >= " << _ptsoft 
         << " && ptmax >= "  << _pthard 
         << " && a selector(" << _selector.description() << ") applied to each decay product";
    return ostr.str();
  }

};

typedef CutsNormal CutsAsymmetric;

/// staggered cuts, where pthard is applied to boson.p1
/// and ptsoft is applied to boson.p2
class CutsStaggered : public CutsBase {
public:
  CutsStaggered(double ptsoft, double pthard, fastjet::Selector selector = fastjet::SelectorIdentity()) : 
    CutsBase(ptsoft, pthard, selector) {}

  bool pass(const Boson & b) const override {
    // three common lines
    if (!pass_other(b)) return false;
    return b.p1.pt() >= _pthard && b.p2.pt() >= _ptsoft;
  }

  std::string description() const override {
    std::ostringstream ostr; 
    ostr << "staggered cuts: pt2 >= " << _ptsoft 
         << " && pt1 >= "  << _pthard 
         << " && a selector(" << _selector.description() << ") applied to each decay product";
    return ostr.str();
  }

};

/// sum cuts where we require 
/// 0.5 * (boson.p1.pt() + boson.p2.pt()) > pthard
/// and min(boson.p1.pt(), boson.p2.pt()) > ptsoft
class CutsSum : public CutsBase {
public:
  CutsSum(double ptsoft, double pthard, fastjet::Selector selector = fastjet::SelectorIdentity()) : 
    CutsBase(ptsoft, pthard, selector) {}

  bool pass(const Boson & b) const override {
    // three common lines
    if (!pass_other(b)) return false;
    double ptmin, ptmax;
    get_ptminmax(b, ptmin, ptmax);

    // then the specifics of what we do
    return ptmin > _ptsoft && 0.5*(ptmin+ptmax) > _pthard;
  }

  std::string description() const override {
    std::ostringstream ostr; 
    ostr << "Sum (balanced) cuts: ptmin >= " << _ptsoft 
         << " && 0.5(ptmin+ptmax) >= "  << _pthard 
         << " && a selector(" << _selector.description() << ") applied to each decay product";
    return ostr.str();
  }
};


/// product cuts where we require 
/// sqrt(boson.p1.pt()*boson.p2.pt()) > pthard
/// and min(boson.p1.pt(), boson.p2.pt()) > ptsoft
class CutsProduct : public CutsBase {
public:
  CutsProduct(double ptsoft, double pthard, fastjet::Selector selector = fastjet::SelectorIdentity()) : 
    CutsBase(ptsoft, pthard, selector) {}

  bool pass(const Boson & b) const override {
    // three common lines
    if (!pass_other(b)) return false;
    double ptmin, ptmax;
    get_ptminmax(b, ptmin, ptmax);

    // then the specifics of what we do
    return ptmin > _ptsoft && ptmin*ptmax > _pthard*_pthard;
  }

  std::string description() const override {
    std::ostringstream ostr; 
    ostr << "product cuts: ptmin >= " << _ptsoft 
         << " && sqrt(ptmin*ptmax) >= "  << _pthard 
         << " && a selector(" << _selector.description() << ") applied to each decay product";
    return ostr.str();
  }

};

/// cuts on the decay pt in the Collins-Soper frame
/// where we require
/// boson.ptCS() > pthard
/// and min(boson.p1.pt(), boson.p2.pt()) > ptsoft
class CutsPtCS : public CutsBase {
public:
  CutsPtCS(double ptsoft, double pthard, fastjet::Selector selector = fastjet::SelectorIdentity()) : 
    CutsBase(ptsoft, pthard, selector) {}

  bool pass(const Boson & b) const override {
    // three common lines
    if (!pass_other(b)) return false;
    double ptmin, ptmax;
    get_ptminmax(b, ptmin, ptmax);

    if (ptmin < _ptsoft) return false;

    return 0.25*b.decay_sinsq_theta()*b.m2() >= _pthard*_pthard;
  }

  std::string description() const override {
    std::ostringstream ostr; 
    ostr << "Cuts using the pt defined in the Collins-Soper frame: ptmin(lab) > " << _ptsoft 
         << " && ptCS > "  << _pthard 
         << " && a selector(" << _selector.description() << ") applied to each decay product";
    return ostr.str();
  }
};


/// Compensating Boost Invariant cuts, where the compensation works on
/// the hardness, but not as concerns rapidity.
///
/// At low Boson transverse momenta (and away from rapidity points where
/// pairs of rapidity and/or ptCS cuts are simultaneously relevant), this
/// class gives results that are identical to CutsPtCS
class CutsCBIPt : public CutsBase {
public:
  CutsCBIPt(double ptsoft, double pthard, fastjet::Selector selector = fastjet::SelectorIdentity()) : 
    CutsBase(ptsoft, pthard, selector) {}

  void set_DY_mode(bool val = true) {_DY_mode = val;}

  bool pass(const Boson & b) const override {
    // three common lines
    if (!pass_other(b)) return false;
    double ptmin, ptmax;
    get_ptminmax(b, ptmin, ptmax);

    if (ptmin < _ptsoft) return false;

    double sinsq_theta = b.decay_sinsq_theta();
    if (_DY_mode) {
      double costheta = sqrt(1 - sinsq_theta);
      double cbar = DY_cbar_from_pthard(b);
      if (costheta < cbar) return false;
    }
    
    bool pass_ptCS = 0.25 * sinsq_theta * b.m2() >= _pthard*_pthard;
    if (pass_ptCS) return true;


    // now try to handle the delicate case where the hard cut does not pass

    // First, we will need the decay azimuth in the Collins-Soper frame
    // (in a range 0..pi/2)
    double phi = b.decay_phi_0topi2();

    // if we are in a region of phi where the soft cut can never play a
    // role at this pt (the corresponding sqrt_contents are negative), 
    // then we pass the cuts
    // GREED1
    double sqrt_contents_phi = sqrt_contents(b,phi);
    if (sqrt_contents_phi < 0) return true;

    // some useful shorthands
    double mB = b.m();
    double mB2 = pow(mB,2);
    double ptB = b.pt();
    double ptB2 = pow(ptB,2);
    double mtB = b.mt();

    // next we check if the second sin(theta) solution is physical; if
    // it is, and if we are below it (i.e. above the corresponding
    // cos(theta) value), then we accept the event
    // GREED2
    double sintheta_soft_cut_soln2 = (2*ptB*mtB*cos(phi) - sqrt(sqrt_contents_phi))
                  / (2*mB2 + ptB2*(1+cos(2*phi)));
    //std::cout << "CBIPt: " << sintheta_soft_cut_soln2 << " " << sqrt(sinsq_theta) << std::endl;
    if (sintheta_soft_cut_soln2 > 0 && sinsq_theta < pow(sintheta_soft_cut_soln2,2)) return true;


    // Preparatory steps for the mirroring
    // when phi>pi/4 larger, we look up how the cut would behave in the
    // opposite "mirror" region
    double mphi = pi/2 - phi;
    // work out where the soft cut would be for the mirror phi
    // [pT2Minussin\[Theta]Soln from boostinv-cuts-extra.nb]
    double sqrt_contents_mphi = sqrt_contents(b,mphi);
    // this is part of GREED1
    if (sqrt_contents_mphi < 0) return true;

    // Next we consider the balancing; when phi is less than pi/4 we
    // don't apply any further correction, because when phi<=pi/4
    // typically we are losing acceptance from the interplay of ptmin
    // and ptCS cuts rather than gaining it, so there is generally
    // nothing to add back.
    if (phi <= pi/4) return false;


    // work out where the (ptCS) hard cut would be in cos(theta)
    double sintheta_hard_cut = 2 * _pthard / b.m();
    double costheta_hard_cut = sqrt(1-pow(sintheta_hard_cut,2));

    // if sqrt_contents_mphi < 0 here (phi>pi/4, mphi < pi/4), 
    // then we already had sqrt_contents_phi < 0 and accepted the event
    // NB: assumes GREED1
    if (sqrt_contents_mphi < 0) {
      std::ostringstream ostr;
      ostr << "Did not expect sqrt_contents_mphi < 0, but it was " << sqrt_contents_mphi;
      throw std::range_error(ostr.str());
    }
    double msintheta_soft_cut = (2* ptB*mtB*cos(mphi) + sqrt(sqrt_contents_mphi))
                  / (2*mB2 + ptB2*(1+cos(2*mphi)));
    // if this cannot be interpreted as a sin(theta), i.e. we have
    // completely excluded high values of sin(theta) (low values of cos)
    // on the mirror side, then we accept the event on the original side
    // (NOT GREED???)
    if (msintheta_soft_cut >= 1) return true;
    // translate this into a cos(theta) cut
    double mcostheta_soft_cut = sqrt(1 - pow(msintheta_soft_cut,2));
    // and then on the original side, introduce a new cos(theta) cut that
    // compensates for missing region on the mirror side
    //std::cout << "CBIPt mcostheta_soft_cut:"<< mcostheta_soft_cut << "\n";
    if (mcostheta_soft_cut < costheta_hard_cut) {
      double new_costheta_cut = costheta_hard_cut + (costheta_hard_cut-mcostheta_soft_cut);
      double costheta = sqrt(1-sinsq_theta);
      //std::cout << "CBIPt new_costheta_cut:"<< new_costheta_cut << "\n";
      return costheta < new_costheta_cut;
    }

    // otherwise, return false
    return false;

  }

  /// returns the contents of the sqrt that appears
  /// in the determination of the sin(theta) associated
  /// with the soft cut at a given phi
  double sqrt_contents(const Boson & b, double phi_in) const {
    double mB2 = b.m2();
    double ptB = b.pt();
    double ptB2 = pow(ptB,2);
    double ptsoft2 = pow(_ptsoft,2);
    return 2*(8*mB2*ptsoft2 - (mB2 - 4*ptsoft2)*ptB2 + (mB2+4*ptsoft2)*ptB2 * cos(2*phi_in));
  }

  std::string description() const override {
    std::ostringstream ostr; 
    ostr << "Compensating Boost Invariant cuts (compensation uses only hardness): ptmin > " << _ptsoft 
         << " && ptCS > "  << _pthard << " (unless compensated)"
         << " && a selector(" << _selector.description() << ") applied to each decay product";
    if (_DY_mode) ostr << ", DY mode turned on";
    return ostr.str();
  }

  /// works out the cbar given the boson and the known
  /// thard cut (or 0 if the c cut is out of range)
  double DY_cbar_from_pthard(const Boson & b) const {
    double s2 = pow(2*_pthard,2) / b.m2();
    if (s2 > 1.0) return 0.0;
    return DY_cbar(sqrt(1-s2));
  }

  /// give a cut cos(theta) < c, returns the cbar value that gives a
  /// zero f^(1) harmonic acceptance when additionally imposing
  /// cos(theta)>cbar. 
  ///
  /// The resulting range is non-zero only for c > sqrt(1/3) (and when 
  /// the range would be negative it returns cbar=c).
  double DY_cbar(double c) const {
    double sqrt_arg = 4 - 3*c*c;
    if (sqrt_arg < 0) return c;
    double cbar = (-c + sqrt(sqrt_arg))/2.0;
    if (cbar > c) return c;
    else          return cbar;
  }

private:
  bool _DY_mode = false;
};


/// Compensating Boost Invariant cuts, where the compensation takes into
/// account both the hardness and rapidity cuts.
///
/// At low Boson transverse momenta (and away from rapidity points where
/// pairs of rapidity and/or ptCS cuts are simultaneously relevant),
/// this class gives results that are identical to CutsPtCS.
///
/// One can optionally transition to a higher ptCS cut beyond some
/// rapidity. This can be necessary in particular if there is a region
/// midway between two rapidity cuts where both cuts affect the
/// acceptance in an identical way (excluding large values of
/// |cos(theta)|): there the compensation has no room to work unless
/// one increases the ptCS cut.
///
class CutsCBIPtRap : public CutsBase {
public:

  typedef RealRange::Segment Segment;

  /// constructor for CBIPtRap. The ptsoft cut is applied to the pt's of
  /// both decay products, while the pthard cut is applied to ptCS and
  /// can be adjusted internally to achieve compensation for the
  /// acceptance. 
  ///
  /// Greed levels specify choices about the behaviour that are relevant
  /// only for ptB > 2*ptsoft. There are four greed levels, with
  /// acceptance increasing as the greed level goes up:
  ///
  ///   0. just does what it can to compensate across the original
  ///   decay_phi and the mirrored values
  ///
  ///   1. for boson+decay kinematics where the ptsoft cut is inactive
  ///   for all decay theta values, for a given decay phi or for any of
  ///   its mirrored values, the pthard cut is removed completely.
  ///
  ///   2. additionally, for boson+decay kinematics where the ptsoft
  ///   cut's upper boundary in cos(theta) is not 1 (call the boundary
  ///   CSU), ignore the allowed region between CSU and 1 in the
  ///   compensation calculations (for this phi and the mirrored ones)
  ///   and accept the event if this boson's cos(theta) is above CSU.
  ///
  ///   3. additionally, for boson+decay kinematics if any of this phi
  ///   or the mirrored phi values has the property that the  upper
  ///   boundary in cos(theta) is not 1, accept the event.
  ///
  /// In the absence of rapidity cuts, the default value of greed 2
  /// gives identical results to the CutsCBIPt class.
  CutsCBIPtRap(double ptsoft, double pthard, int greed = default_greed) 
            : CutsBase(ptsoft,pthard), _greed(greed) {
    // a negative greed value 
    if (greed < 0) _greed = default_greed;
  }

  static constexpr int default_greed = 2;

  /// override the base class functions for registering because we need to store the
  /// info in a more granular form 
  void add_upper_rap_cut(double rap) override {_upper_rap_cut = rap; _selector &= fastjet::SelectorRapMax(rap);}
  void add_lower_rap_cut(double rap) override {_lower_rap_cut = rap; _selector &= fastjet::SelectorRapMin(rap);}
  void add_rap_exclusion(double rap1, double rap2) override {
    assert(rap1 < rap2);
    _rap_exclusions.push_back({rap1,rap2});
    _selector &= !(fastjet::SelectorRapRange(rap1,rap2));
  }

  /// adds an optional higher pt cut 
  void add_higher_pthard(double from_abs_rap, double higher_pthard) {
    _abs_rap_for_higher_pthard = from_abs_rap;
    _higher_pthard = higher_pthard;
  }

  bool pass(const Boson & b) const override {
    // first apply the non-negotiable cuts
    if (!pass_other(b)) return false;
    double ptmin = std::min(b.p1.pt(), b.p2.pt());
    if (ptmin < _ptsoft) return false;

    double pt_CS = b.pt_CS();
    bool pass_pthard =  pthard_min_at_yB(b.rap()) <= pt_CS && pt_CS <= pthard_max_at_yB(b.rap());
    RealRange extra_veto, extra_allowed;
    get_ptrap_extra_veto_or_allowed(b, extra_veto, extra_allowed);

    //cout << pass_pthard << " " << extra_veto << " " << extra_allowed << endl;
    if (extra_veto.size() > 0) {
      return pass_pthard && !extra_veto.contains(b.decay_cos_theta());
    } else if (extra_allowed.size() > 0) {
      return pass_pthard || extra_allowed.contains(b.decay_cos_theta());
    } else {
      return pass_pthard;
    }
  }


  /// returns the pthard value (ie. the cut on ptCS) to be used at this boson rapidity
  virtual double pthard_min_at_yB(double yB) const {
    if (std::abs(yB) < _abs_rap_for_higher_pthard || _higher_pthard < 0.0) return _pthard;
    else return _higher_pthard;
  }

  void get_ptrap_extra_veto_or_allowed(const Boson & b, 
                             RealRange & extra_veto, RealRange & extra_allowed) const {
    extra_veto.reset();
    extra_allowed.reset();

    // initial setup
    double yB  = b.rap();
    double ptB = b.pt();
    double mB  = b.m();
    double phi = b.decay_phi_0topi();
    double pi_phi = pi - phi;
    // the four phi values we will investigate
    // (i.e. the original one and the reflected ones)
    std::array<double,4> phivals {{phi, pi_phi, pi/2 - std::min(phi,pi_phi), pi/2 + std::min(phi,pi_phi)}};
    // and the fixed allowed regions that are independent of phi
    double sinhard_max = 2 * (pthard_max_at_yB(yB) / mB);
    double coshard_min = sinhard_max <= 1.0 ? sqrt(1-pow(sinhard_max,2)) : 0.0;
    RealRange pthard_allowed = Segment(coshard_min, sqrt(1 - pow(2*pthard_min_at_yB(yB)/mB,2)));;

    // pt0_allowed is the allowed extent for pt = 0, with all cuts included
    RealRange pt0_allowed = rap_allowed_costheta(yB, 0.0, mB, phivals[0]) && pthard_allowed;
    double pt0_extent = pt0_allowed.extent();
    //cout << "pt0_extent=" << pt0_extent << endl;

    // next we identify the allowed regions at the four phi values, with
    // and without the pthard cut
    std::array<RealRange,4> rap_ptmin_allowed, rap_ptminmax_allowed;
    std::array<double,4> rap_ptmin_extra, rap_ptminmax_extra;
    double total_rap_ptminmax_extra = 0.0;
    for (unsigned int i = 0; i < 4; i++) {
      // first understand the ptsoft_vetoed region; 
      Segment ptsoft_vetoed = ptsoft_vetoed_costheta(_ptsoft, ptB, mB, phivals[i]);
      //cout << i << " " << phivals[i] << " " << ptsoft_vetoed << " " << _ptsoft << endl;

      // Greed levels only affect the behaviour for ptB >= 2 psoft
      //
      // With greed-level 1 or higher, if, for any of the phi values,
      // the extent is zero, then we don't apply rapidity compensation
      // here and instead just accept the event
      if (ptsoft_vetoed.extent() == 0 && _greed >= 1) {
        extra_allowed = RealRange(0,1);
        return;
      }
      if (ptsoft_vetoed.extent() != 0 && ptsoft_vetoed.hi() < 1 && _greed >= 2) {
        // if the upper end of the vetoed region is not at 1, then
        // greed-level 3 just accepts the event, regardless of whether
        // we trigger this for the original phi or a reflected one
        if (_greed == 3) {extra_allowed = RealRange(0,1); return;}

        // otherwise (the case of greed == 2), always accept the
        // cos(theta) range above ptsoft_vetoed.hi(), regardless of how
        // it enters into the calculations below. [NB: but calling
        // routine ignores this if there is something to veto... ]
        assert(_greed == 2);
        //cout << "AG " << phivals[i]/pi << " " << b.decay_cos_theta()  << " ptsoft_vetoed=" << ptsoft_vetoed << endl;
        if (i == 0 && b.decay_cos_theta() > ptsoft_vetoed.hi()) {
          extra_allowed = RealRange(0,1); return;
        }
      }

      /// The following line ensures that for greed-level 2, we reproduce
      /// the behaviour of the CBIpt cuts
      if (_greed >= 2) ptsoft_vetoed.reset(ptsoft_vetoed.lo(), 1.0);

      // for the purposes of rapidity compensation only use the lower part?
      //RealRange ptsoft_allowed = !RealRange(Segment(ptsoft_vetoed.lo(), 1.0));
      RealRange ptsoft_allowed = !RealRange(Segment(ptsoft_vetoed));

      rap_ptmin_allowed[i]    = rap_allowed_costheta(yB, ptB, mB, phivals[i]) && ptsoft_allowed;
      rap_ptminmax_allowed[i] = rap_ptmin_allowed[i] && pthard_allowed;
      rap_ptmin_extra[i]      = rap_ptmin_allowed[i]   .extent() - pt0_extent;
      rap_ptminmax_extra[i]   = rap_ptminmax_allowed[i].extent() - pt0_extent;
      total_rap_ptminmax_extra += rap_ptminmax_extra[i];
      //cout << i << " phi=" << phivals[i] << " " << rap_ptmin_allowed[i] << " " << rap_ptminmax_allowed[i] 
      //    << " extra=" << rap_ptminmax_extra[i] << endl;
    }  
    //cout << "phi,total_rap_ptminmax_extra=" << phi << " " << total_rap_ptminmax_extra << " " << rap_ptmin_allowed[0] << endl;

    //cout << "AG total_rap_ptminmax_extra = " << total_rap_ptminmax_extra << endl;
    if (total_rap_ptminmax_extra > 0) {
      // we need to introduce vetoes and we will distribute the vetoes
      // accoring to the relative amount of extra allowed phase space
      // relative to pt=0 at each phi value (ignoring those with no
      // extra allowed phase space)
      double denom = 0.0;
      for (unsigned i = 0; i < 4; i++) denom += std::max(rap_ptminmax_extra[i],0.0);
      double thiswgt = std::max(rap_ptminmax_extra[0],0.0) / denom;
      double amount_to_veto = thiswgt * total_rap_ptminmax_extra;
      if (amount_to_veto >0) extra_veto = rap_ptminmax_allowed[0].veto_amount_from_top(amount_to_veto);
    } else {
      // we need to recover whatever we can
      double amount_to_recover = -total_rap_ptminmax_extra;
      // first work out the maximum that can be recovered from each phi (and their total)
      double total_recoverable = 0.0;
      std::array<double,4> recoverable;
      for (unsigned i = 0; i < 4; i++) {
        recoverable[i] = rap_ptmin_extra[i] - rap_ptminmax_extra[i];
        total_recoverable += recoverable[i];
        //cout << "AG recoverable " << i << " " << recoverable[i] << endl;
      }
      double to_recover_here;
      // // diagnostics, that we probably want to promote to a proper test suite
      // cout << "phi, recoverable, to recover = " << phi << " " << total_recoverable << "(";
      // for (auto r: recoverable)  cout << r << ",";
      // cout << ") " <<  amount_to_recover  << "(";
      // for (auto r: rap_ptminmax_extra)  cout << r << ",";
      // cout << ")" << endl;
      if (total_recoverable < amount_to_recover) {
        // if we can't recover everything we need, then just recover the max
        // from this phi
        to_recover_here = recoverable[0];
      } else {
        // otherwise recover in proportion to the amount that is recoverable
        // at each phi
        to_recover_here = amount_to_recover * recoverable[0] / total_recoverable;
      }
      if (to_recover_here > 0.0) {
        RealRange range_for_recovery = rap_ptmin_allowed[0] && (!pthard_allowed);
        extra_allowed |= range_for_recovery.extra_amount_from_bottom(to_recover_here);
      }
      //std::cout << "AG total_recoverable, extra_allowed = " << total_recoverable << " " << extra_allowed << " " << amount_to_recover<< endl;
    }

  }

  /// return a default value that is large enough that unlikely ever to
  /// matter physically, but not so large that we can't safely square it, etc.
  virtual double pthard_max_at_yB(double yB) const {return _pthard_max;}

  void set_pthard_max(double value) {_pthard_max = value;}

  /// returns the segment of cos(theta) values in which the ptsoft cut vetoes decays
  /// of a boson with the specified kinematic variables
  static Segment ptsoft_vetoed_costheta(double ptsoft, double ptB, double mB, double phi) {
    double mB2 = pow(mB,2);
    double ptB2 = pow(ptB,2);
    double mtB = sqrt(mB2 + ptB2);
    double ptsoft2 = pow(ptsoft,2);
    double cos2phi = cos(2*phi);

    double sqrt_contents = 2*(8*mB2*ptsoft2 - (mB2 - 4*ptsoft2)*ptB2 + (mB2+4*ptsoft2)*ptB2 * cos2phi);
    // no veto when the sqrt contents are negative
    if (sqrt_contents < 0) return Segment(1,1);
    double sqrt_result = sqrt(sqrt_contents);

    double denom = (2*mB2 + ptB2*(1+cos2phi));
    // NB abs(cos(phi)) ensures that we don't need to care whether phi is in the
    // range 0..pi/2 or pi/2..pi; or equivalently, it means we are always cutting on the 
    // softer particle
    double num_extra = 2*ptB*mtB*abs(cos(phi));
    double sintheta_plus  = (num_extra + sqrt_result)/denom;
    double sintheta_minus = (num_extra - sqrt_result)/denom;

    // now convert the sintheta solutions into sensible costheta solutions
    double costheta_lower = sintheta_plus < 1.0 ? sqrt(1-pow(sintheta_plus,2))  : 0;
    double costheta_upper;
    if      (sintheta_minus < 0.0) costheta_upper = 1.0;
    else if (sintheta_minus > 1.0) costheta_upper = 0.0;
    else                           costheta_upper = sqrt(1-pow(sintheta_minus,2));

    return Segment(costheta_lower, costheta_upper);
  }

  /// returns the allowed costheta range after applying the rapidity cuts
  RealRange rap_allowed_costheta(double yB, double ptB, double mB, double phi) const {
    RealRange result(0,1);
    if (_upper_rap_cut < _no_rap_cut) {
      if (yB >= _upper_rap_cut) return RealRange();
      result &= rap_ymax_allowed_costheta(_upper_rap_cut, yB, ptB, mB, phi);
    }
    if (_lower_rap_cut > -_no_rap_cut) {
      if (yB <= _lower_rap_cut) return RealRange();
      result &= rap_ymax_allowed_costheta(_lower_rap_cut, yB, ptB, mB, phi);
    }

    for (const auto & excl: _rap_exclusions) {
      double ylo = excl.first;
      double yhi = excl.second;
      if (yB < ylo) {
        result &=  rap_ymax_allowed_costheta(ylo, yB, ptB, mB, phi) 
              || (!rap_ymax_allowed_costheta(yhi, yB, ptB, mB, phi) );
      } else if (yB < yhi) {
        result &= !(rap_ymax_allowed_costheta(ylo, yB, ptB, mB, phi))
               && !(rap_ymax_allowed_costheta(yhi, yB, ptB, mB, phi));
      } else {
        result &= (!rap_ymax_allowed_costheta(ylo, yB, ptB, mB, phi) )
              ||   rap_ymax_allowed_costheta(yhi, yB, ptB, mB, phi);
      }
    }
    return result;
  }

  /// returns a RealRange for the allowed costheta when forcing a decay
  /// product to have 
  /// - y < ymax (y > ymax), 
  /// - given a boson at rapidity yB < y (yB > y), 
  /// - with pt = ptB mass=mB, and a CS decay phi.
  /// 
  static RealRange rap_ymax_allowed_costheta(double ymax, double yB, double ptB, double mB, double phi) {

    // set up a whole bunch of shorthands
    double ycut = ymax - yB;
    double csch = 1.0 / sinh(ycut);
    double cos2phi = cos(2*phi);
    double cosphi  = cos(phi);
    double ptB2  = ptB * ptB;
    double mB2   = mB * mB;
    double csch2 = csch * csch; 

    // sort out signs for case where yB < ymax (+ special yB=ymax case)
    if      (ycut <  0) cosphi = -cosphi;
    else if (ycut == 0) return RealRange();

    // now start the calculations
    double sqrt_arg = ptB2 * (cos2phi - 1.0) + 2* mB2 * csch2;
    sqrt_arg *= 2 * mB2 * (1 + csch2);

    // full range is allowed when the sqrt_arg < 0 
    // (because there is no costheta value that intersects
    // our region)
    //cout << "sqrt_arg ="  << sqrt_arg << endl;
    if (sqrt_arg <= 0) return RealRange(0.0, 1.0);

    double sqrt_res = sqrt(sqrt_arg);
    double denom = ptB2 * (1 + cos2phi) + 2*mB2 * (1 + csch2);
    double num1  = -2*ptB*sqrt(mB2 + ptB2) * cosphi;

    // get the sintheta solutions
    double sintheta1 = (num1 - sqrt_res)/denom;
    double sintheta2 = (num1 + sqrt_res)/denom;

    // put them in a sensible range
    //cout << "sintheta1,2 = " << sintheta1 << " " << sintheta2 << endl;
    if (sintheta1 < 0) sintheta1 = 0;
    if (sintheta2 < 0) sintheta2 = 0;
    if (sintheta1 > 1) sintheta1 = 1;
    if (sintheta2 > 1) sintheta2 = 1;

    // convert to cos theta
    double costheta1 = sqrt(1.0 - pow(sintheta1,2));
    double costheta2 = sqrt(1.0 - pow(sintheta2,2));

    return !RealRange(costheta2, costheta1);

  }



  std::string description() const override {
    std::ostringstream ostr;
    ostr << "Compensating Boost Invariant cuts (compensation accounts for hardness and rapidity) with ptmin > " << _ptsoft
         << ", ptCS > " << _pthard << "(flexible for compensation)";
    if (_higher_pthard >= 0) ostr << " (raised to " << _higher_pthard
                                  << " for |yB| > " << _abs_rap_for_higher_pthard << ")";
    if (_pthard_max != _pthard_max_default) ostr << ", ptCS < " << _pthard_max;
    ostr << ", rapidity cuts = " << _selector.description() 
         << ", all supplemented with pt and rapidity compensation"
         << " with greed level " << _greed;
    return ostr.str();
  }

private:
  int _greed;
  static constexpr double _pthard_max_default = 1e100;
  double _pthard_max = _pthard_max_default;

  static constexpr double _no_rap_cut = 1e20;

  /// our internal storage of rapidity cuts
  double _upper_rap_cut =   _no_rap_cut;
  double _lower_rap_cut = - _no_rap_cut;
  std::vector<std::pair<double, double> > _rap_exclusions;

  /// absolute rapidity from which we apply a higher pthard cut
  double _abs_rap_for_higher_pthard=_no_rap_cut;
  /// value of that higher pthard cut; negative values turn off the functionality
  double _higher_pthard = -1.0;

};

/// A class that implements cuts on the pt of the softer boson decay
/// product and on the rapidities of both decay products. It provides
/// information on the acceptance, which can be used to defiducialise
/// the acceptance for a scalar boson decay.
///
/// This class should be used as follows
///
///     CutsDefid cuts_defid(ptsoft_min); 
///     // [... add rapidity cuts usual ...] 
///     // Decide whether you want to normalise to the Born acceptance. 
///     // (see scalar_acceptance_four_phi(...) documentation for details)
///     bool normalise_to_Born = true;
///
///     // then, inside event loop: 
///     bool accept = cuts_defid.pass(boson); 
///     if (accept) {
///       double weight = 1.0 / cuts_defid.scalar_acceptance_four_phi(boson, normalise_to_Born);
///       // then bin the event with the given weight
///     }
///
/// If you use the idea of defiducialisation, you should refer also to
/// A. Glazov, arXiv:2001.02933 (Eur.Phys.J.C 80 (2020) 9, 875).
///
/// The idea of normalising to a rapidity-dependent result was suggested
/// by the referee of arXiv:2106.08329.
///
class CutsDefid4Phi : public CutsCBIPtRap {
public:

  CutsDefid4Phi(double ptsoft) : CutsCBIPtRap(ptsoft, 0.0, CutsCBIPtRap::default_greed) {}

  bool pass(const Boson & b) const override {
    // first apply the non-negotiable cuts
    if (!pass_other(b)) return false;
    double ptmin = std::min(b.p1.pt(), b.p2.pt());
    if (ptmin < _ptsoft) return false;

    return true;
  }

  std::string description() const override {
    std::ostringstream ostr; 
    ostr << "Defiducialisable (4-phi) symmetric cut with ptmin >= " << _ptsoft 
         << " and rapidity cuts = " << _selector.description();
    return ostr.str();
  }

  /// Returns the acceptance for this value of the boson pt, rapidity,
  /// averaged across Collins-Soper polar decay angle at this and the 3
  /// other mirror phi values. 
  ///
  /// If normalise_to_Born is false, then the result is such that for a
  /// scalar decay, after integration of the decay angles, the events
  /// weights (calculated as one over the return value of this function)
  /// exactly cancel the acceptance. This can be used for obtaining
  /// total Higgs cross section in a given Higgs rapidity window,
  /// without the need to know anything about the rapidity distribution
  /// of the Higgs boson. Beware, this should not be used for Higgs
  /// bosons close to the maximum allowed photon pseudorapidity, because
  /// the acceptance tends to zero, and so the event weights would
  /// diverge
  ///
  /// If normalise_to_Born is true, the result is returned normalised to
  /// the Born acceptance for that rapidity. This can be useful so as to
  /// avoid giving a very large weight to events that are close to the
  /// edge of the rapidity region. 
  ///
  double scalar_acceptance_four_phi(const Boson & b, bool normalise_to_Born = false) const {
    double yB  = b.rap();
    double ptB = b.pt();
    double mB  = b.m();
    double phi = b.decay_phi_0topi();
    double pi_phi = pi - phi;
    // the four phi values we will investigate
    // (i.e. the original one and the reflected ones)
    constexpr unsigned nphi = 4;
    std::array<double,nphi> phivals {{phi, pi_phi, pi/2 - std::min(phi,pi_phi), pi/2 + std::min(phi,pi_phi)}};

    double acceptance_sum = 0.0;
    for (unsigned i = 0; i < nphi; i++) {
      RealRange ptmin_allowed = !ptsoft_vetoed_costheta(_ptsoft, ptB, mB, phivals[i]);
      RealRange ptmin_rap_allowed = rap_allowed_costheta(yB, ptB, mB, phivals[i]) && ptmin_allowed;

      acceptance_sum += ptmin_rap_allowed.extent();
    }

    if (! normalise_to_Born) {
      return acceptance_sum / nphi;
    } 

    // to get a result normalised to the Born acceptance, we can use the fact that
    // the Born acceptance does not depend on phi, so we need only consider one 
    // of the phi values
    RealRange born_ptmin_allowed = !ptsoft_vetoed_costheta(_ptsoft, 0.0, mB, phivals[0]);
    RealRange born_ptmin_rap_allowed = rap_allowed_costheta(yB, 0.0, mB, phivals[0]) && born_ptmin_allowed;
    double born_acceptance_sum = born_ptmin_rap_allowed.extent();
    return acceptance_sum / (4*born_acceptance_sum);
  }

  /// returns the acceptance for a scalar decay for a "Born" boson, i.e.
  /// one with zero pt but the same rapidity as this boson
  double scalar_acceptance_Born(const Boson & b) const {
    double yB  = b.rap();
    double mB  = b.m();
    double phi = 0.0;
    RealRange born_ptmin_allowed = !ptsoft_vetoed_costheta(_ptsoft, 0.0, mB, phi);
    RealRange born_ptmin_rap_allowed = rap_allowed_costheta(yB, 0.0, mB, phi) && born_ptmin_allowed;
    return born_ptmin_rap_allowed.extent();
  }


};

} // end of tbc namespace

#endif // __TWOBODYCUTS_HH__
