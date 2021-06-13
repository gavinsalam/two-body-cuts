#ifndef __BOSON_HH__
#define __BOSON_HH__
#include "fastjet/PseudoJet.hh"


class Boson;
std::ostream & operator<<(std::ostream & ostr, const Boson & b);


class Boson : public fastjet::PseudoJet {
public:
  /// create a boson at zero pt and zero rapidity
  /// with decay kinematics given by theta and phi,
  /// defined in the Collins-Soper frame
  Boson(double mass, double theta, double phi) {
    init(mass, theta, phi);
  }  

  /// create a boson the specified pt (in the x direction) and rapidity
  /// as well as decay kinematics
  Boson(double mass, double theta, double phi, double ptx, double rap) {
    init(mass, theta, phi);
    fastjet::PseudoJet boost_pt = fastjet::PtYPhiM(ptx,   0, 0, mass);
    fastjet::PseudoJet boost_y  = fastjet::PtYPhiM(0  , rap, 0, mass);
    boost(boost_pt);
    boost(boost_y);
  }  

  /// create a boson from the two input fastjet::PseudoJets
  Boson(const fastjet::PseudoJet & p1_in, const fastjet::PseudoJet & p2_in) {
    p1 = p1_in;
    p2 = p2_in;
    reset_momentum(sum());
  }

  const fastjet::PseudoJet & operator()() const {
    return *this;
  }

  Boson & boost(const fastjet::PseudoJet & boost_vector) {
    p1.boost(boost_vector);
    p2.boost(boost_vector);
    reset_momentum(sum());
    return *this;
  }
  Boson & unboost(const fastjet::PseudoJet & boost_vector) {
    p1.unboost(boost_vector);
    p2.unboost(boost_vector);
    reset_momentum(sum());
    return *this;
  }

  fastjet::PseudoJet sum() {
    return p1+p2;
  }

//  Boson z_mirrored_copy() const {
//    // determine the boosts
//    fastjet::PseudoJet rap_boost = PtYPhiM(0, rap(), 0, 1);
//    fastjet::PseudoJet pt_boost = fastjet::PseudoJet(px()/m(), py()/m(),0,1.0);
//    // take a copy and get it into the rest frame
//    Boson copy1 = *this;
//    copy1.unboost(rap_boost);
//    copy1.unboost(pt_boost);
//    // get a 2nd copy where we have mirrored the z components
//    Boson copy2(
//      fastjet::PseudoJet(copy1.p1.px(), copy1.p1.py(), -copy1.p1.pz(), copy1.p1.E()),
//      fastjet::PseudoJet(copy1.p2.px(), copy1.p2.py(), -copy1.p2.pz(), copy1.p2.E())
//    );
//    // boost it back to the proper pt and rapidity
//    copy2.boost(pt_boost);
//    copy2.boost(rap_boost);
//    return copy2;
//  }
//
//  Boson zcos_mirrored_copy(double abs_cos_threshold) const {
//    // determine the boosts
//    fastjet::PseudoJet rap_boost = PtYPhiM(0, rap(), 0, 1);
//    // take a copy and get it into the rest frame
//    Boson copy1 = *this;
//    copy1.unboost(rap_boost);
//    fastjet::PseudoJet pt_boost = copy1;
//    copy1.unboost(pt_boost);
//    double rest_costheta = copy1.p1.pz()/copy1.p1.modp();
//    double rest_phi      = atan2(copy1.p1.py(), copy1.p1.px());
//    // first mirror cos around abs_cos_threshold
//    double new_costheta = rest_costheta;
//    if (std::abs(rest_costheta) > abs_cos_threshold) {
//      // watch out when abs_cos_threshold < 1/2
//      new_costheta = abs_cos_threshold - (std::abs(rest_costheta) - abs_cos_threshold);
//      if (rest_costheta < 0) new_costheta = -new_costheta;
//    } 
//    // and now do the z mirroring
//    new_costheta = -new_costheta;
//
//    // get a 2nd copy where we have mirrored the z components
//    Boson copy2(copy1.m(), acos(new_costheta), rest_phi);
//    copy2.boost(pt_boost);
//    copy2.boost(rap_boost);
//    // cout << copy1 << endl;
//    // cout << "---- " << rest_costheta << " " << sqrt(1-decay_sinsq_theta()) << " " << sqrt(1-copy1.decay_sinsq_theta()) << " " << new_costheta << " " << sqrt(1-copy2.decay_sinsq_theta()) << endl;
//    return copy2;
//    //   fastjet::PseudoJet(copy1.p1.px(), copy1.p1.py(), -copy1.p1.pz(), copy1.p1.E()),
//    //   fastjet::PseudoJet(copy1.p2.px(), copy1.p2.py(), -copy1.p2.pz(), copy1.p2.E())
//    // );
//    // // boost it back to the proper pt and rapidity
//    // copy2.boost(pt_boost);
//    // copy2.boost(rap_boost);
//    // return copy2;
//  }
//
//
//  Boson xy_mirrored_copy() const {
//    // determine the boosts
//    fastjet::PseudoJet rap_boost = PtYPhiM(0, rap(), 0, 1);
//    fastjet::PseudoJet pt_boost = fastjet::PseudoJet(px()/m(), py()/m(),0,1.0);
//    // take a copy and get it into the rest frame
//    Boson copy1 = *this;
//    copy1.unboost(rap_boost);
//    copy1.unboost(pt_boost);
//    // get a 2nd copy where we have mirrored the z components
//    Boson copy2(
//      fastjet::PseudoJet(-copy1.p1.px(), -copy1.p1.py(), +copy1.p1.pz(), copy1.p1.E()),
//      fastjet::PseudoJet(-copy1.p2.px(), -copy1.p2.py(), +copy1.p2.pz(), copy1.p2.E())
//    );
//    // boost it back to the proper pt and rapidity
//    copy2.boost(pt_boost);
//    copy2.boost(rap_boost);
//    return copy2;
//  }
//

  // enquiry functions
  /// returns the sin^2(theta) where theta is the decay
  /// angle in the boson rest frame 
  /// (check Collins-Soper angle correspondence)
  double decay_sinsq_theta() const {
    fastjet::PseudoJet Delta_t = fastjet::PseudoJet(p1.px()-p2.px(), p1.py()-p2.py(), 0.0, 0.0);
    fastjet::PseudoJet pinv = Delta_t;
    if (pt() != 0) {
      //pinv -= b * (dot_product(b,Delta_t)/pt2() * (m()/mt()-1));
      pinv -= (*this) * (dot_product(*this,Delta_t)/pt2() * (sqrt(m2()/mt2())-1));
    }
    return std::min(1.0, pinv.pt2() / m2());
  }

  /// returns sin_theta in the range 0 .. 1
  double decay_sin_theta() const {
    return sqrt(decay_sinsq_theta());
  }
  /// returns sin_theta in the range 0 .. pi/2
  double decay_theta() const {
    return asin(sqrt(decay_sinsq_theta()));
  }

  double decay_cos_theta() const {
    return sqrt(1.0 - decay_sinsq_theta());
  }

  double pt_CS() const {
    return 0.5 * sqrt(m2() * decay_sinsq_theta());
  }

  /// returns the phi value in the range 0 .. pi/2.
  /// For non-zero pt of the boson; phi=0 means that 
  /// the decay is in the plane formed by the beam
  /// and the boson direction
  double decay_phi_0topi2() const {
    const Boson & b = *this;
    // construct a transverse vector parallel to the boson transverse
    // direction and with pt=1
    if (b.pt2() == 0) {
      double this_phi = std::abs(atan2(std::abs(p1.py()), std::abs(p1.px())));
      if (this_phi <= fastjet::pi/2) return this_phi;
      else                           return fastjet::pi - this_phi;
      
    }
    fastjet::PseudoJet boson_perp_norm(b.px(), b.py(), 0, 0);
    boson_perp_norm /= b.pt();
    double  sintheta_cosphi = dot_product(b.p1 - b.p2, boson_perp_norm) / b.mt();
    return acos(std::min(1.0, std::abs(sintheta_cosphi) / decay_sin_theta()));
  }

  // same as decay_phi_0topi2() (retained for legacy reasons)
  double decay_phi() const {return decay_phi_0topi2();}


  /// returns the decay angle phi in the range 0 .. pi.
  ///
  /// The distinction between the range 0..pi/2 and pi/2..pi is only
  /// well defined when the boson has non-zero pt and the pi/2...pi
  /// range means that the particle at more positive rapidity has the
  /// lower pt
  double decay_phi_0topi() const  {
    double this_decay_phi = decay_phi_0topi2();
    if ( (p1.pt2() > p2.pt2()) != (p1.rap() > p2.rap()) ) return fastjet::pi - this_decay_phi;
    else return this_decay_phi;
  }

  fastjet::PseudoJet p1,p2;

protected:
  void init(double mass, double theta, double phi) {
    double m2=mass/2;
    p1.reset_momentum(m2 * cos(phi)*sin(theta),
                      m2 * sin(phi)*sin(theta),
                      m2 * cos(theta),
                      m2);
    p2.reset_momentum(-p1.px(), -p1.py(), -p1.pz(), p1.E());
    reset_momentum(sum());
  }
};

std::ostream & operator<<(std::ostream & ostr, const fastjet::PseudoJet & p) {
  ostr << p.px() << " "
       << p.py() << " "
       << p.pz() << " "
       << p.E() << " "
       << "pt=" << p.pt() << " "
       << "m="  << p.m() << " "
       << "y="  << p.rap() << " ";
  return ostr;
}
std::ostream & operator<<(std::ostream & ostr, const Boson & b) {
  ostr << "boson: " << b() 
       << " ptCS=" << b.pt_CS()
       << " cos(theta)=" << b.decay_cos_theta()
       << " phi/pi=" << b.decay_phi_0topi()/fastjet::pi
       << "\n"
       << "   p1: " << b.p1 << "\n"
       << "   p2: " << b.p2 << std::endl;
  return ostr;
}


#endif // __BOSON_HH__