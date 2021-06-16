// This file is part of the two-body-cuts project
// Copyright [2021] [Gavin Salam and Emma Slade]
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

#ifndef __REALRANGE_HH__
#define __REALRANGE_HH__
#include <vector>
#include <cassert>
#include <iostream>
#include "ArrVec.hh"

namespace tbc {

// On MacOS (tested with 10.15.7 and system c++, Apple clang version
// 12.0.0) there is a significant speed-up (x2-3) to be had with the
// ArrVec class, which avoids the need for vector allocations for small
// sizes, instead using (part of) a fixed size array. 
//
// On linux it seems that the gain is much more modest, at the 10-15%
// level.

//template <class T> using lclvector = std::vector<T>; 
//template <class T> using lclvector = ArrVec<T,10>;
template <class T> using lclvector = ArrVec<T,6>;

/// class to hold a range, which means a (possibly non-continuous)
/// subset of the interval 0 to 1 (those values are represented as LO
/// and HI in the class).
///
/// One of the main purposes of a RealRange is that it implements
/// logical operations between RealRanges, making it easy to put
/// together a range of constraints as long as the can easily be mapped
/// to a common 0..1 region.
///
/// Internally, the range is held as an ordered vector of
/// non-overlapping RealRange:Segment objects, each of which contains a
/// continuous subset of the interval 0 to 1.
class RealRange {
public:
  class Segment {
  public:
    /// default constructor returns an empty segment
    Segment() : lo_(LO), hi_(LO) {}
    Segment (double lo_in, double hi_in) : lo_(lo_in), hi_(hi_in) {
      if (lo_ > hi_) std::swap(lo_, hi_);
      lo_ = std::max(LO, lo_);
      hi_ = std::min(HI, hi_);      
    }

    Segment & reset() {lo_ = LO; hi_ = LO; return *this;}
    Segment & reset(double lo, double hi) {
      lo_ = lo; hi_ = hi;
      return *this;
    }

    double lo() const {return lo_;}
    double hi() const {return hi_;}
    double extent() const {return hi() - lo();}
    double is_empty() const {return lo_ == hi_;}

    /// returns true iff there is no overlap between this
    /// and the other segment
    bool distinct(const Segment & other) const {
      return other.hi() < lo() || hi() < other.lo();
    }
    /// returns true iff this segment has some overlap
    /// with the other
    bool overlaps(const Segment & other) const {
      return !distinct(other);
    }

    /// returns true if the segment contains this value
    bool contains(double value) const {
      return value >= lo_ && value <= hi_;
    }

    /// merges the other segment with this one (i.e. logical or),
    /// assuming they are not distinct
    Segment & merge(const Segment & other) {
      assert(!distinct(other));
      lo_ = std::min(lo_, other.lo_);
      hi_ = std::max(hi_, other.hi_);
      return *this;
    }

    /// intersects the other segment with this one (i.e. logical and),
    /// assuming they are not distinct
    Segment & intersect(const Segment & other) {
      assert(!distinct(other));
      lo_ = std::max(lo_, other.lo_);
      hi_ = std::min(hi_, other.hi_);
      return *this;
    }

    bool operator==(const Segment & other) const {
      return (other.lo() == lo() && other.hi() == hi()) || (other.extent() == 0 && extent() == 0);
    }
    bool operator!=(const Segment & other) const {return ! ((*this)==other);}


  private:
    double lo_, hi_;
  public:
    constexpr static double LO = 0.0, HI = 1.0;
  };
  
  RealRange() {}
  RealRange(double lo, double hi) : segments_(1, Segment(lo, hi)) {}
  RealRange(const lclvector<Segment> & segments) : segments_(segments) {}
  RealRange(const Segment & segment) : segments_(1,segment) {}

  /// sets the range back to an empty range
  RealRange & reset() {segments_.resize(0); return *this;}

  /// returns the segments of which this range is composed.
  const lclvector<Segment> & segments() const {return segments_;}

  /// returns the total extent of the range, i.e. the sum of the extents 
  /// of the segments of which it is composed.
  double extent() const {
    double result = 0;
    for (const auto & s: segments()) result += s.extent();
    return result;
  }

  /// returns true if the range contains the specified value
  bool contains(double value) const {
    for (const auto & s: segments()) {
      if (s.contains(value)) return true;
    }
    return false;
  }

  /// returns the number of separate segments of which this range
  /// is composed.
  unsigned size() const {return segments().size();}

  /// returns true if the range is empty (i.e. not point is contained
  /// within it)
  bool is_empty() const {return segments().size() == 0;}

  /// performs (local) logical or of this range with the other
  RealRange & operator|=(const RealRange & other);
  /// performs (local) logical and of this range with the other
  RealRange & operator&=(const RealRange & other);

  /// returns true if this range is identical to the other
  bool operator==(const RealRange & other) const {
    if (size() != other.size()) return false;
    for (unsigned i = 0; i < size(); i++) {
      if (segments_[i] != other.segments_[i]) return false;
    }
    return true;
  }

  /// returns true if this range differs to the other
  bool operator!=(const RealRange & other) const {return !(*this == other);}

  ~RealRange() {}

  /// returns a veto range such that 
  /// (*this && veto).extent() - this->extent() == amount_to_veto
  /// The veto is contructed progressively from the top of the range.
  /// The amount_to_veto is set to max(0, amount_to_veto - this->extent()) 
  /// after the operation
  RealRange veto_amount_from_top(double & amount_to_veto) const {
    RealRange extra_veto;
    if (amount_to_veto < 0) return extra_veto;
    for (int i = int(size()) - 1; i >= 0; i--) {
      const auto & s = segments()[i];
      double extent = s.extent();
      if (extent > amount_to_veto) {
        extra_veto |= RealRange(s.hi() - amount_to_veto, s.hi());
        amount_to_veto = 0;
        break;
      } else {
        extra_veto |= RealRange(s);
        amount_to_veto -= s.extent();
      }
    }
    return extra_veto;
  }

  /// returns an extra range such that (*this && extra).extent() == amount_to_add
  /// the region is added progressively from the bottom of the range.
  /// The amount_to_add is set to max(0, amount_to_add - this->extent()) 
  /// after the operation
  RealRange extra_amount_from_bottom(double & amount_to_add) const {    
    RealRange extra_allowed;
    if (amount_to_add < 0) return extra_allowed;
    for (const auto & s: segments()) {
      double extent = s.extent();
      if (extent > amount_to_add) {
        extra_allowed |= RealRange(s.lo(), s.lo() + amount_to_add);
        amount_to_add = 0;
        break;
      } else {
        extra_allowed |= RealRange(s);
        amount_to_add -= s.extent();
      }
    }
    return extra_allowed;
  }

private:
  /// segments, to be held ordered
  lclvector<Segment> segments_;
};

std::ostream & operator<< (std::ostream & ostr, const RealRange::Segment & segment) {
  ostr << "[" << segment.lo() << "," << segment.hi() << "]";
  return ostr;
}
std::ostream & operator<< (std::ostream & ostr, const RealRange & frange) {
  for(const auto & segment: frange.segments() ) {
    ostr << segment;
  }
  ostr << "(tot=" << frange.extent() << ")";
  return ostr;
}


/// returns the negated range, e.g. [0.2,0.8] becomes [0.0,0.2][0.8,1.0]
RealRange operator!(const RealRange & other) {
  const auto & segments = other.segments();
  typedef RealRange::Segment Segment;
  // if our range is empty, then it's easy
  if (other.is_empty()) return RealRange(Segment::LO,Segment::HI);

  // otherwise construct things
  lclvector<Segment> new_segments;

  // special treatment for the lower end
  if (segments[0].lo() != Segment::LO) {
    new_segments.push_back(Segment(Segment::LO, segments[0].lo()));
  }

  // easy treatment for everything in between
  for (unsigned i = 0; i+1 < segments.size(); ++i) {
    new_segments.push_back(Segment(segments[i].hi(), segments[i+1].lo()));
  }

  // special treatment for the last one11
  if (segments.back().hi() != Segment::HI) {
    new_segments.push_back(Segment(segments.back().hi(), Segment::HI));
  }

  return RealRange(new_segments);
}

/// Returns the result of the logical or operation between ranges a and
/// b. E.g. if a=[0.0,0.2][0.4,0.6] and b=[0.3,0.5], the result is 
/// [0.0,0.2][0.3,0.6]
RealRange operator||(const RealRange & a, const RealRange & b) {
  typedef RealRange::Segment Segment;
  if (a.is_empty()) return b;
  if (b.is_empty()) return a;

  /// start off with an initial empty segment
  RealRange::Segment current_segment;
  lclvector<RealRange::Segment> new_segments;
  unsigned ia = 0, ib = 0;
  const Segment * seg_a, * seg_b;
  while (true) {    
    seg_a = ia < a.segments().size() ? &(a.segments()[ia]) : 0;
    seg_b = ib < b.segments().size() ? &(b.segments()[ib]) : 0;
    const Segment * next;
    if (!seg_a && !seg_b) break;
    if      (!seg_a) {next = seg_b; ib++;}
    else if (!seg_b) {next = seg_a; ia++;}
    else {
      if (seg_a->lo() < seg_b->lo()) {next = seg_a; ++ia;}
      else                           {next = seg_b; ++ib;} 
    }
    if (current_segment.is_empty()) {
      current_segment = *next;
    } else if (next->overlaps(current_segment)) {
      current_segment.merge(*next);
    } else {
      new_segments.push_back(current_segment);
      current_segment = *next;
    }
  }
  if (!current_segment.is_empty()) {
    new_segments.push_back(current_segment);
  }
  return RealRange(new_segments);
}


/// Returns the result of the logical and operation between ranges a and
/// b. E.g. if a=[0.0,0.2][0.4,0.6] and b=[0.3,0.5], the result is 
/// [0.4,0.5]
RealRange operator&&(const RealRange & a, const RealRange & b) {
  typedef RealRange::Segment Segment;
  if (a.is_empty() || b.is_empty()) return RealRange();

  /// start off with an initial empty segment
  RealRange::Segment current_a, current_b;
  lclvector<RealRange::Segment> new_segments;
  unsigned ia = 0, ib = 0;
  while (true) {
    if (current_a.is_empty() && ia < a.segments().size()) current_a = a.segments()[ia++];
    if (current_b.is_empty() && ib < b.segments().size()) current_b = b.segments()[ib++];

    if (current_a.is_empty() || current_b.is_empty()) break;

    //std::cout << ia << "," << ib << ": " << current_a << " " << current_b << std::endl;
    bool distinct = current_a.distinct(current_b);
    if (current_a.lo() < current_b.lo()) {
      if (distinct) {current_a.reset(); continue;}
      Segment new_seg = current_a;
      new_seg.intersect(current_b);
      new_segments.push_back(new_seg);
      if (current_b.hi() < current_a.hi()) {
        current_a.reset(current_b.hi(), current_a.hi());
        current_b.reset();
      } else {
        current_a.reset();
        current_b.reset(current_a.hi(), current_b.hi());
      }
    } else {
      if (distinct) {current_b.reset(); continue;}
      Segment new_seg = current_b;
      new_seg.intersect(current_a);
      new_segments.push_back(new_seg);
    }
    if (current_b.hi() < current_a.hi()) {
      current_a.reset(current_b.hi(), current_a.hi());
      current_b.reset();
    } else {
      current_a.reset();
      current_b.reset(current_a.hi(), current_b.hi());
    }

  }
  
  return RealRange(new_segments);
}



RealRange & RealRange::operator|=(const RealRange & other) {
  *this = (*this) || other;
  return *this;
}
RealRange & RealRange::operator&=(const RealRange & other) {
  *this = (*this) && other;
  return *this;
}

} // end of tbc namespace

#endif // __REALRANGE_HH__
