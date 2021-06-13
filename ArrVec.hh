#include<vector>
#include<array>
#include<algorithm>

template<class T, unsigned int maxarrsize = 10>
class ArrVec {
public:

  typedef T value_type;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;

  ArrVec() : size_(0) {}

  ArrVec(size_type n) {
    size_ = n;
    if (n > maxarrsize) vector_.resize(n);
  }
  ArrVec(size_type n, const T& value) {
    size_ = n;
    if (n > maxarrsize) vector_.resize(n);
    T * loc = begin();
    T * endloc = loc + n;
    while (loc < endloc) {
      *loc = value;
      ++loc;
    }
  }
  const T & operator[](unsigned int i) const {
    if (vector_.size() == 0) return array_[i];
    else                     return vector_[i];
  }
  T & operator[](unsigned int i) {
    if (vector_.size() == 0) return array_[i];
    else                     return vector_[i];
  }

  T * begin() {return vector_.size()>0 ? &(vector_[0]) : &(array_[0]);}
  T * end()   {return begin() + size();}

  const T * begin() const {return vector_.size()>0 ? &(vector_[0]) : &(array_[0]);}
  const T * end()   const {return begin() + size();}

  void resize(size_type n) {
    size_type oldsize = size_;
    size_ = n;
    if (vector_.size() == 0) {
      if (n > maxarrsize) {
        vector_.resize(n);
        std::copy(array_.begin(), array_.begin()+oldsize, vector_.begin());
      }
    } else {
      vector_.resize(n);
    }    
  }

  void resize(size_type n, const T & v) {
    size_type oldsize = size_;
    size_ = n;
    if (vector_.size() == 0) {
      if (n > maxarrsize) {
        vector_.resize(n);
        std::copy(array_.begin(), array_.begin()+oldsize, vector_.begin());
        std::fill(vector_.begin()+oldsize, vector_.begin()+size_, v);
      } else {
        std::fill(array_.begin()+oldsize, array_.begin()+size_, v);
      }
    } else {
      vector_.resize(n, v);
    }    
  }


  void push_back(const T & value) {
    if (vector_.size() != 0) {
      vector_.push_back(value);
      ++size_;
    } else if (size_ < maxarrsize) {
      array_[size_] = value;
      ++size_;
    } else {
      resize(size_+1, value);
    }
  }

  const T & front() const {
    if (vector_.size() != 0) {return vector_.front();}
    else                     {return array_.front();}
  }

  T & front() {
    if (vector_.size() != 0) {return vector_.front();}
    else                     {return array_.front();}
  }

  const T & back() const {
    if (vector_.size() != 0) {return vector_.back();}
    else                     {return array_[size_-1];}
  }

  T & back() {
    if (vector_.size() != 0) {return vector_.back();}
    else                     {return array_[size_-1];}
  }


  size_t size() const {return size_;}

private:
  std::array<T,maxarrsize> array_;
  std::vector<T> vector_;
  size_t size_ = 0;
};
