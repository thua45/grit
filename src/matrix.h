// A 2D matrix class with the functionality I need for PSSMs
/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#ifndef MCF_MATRIX_H
#define MCF_MATRIX_H

#include <algorithm>
#include <vector>
//#include "MCFgen.hpp"  // is_reverse

namespace mcf {
  // perl-style die
  void mcf_die(const std::string & message);

  // return true if the range equals its reverse
  template<class It> bool is_reverse(It start, It end);
  // normalize a range (make it sum to 1)
  // no guarantee that the answer will be exactly 1
  // the range had better contain floating-point values!
  template<class It> void normalize(It start, It end);

  template<class T> class matrix {
  private:
    unsigned nrows;
    unsigned ncols;
    std::vector<T> data;

  public:
    matrix()
      : nrows(0), ncols(0) {}

    matrix(unsigned r, unsigned c)
      : nrows(r), ncols(c), data(r*c) {}
    
    unsigned rows() const { return nrows; }
    
    unsigned cols() const { return ncols; }
    
    // assumes i+ncols is a valid iterator:
    template<class It> void push_row(It i)
    { data.insert(data.end(), i, i+ncols); ++nrows; }

    // Not sure why "typename" is needed:

    typename std::vector<T>::iterator operator[] (unsigned r)
    { return data.begin() + r * ncols; }
    //    T * operator[] (unsigned r)
    //    { return &data[0] + r * ncols; }

    typename std::vector<T>::const_iterator operator[] (unsigned r) const
    { return data.begin() + r * ncols; }
    //    const T * operator[] (unsigned r) const
    //    { return &data[0] + r * ncols; }

    void rotate180()  // same as reverse complement
    { std::reverse(data.begin(), data.end()); }

    bool is_rotate180() const  // is it palindromic?
    { return mcf::is_reverse(data.begin(), data.end()); }
  };
}

template<class It> void mcf::normalize(It start, It end)
{
  typename std::iterator_traits<It>::value_type tot =
    std::accumulate(start, end,  // added typename 16-2-2005:
		    typename std::iterator_traits<It>::value_type(0));

  assert(tot != 0);  // doesn't like being prefixed by std::

  tot = 1/tot;
  for (; start < end; ++start)
    *start *= tot;
}

inline void mcf::mcf_die(const std::string & message)
{
  std::cerr << message << std::endl;
  exit(EXIT_FAILURE);
}

#endif
