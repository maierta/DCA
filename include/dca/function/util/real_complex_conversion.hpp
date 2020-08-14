// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Giovanni Balduzzi (gbalduzz@gitp.phys.ethz.ch)
//
// This file provides utility methods to convert functions from real to complex and vice versa.

#ifndef DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP
#define DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP

#include <complex>
#include <limits>
#include <stdexcept>

#include "dca/function/function.hpp"

#ifdef DCA_HAVE_MPI
#include <mpi.h>
#endif

namespace dca {
namespace func {
namespace util {
// dca::func::util::

// Returns a complex valued function whose real part is equal to f.
template <typename Scalartype, typename Dmn>
auto complex(const function<Scalartype, Dmn>& f) {
  function<std::complex<Scalartype>, Dmn> f_complex;

  for (int i = 0; i < f_complex.size(); ++i)
    f_complex(i) = std::complex<Scalartype>(f(i), 0);

  return f_complex;
}

// Returns a real valued function that is equal to the real part of f.
// If check_imaginary = true, the method checks whether the imaginary part of f is zero and throws a
// std::logic_error if this is not the case.
template <typename Scalartype, typename Dmn>
auto real(const function<std::complex<Scalartype>, Dmn>& f, const bool check_imaginary = false) {
  function<Scalartype, Dmn> f_real;
  Scalartype max_img = 0, rel_max = 0;

  for (int i = 0; i < f_real.size(); ++i) {
    const auto val = std::abs(f(i).imag());
    if (val > max_img) {
      max_img = val;
      rel_max = val / std::abs(f(i).real());
    }

    f_real(i) = f(i).real();
  }

  if (check_imaginary && max_img > 1e3 * std::numeric_limits<Scalartype>::epsilon()) {
    // throw(std::logic_error("The function is not purely real."));

    int id = 0;
#ifdef DCA_HAVE_MPI
    int initialized;
    MPI_Initialized(&initialized);
    if (initialized) {
      MPI_Comm_rank(MPI_COMM_WORLD, &id);
    }
#endif  // DCA_HAVE_MPI
    if (id == 0) {
      std::cerr << "WARNING: The function" << f.get_name() << " is not purely real. Max img "
                << max_img << " relative of max " << rel_max << std::endl;
    }
  }

  return f_real;
}

}  // namespace util
}  // namespace func
}  // namespace dca

#endif  // DCA_FUNCTION_UTIL_REAL_COMPLEX_CONVERSION_HPP
