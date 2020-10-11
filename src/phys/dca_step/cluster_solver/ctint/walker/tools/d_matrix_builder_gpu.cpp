// Copyright (C) 2018 ETH Zurich
// Copyright (C) 2018 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Authors: Giovanni Balduzzi (gbalduzz@itp.phys.ethz.ch)
//          Jérémie Bouquet (bouquetj@gmail.com).
//
// GPU implementation of d_matrix_builder.hpp

#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/d_matrix_builder_gpu.hpp"

#include <complex>

#include "dca/linalg/matrix_view.hpp"
#include "dca/phys/dca_step/cluster_solver/ctint/walker/tools/kernels_interface.hpp"

namespace dca {
namespace phys {
namespace solver {
namespace ctint {
// dca::phys::solver::ctint::

using namespace dca::linalg;

template <typename Real>
DMatrixBuilder<linalg::GPU, Real>::DMatrixBuilder(const G0Interpolation<GPU, Real>& g0,
                                                  const linalg::Matrix<int, linalg::CPU>& site_add,
                                                  const linalg::Matrix<int, linalg::CPU>& site_diff,
                                                  const int nb, const int r0)
    : BaseClass(g0, site_diff, nb), g0_ref_(g0) {
  assert(site_add.size() == site_diff.size());
  SolverHelper::set(site_add.ptr(), site_add.leadingDimension(), site_diff.ptr(),
                    site_diff.leadingDimension(), nb, site_add.nrRows(), r0);
}

template <typename Scalar>
void DMatrixBuilder<GPU, Scalar>::computeG0(Matrix& G0,
                                            const details::DeviceConfiguration& configuration,
                                            const int n_init, bool right_section,
                                            cudaStream_t stream) const {
  if (G0.nrRows() * G0.nrCols() == 0)
    return;

  dca::linalg::MatrixView<Scalar, linalg::GPU> g0_view(G0);
  details::buildG0Matrix(g0_view, n_init, right_section, configuration, g0_ref_, stream);
}

// Instantation.
template class DMatrixBuilder<GPU, float>;
template class DMatrixBuilder<GPU, double>;
template class DMatrixBuilder<GPU, std::complex<float>>;
template class DMatrixBuilder<GPU, std::complex<double>>;

}  // namespace ctint
}  // namespace solver
}  // namespace phys
}  // namespace dca
