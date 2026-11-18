// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Three-dimensional periodic Anderson model on a cubic lattice.

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_PERIODIC_ANDERSON_MODEL_3D_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_PERIODIC_ANDERSON_MODEL_3D_HPP

#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"

namespace dca {
namespace phys {
namespace models {

// The band indices 0 and 1 denote the conduction (c) and localized (f) orbitals,
// respectively.  The Hamiltonian is
//
//   H_0(k) = [ -2t(cos(k_x) + cos(k_y) + cos(k_z))   V         ]
//            [ V                                      epsilon_f ],
//
// with an on-site interaction U n_{f,up} n_{f,down} on the f orbital.
template <typename PointGroup>
class PeriodicAndersonModel3D {
public:
  static constexpr bool complex_g0 = false;
  static constexpr bool spin_symmetric = true;

  using LDA_point_group = domains::no_symmetry<3>;
  using DCA_point_group = PointGroup;

  static constexpr ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  static constexpr ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  static constexpr int DIMENSION = 3;
  static constexpr int BANDS = 2;
  static constexpr int CONDUCTION_BAND = 0;
  static constexpr int LOCALIZED_BAND = 1;

  static constexpr int transformationSignOfR(int, int, int) {
    return 1;
  }
  static constexpr int transformationSignOfK(int, int, int) {
    return 1;
  }

  static const double* initializeRDCABasis();
  static const double* initializeKDCABasis();
  static const double* initializeRLDABasis();
  static const double* initializeKLDABasis();

  static std::vector<int> flavors();
  static std::vector<std::vector<double>> aVectors();
  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  template <typename BandDmn, typename SpinDmn, typename RDmn, typename Parameters>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& interaction,
      const Parameters& parameters);

  template <class Domain>
  static void initializeHSymmetry(func::function<int, Domain>& symmetries);

  template <typename Parameters, typename Scalar, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0(
      const Parameters& parameters,
      func::function<Scalar, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);
};

template <typename PointGroup>
const double* PeriodicAndersonModel3D<PointGroup>::initializeRDCABasis() {
  static const std::array<double, 9> basis{1., 0., 0., 0., 1., 0., 0., 0., 1.};
  return basis.data();
}

template <typename PointGroup>
const double* PeriodicAndersonModel3D<PointGroup>::initializeKDCABasis() {
  static const std::array<double, 9> basis{2. * M_PI, 0., 0., 0., 2. * M_PI, 0., 0., 0., 2. * M_PI};
  return basis.data();
}

template <typename PointGroup>
const double* PeriodicAndersonModel3D<PointGroup>::initializeRLDABasis() {
  return initializeRDCABasis();
}

template <typename PointGroup>
const double* PeriodicAndersonModel3D<PointGroup>::initializeKLDABasis() {
  return initializeKDCABasis();
}

template <typename PointGroup>
std::vector<int> PeriodicAndersonModel3D<PointGroup>::flavors() {
  return {CONDUCTION_BAND, LOCALIZED_BAND};
}

template <typename PointGroup>
std::vector<std::vector<double>> PeriodicAndersonModel3D<PointGroup>::aVectors() {
  return {{0., 0., 0.}, {0., 0., 0.}};
}

template <typename PointGroup>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> PeriodicAndersonModel3D<
    PointGroup>::orbitalPermutations() {
  return {};
}

template <typename PointGroup>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename Parameters>
void PeriodicAndersonModel3D<PointGroup>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& interaction,
    const Parameters& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("The three-dimensional periodic Anderson model has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  interaction = 0.;
  const int origin = RDmn::parameter_type::origin_index();
  const double U = parameters.get_U();

  interaction(LOCALIZED_BAND, 0, LOCALIZED_BAND, 1, origin) = U;
  interaction(LOCALIZED_BAND, 1, LOCALIZED_BAND, 0, origin) = U;
}

template <typename PointGroup>
template <class Domain>
void PeriodicAndersonModel3D<PointGroup>::initializeHSymmetry(func::function<int, Domain>& symmetries) {
  symmetries = -1;

  for (int spin = 0; spin < 2; ++spin) {
    symmetries(CONDUCTION_BAND, spin, CONDUCTION_BAND, spin) = CONDUCTION_BAND;
    symmetries(LOCALIZED_BAND, spin, LOCALIZED_BAND, spin) = LOCALIZED_BAND;
  }
}

template <typename PointGroup>
template <typename Parameters, typename Scalar, typename BandDmn, typename SpinDmn, typename KDmn>
void PeriodicAndersonModel3D<PointGroup>::initializeH0(
    const Parameters& parameters,
    func::function<Scalar, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("The three-dimensional periodic Anderson model has two bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  H_0 = Scalar(0);

  const double t = parameters.get_t();
  const double epsilon_f = parameters.get_epsilon_f();
  const double hybridization = parameters.get_V();
  const auto& momenta = KDmn::get_elements();

  for (int k_ind = 0; k_ind < KDmn::dmn_size(); ++k_ind) {
    const auto& k = momenta[k_ind];
    const double epsilon_k = -2. * t * (std::cos(k[0]) + std::cos(k[1]) + std::cos(k[2]));

    for (int spin = 0; spin < SpinDmn::dmn_size(); ++spin) {
      H_0(CONDUCTION_BAND, spin, CONDUCTION_BAND, spin, k_ind) = epsilon_k;
      H_0(LOCALIZED_BAND, spin, LOCALIZED_BAND, spin, k_ind) = epsilon_f;
      H_0(CONDUCTION_BAND, spin, LOCALIZED_BAND, spin, k_ind) = hybridization;
      H_0(LOCALIZED_BAND, spin, CONDUCTION_BAND, spin, k_ind) = hybridization;
    }
  }
}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_PERIODIC_ANDERSON_MODEL_3D_HPP
