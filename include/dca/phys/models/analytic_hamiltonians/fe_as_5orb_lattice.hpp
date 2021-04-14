// Copyright (C) 2021 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE.txt for terms of usage.
//  See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Author: Thomas Maier (maierta@ornl.gov)
//
// 5-orbital FeAs lattice (see Graser et al., NJP 11, 25016 (2009).

#ifndef DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FE_AS_5ORB_HPP
#define DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_FE_AS_5ORB_HPP

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/2d/2d_square.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/analytic_hamiltonians/cluster_shape_type.hpp"
#include "dca/util/type_list.hpp"

namespace dca {
namespace phys {
namespace models {
// dca::phys::models::


template <typename PointGroup>
class FeAs5Orb {
public:
  typedef domains::no_symmetry<2> LDA_point_group;
  typedef PointGroup DCA_point_group;

  const static ClusterShapeType DCA_cluster_shape = BETT_CLUSTER;
  const static ClusterShapeType LDA_cluster_shape = PARALLELEPIPED;

  const static int DIMENSION = 2;
  const static int BANDS = 5;

  // Rotations of pi/2 are an anti-symmetry on the band off-diagonal.
  static int transformationSign(int b1, int b2, int s);

  static int transformationSignOfR(int b1, int b2, int s) {
    return transformationSign(b1, b2, s);
  }
  static int transformationSignOfK(int b1, int b2, int s) {
    return transformationSign(b1, b2, s);
  }

  static double* initializeRDCABasis();

  static double* initializeKDCABasis();

  static double* initializeRLDABasis();

  static double* initializeKLDABasis();

  static std::vector<int> flavors();

  static std::vector<std::vector<double>> aVectors();

  static std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> orbitalPermutations();

  // Initializes the interaction Hamiltonian in real space.
  template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
  static void initializeHInteraction(
      func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
      const parameters_type& parameters);

  template <class HType, class Parameters>
  static void initializeNonDensityInteraction(HType& interaction, const Parameters& pars);

  template <class domain>
  static void initializeHSymmetry(func::function<int, domain>& H_symmetry);

  // Initializes the tight-binding (non-interacting) part of the momentum space Hamiltonian.
  // Preconditions: The elements of KDmn are two-dimensional (access through index 0 and 1).
  template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
  static void initializeH0(
      const ParametersType& parameters,
      func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                    func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0);
};

template <typename PointGroup>
int FeAs5Orb<PointGroup>::transformationSign(int b1, int b2, int s) {
  if (!std::is_same<PointGroup, domains::D4>::value)
    return 1;

  if (b1 == b2)
    return 1;

  //  using List = typename PointGroup::PointGroup_list;
  using dca::util::IndexOf;

  constexpr bool symmetrize_off_diagonal = false;

  if (symmetrize_off_diagonal) {
    //    const bool is_odd_rotation = IndexOf<domains::Cn_2D<1, 4>, List>::value == s ||
    //                                 IndexOf<domains::Cn_2D<3, 4>, List>::value == s;
    // TODO: generalize.
    const bool is_odd_rotation = (s % 2) == 1;

    return is_odd_rotation ? -1 : 1;
  }

  else {
    // TODO: generalize.
    const bool is_identity = s == 0;
    return is_identity ? 1 : 0;
  }
}

template <typename PointGroup>
double* FeAs5Orb<PointGroup>::initializeRDCABasis() {
  static std::array<double, 4> r_DCA{1, 0, 0, 1};
  return r_DCA.data();
}

template <typename PointGroup>
double* FeAs5Orb<PointGroup>::initializeKDCABasis() {
  static std::array<double, 4> k_DCA{2 * M_PI, 0, 0, 2 * M_PI};
  return k_DCA.data();
}

template <typename PointGroup>
double* FeAs5Orb<PointGroup>::initializeRLDABasis() {
  static std::array<double, 4> basis{1, 0, 0, 1};
  return basis.data();
}

template <typename PointGroup>
double* FeAs5Orb<PointGroup>::initializeKLDABasis() {
  static std::array<double, 4> basis{2 * M_PI, 0, 0, 2 * M_PI};
  return basis.data();
}

template <typename PointGroup>
std::vector<int> FeAs5Orb<PointGroup>::flavors() {
  return {0, 1, 2, 3, 4};
}

template <typename PointGroup>
std::vector<std::vector<double>> FeAs5Orb<PointGroup>::aVectors() {
  return {{0, 0}, {0, 0}, {0, 0}, {0, 0}, {0, 0}};
}

template <typename PointGroup>
std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> FeAs5Orb<PointGroup>::orbitalPermutations() {
  return {};
}

template <typename PointGroup>
template <typename BandDmn, typename SpinDmn, typename RDmn, typename parameters_type>
void FeAs5Orb<PointGroup>::initializeHInteraction(
    func::function<double, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                              func::dmn_variadic<BandDmn, SpinDmn>, RDmn>>& H_interaction,
    const parameters_type& parameters) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("5 orbital lattice has five bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const int origin = RDmn::parameter_type::origin_index();

  const double U = parameters.get_U();  // Intra-orbital interaction.
  const double V = parameters.get_V();  // Inter-orbital interaction.
  const double J = parameters.get_J();  // Spin-spin term.
  // const double V_prime = parameters.get_V_prime();  // Different orbital, same spin.

  H_interaction = 0.;

  // Note: the solvers expect a double counted Hamiltonian.
  for (int b1 = 0; b1 < BANDS; b1++) {
    for (int s1 = 0; s1 < 2; s1++) {
      for (int b2 = 0; b2 < BANDS; b2++) {
        for (int s2 = 0; s2 < 2; s2++) {
          if (b1 == b2 && s1 != s2)
            H_interaction(b1, s1, b2, s2, origin) = U;

          if (b1 != b2)
            H_interaction(b1, s1, b2, s2, origin) = V - (s1 == s2) * J;
        }
      }
    }
  }
}

template <typename PointGroup>
template <class HType, class Parameters>
void FeAs5Orb<PointGroup>::initializeNonDensityInteraction(HType& interaction,
                                                            const Parameters& pars) {
  const double J = pars.get_J();
  const double Jp = pars.get_Jp();

  constexpr int up(0), down(1);

  interaction = 0.;
  for (int b1 = 0; b1 < BANDS; b1++)
    for (int b2 = 0; b2 < BANDS; b2++) {
      if (b1 == b2)
        continue;
      interaction(b1, up, b2, up, b2, down, b1, down, 0) = J;
      interaction(b1, up, b2, up, b1, down, b2, down, 0) = Jp;
    }
}
template <typename PointGroup>
template <class domain>
void FeAs5Orb<PointGroup>::initializeHSymmetry(func::function<int, domain>& H_symmetries) {
  H_symmetries = -1;

  H_symmetries(0, 0, 0, 0) = 0;
  H_symmetries(0, 1, 0, 1) = 0;

  H_symmetries(1, 0, 1, 0) = 1;
  H_symmetries(1, 1, 1, 1) = 1;
}

template <typename PointGroup>
template <typename ParametersType, typename ScalarType, typename BandDmn, typename SpinDmn, typename KDmn>
void FeAs5Orb<PointGroup>::initializeH0(
    const ParametersType& parameters,
    func::function<ScalarType, func::dmn_variadic<func::dmn_variadic<BandDmn, SpinDmn>,
                                                  func::dmn_variadic<BandDmn, SpinDmn>, KDmn>>& H_0) {
  if (BandDmn::dmn_size() != BANDS)
    throw std::logic_error("FeAs5Orb lattice has 5 bands.");
  if (SpinDmn::dmn_size() != 2)
    throw std::logic_error("Spin domain size must be 2.");

  const auto& k_vecs = KDmn::get_elements();

  const auto t11x = -0.14;
  const auto t33x =  0.35;
  const auto t44x =  0.23;
  const auto t55x = -0.10;

  const auto t11y = -0.40;
  
  const auto t11xy =  0.28;
  const auto t33xy = -0.105;
  const auto t44xy =  0.15;
  
  const auto t11xx =  0.02;
  const auto t33xx = -0.02;
  const auto t44xx = -0.03;
  const auto t55xx = -0.04;
  
  const auto t11xxy = -0.035;
  const auto t44xxy = -0.03;
  const auto t55xxy =  0.02;

  const auto t11xyy =  0.005;
  
  const auto t11xxyy =  0.035;
  const auto t44xxyy = -0.03;
  const auto t55xxyy = -0.01;
  
  const auto t13x = -0.354;
  const auto t14x =  0.339;
  const auto t15x = -0.198;
  const auto t35x = -0.3;

  const auto t12xy =  0.05;
  const auto t13xy =  0.099;
  const auto t14xy =  0.014;
  const auto t15xy = -0.085;
  const auto t45xy = -0.15;

  const auto t12xxy = -0.015;
  const auto t13xxy =  0.021;
  const auto t14xxy =  0.028;
  const auto t34xxy = -0.01;
  const auto t35xxy = -0.02;
  
  const auto t12xxyy =  0.035;
  const auto t15xxyy = -0.014;
  const auto t45xxyy =  0.01;

  H_0 = 0;

  const ScalarType I(0, 1.);

  for (int k_id = 0; k_id < k_vecs.size(); ++k_id) {
    const auto& k = k_vecs[k_id];
    const auto cx = cos(k[0]);
    const auto cy = cos(k[1]);
    const auto c2x = cos(2*k[0]);
    const auto c2y = cos(2*k[1]);
    const auto sx = sin(k[0]);
    const auto sy = sin(k[1]);
    const auto s2x = sin(k[0]);
    const auto s2y = sin(k[1]);

    const auto xi11 = 2*t11x*cx + 2*t11y*cy + 4*t11xy*cx*cy
                    + 2*t11xx*(c2x-c2y) + 4*t11xxy*c2x*cy
                    + 4*t11xyy*c2y*cx + 4*t11xxyy*c2x*c2y;
    
    const auto xi22 = 2*t11y*cx + 2*t11x*cy + 4*t11xy*cx*cy
                    - 2*t11xx*(c2x-c2y) + 4*t11xyy*c2x*cy
                    + 4*t11xxy*c2y*cx + 4*t11xxyy*c2x*c2y;
    
    const auto xi33 = 2*t33x*(cx+cy) + 4*t33xy*cx*cy + 2*t33xx*(c2x+c2y);
    
    const auto xi44 = 2*t44x*(cx+cy) + 4*t44xy*cx*cy + 2*t44xx*(c2x+c2y)
                    + 4*t44xxy*(c2x*cy+c2y*cx) + 4*t44xxyy*c2x*c2y;
    
    const auto xi55 = 2*t55x*(cx+cy) + 2*t55xx*(c2x+c2y)   
                    + 4*t55xxy*(c2x*cy+c2y*cx) + 2*t55xxyy*c2x*c2y;

    const auto xi12 = -4*t12xy*sx*sy - 4*t12xxy*(s2x*sy+s2y*sx) - 4*t12xxyy*s2x*s2y;

    const auto xi13 = 2.*I*t13x*sy + 4.*I*t13xy*sy*cx
                    - 4.*I*t13xxy*(s2y*cx-c2x*sy);

    const auto xi23 = -2.*I*t13x*sx - 4.*I*t13xy*sx*cy;
                    + 4.*I*t13xxy*(s2x*cy-c2y*sx);

    const auto xi14 = 2.*I*t14x*sx + 4.*I*t14xy*cy*sx + 4.*I*t14xxy*s2x*cy;
    
    const auto xi24 = 2.*I*t14x*sy + 4.*I*t14xy*cx*sy + 4.*I*t14xxy*s2y*cx;

    const auto xi15 = 2.*I*t15x*sy - 4.*I*t15xy*sy*cx - 4.*I*t15xxyy*s2y*c2x;

    const auto xi25 = 2.*I*t15x*sx - 4.*I*t15xy*sx*cy - 4.*I*t15xxyy*s2x*c2y;
    
    const auto xi34 = 4*t34xxy*(s2y*sx-s2x*sy);

    const auto xi35 = 2*t35x*(cx-cy) + 4*t35xxy*(c2x*cy-c2y*cx);

    const auto xi45 = 4*t45xy*sx*sy + 4*t45xxyy*s2x*s2y;

    const auto e1 =  0.13;
    const auto e3 = -0.22;
    const auto e4 =  0.30;
    const auto e5 = -0.211;

    for (int s = 0; s < 2; ++s) {
      H_0(0, s, 0, s, k_id) = xi11 + e1;
      H_0(1, s, 1, s, k_id) = xi22 + e1;
      H_0(2, s, 2, s, k_id) = xi33 + e3;
      H_0(3, s, 3, s, k_id) = xi44 + e4;
      H_0(4, s, 4, s, k_id) = xi55 + e5;

      H_0(0, s, 1, s, k_id) = xi12;
      H_0(1, s, 0, s, k_id) = std::conj(xi12);

      H_0(0, s, 2, s, k_id) = xi13;
      H_0(2, s, 0, s, k_id) = std::conj(xi13);

      H_0(1, s, 2, s, k_id) = xi23;
      H_0(2, s, 1, s, k_id) = std::conj(xi23);

      H_0(0, s, 3, s, k_id) = xi14;
      H_0(3, s, 0, s, k_id) = std::conj(xi14);

      H_0(1, s, 3, s, k_id) = xi24;
      H_0(3, s, 1, s, k_id) = std::conj(xi24);

      H_0(0, s, 4, s, k_id) = xi15;
      H_0(4, s, 0, s, k_id) = std::conj(xi15);

      H_0(1, s, 4, s, k_id) = xi25;
      H_0(4, s, 1, s, k_id) = std::conj(xi25);

      H_0(2, s, 3, s, k_id) = xi34;
      H_0(3, s, 2, s, k_id) = std::conj(xi34);

      H_0(2, s, 4, s, k_id) = xi35;
      H_0(4, s, 2, s, k_id) = std::conj(xi35);

      H_0(3, s, 4, s, k_id) = xi45;
      H_0(4, s, 3, s, k_id) = std::conj(xi45);

    }
  }
}

}  // namespace models
}  // namespace phys
}  // namespace dca

#endif  // DCA_PHYS_MODELS_ANALYTIC_HAMILTONIANS_BILAYER_LATTICE_HPP
