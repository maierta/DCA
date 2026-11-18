// Copyright (C) 2026 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.

#include "dca/platform/dca_gpu.h"
#include "dca/phys/models/analytic_hamiltonians/periodic_anderson_model_3d.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "dca/function/domains.hpp"
#include "dca/function/function.hpp"
#include "dca/io/json/json_reader.hpp"
#include "dca/phys/domains/cluster/cluster_domain_aliases.hpp"
#include "dca/phys/domains/cluster/cluster_domain_initializer.hpp"
#include "dca/phys/domains/cluster/symmetries/point_groups/no_symmetry.hpp"
#include "dca/phys/models/tight_binding_model.hpp"
#include "dca/phys/parameters/model_parameters.hpp"
#include "dca/testing/gtest_h_w_warning_blocking.h"

namespace {

using PointGroup = dca::phys::domains::no_symmetry<3>;
using Lattice = dca::phys::models::PeriodicAndersonModel3D<PointGroup>;
using Model = dca::phys::models::TightBindingModel<Lattice>;
using Parameters = dca::phys::params::ModelParameters<Model>;

using BandDmn = dca::func::dmn_0<dca::func::dmn<2, int>>;
using SpinDmn = dca::func::dmn_0<dca::func::dmn<2, int>>;
using BandSpinDmn = dca::func::dmn_variadic<BandDmn, SpinDmn>;

}  // namespace

TEST(PeriodicAndersonModel3DTest, BasisAndModelParameters) {
  EXPECT_EQ(3, Lattice::DIMENSION);
  EXPECT_EQ(2, Lattice::BANDS);

  const double* real_basis = Lattice::initializeRDCABasis();
  const double* momentum_basis = Lattice::initializeKDCABasis();
  for (int row = 0; row < Lattice::DIMENSION; ++row)
    for (int column = 0; column < Lattice::DIMENSION; ++column) {
      const double identity = row == column ? 1. : 0.;
      EXPECT_DOUBLE_EQ(identity, real_basis[row * Lattice::DIMENSION + column]);
      EXPECT_DOUBLE_EQ(2. * M_PI * identity, momentum_basis[row * Lattice::DIMENSION + column]);
    }

  dca::io::JSONReader reader;
  Parameters parameters;
  reader.open_file(
      DCA_SOURCE_DIR
      "/test/unit/phys/models/analytic_hamiltonians/periodic_anderson_model_input.json");
  parameters.readWrite(reader);
  reader.close_file();

  EXPECT_DOUBLE_EQ(1.25, parameters.get_t());
  EXPECT_DOUBLE_EQ(-2.5, parameters.get_epsilon_f());
  EXPECT_DOUBLE_EQ(0.75, parameters.get_V());
  EXPECT_DOUBLE_EQ(6., parameters.get_U());
}

TEST(PeriodicAndersonModel3DTest, InitializeH0) {
  using KDmn = dca::func::dmn_0<dca::func::dmn<4, std::vector<double>>>;
  KDmn::parameter_type::set_elements(
      {{0., 0., 0.}, {M_PI, 0., 0.}, {M_PI, M_PI, 0.}, {M_PI, M_PI, M_PI}});

  dca::func::function<std::complex<double>, dca::func::dmn_variadic<BandSpinDmn, BandSpinDmn, KDmn>> H_0;

  Parameters parameters;
  parameters.set_t(1.);
  parameters.set_epsilon_f(-2.);
  parameters.set_V(0.5);
  Lattice::initializeH0(parameters, H_0);

  const double conduction_energies[] = {-6., -2., 2., 6.};
  for (int k = 0; k < KDmn::dmn_size(); ++k) {
    for (int spin = 0; spin < SpinDmn::dmn_size(); ++spin) {
      EXPECT_DOUBLE_EQ(conduction_energies[k], H_0(0, spin, 0, spin, k).real());
      EXPECT_DOUBLE_EQ(-2., H_0(1, spin, 1, spin, k).real());
      EXPECT_DOUBLE_EQ(0.5, H_0(0, spin, 1, spin, k).real());
      EXPECT_DOUBLE_EQ(0.5, H_0(1, spin, 0, spin, k).real());
    }

    for (int band_1 = 0; band_1 < BandDmn::dmn_size(); ++band_1)
      for (int band_2 = 0; band_2 < BandDmn::dmn_size(); ++band_2) {
        EXPECT_DOUBLE_EQ(0., H_0(band_1, 0, band_2, 1, k).real());
        EXPECT_DOUBLE_EQ(0., H_0(band_1, 1, band_2, 0, k).real());
        for (int spin = 0; spin < SpinDmn::dmn_size(); ++spin)
          EXPECT_DOUBLE_EQ(0., H_0(band_1, spin, band_2, spin, k).imag());
      }
  }
}

TEST(PeriodicAndersonModel3DTest, InteractionActsOnlyOnLocalizedOrbital) {
  using CDA = dca::phys::ClusterDomainAliases<Lattice::DIMENSION>;
  using RClusterDmn = typename CDA::RClusterDmn;

  const std::vector<std::vector<int>> cluster{{2, 0, 0}, {0, 2, 0}, {0, 0, 2}};
  dca::phys::domains::cluster_domain_initializer<RClusterDmn>::execute(
      Lattice::initializeRDCABasis(), cluster);

  dca::func::function<double, dca::func::dmn_variadic<BandSpinDmn, BandSpinDmn, RClusterDmn>> interaction;
  Parameters parameters;
  parameters.set_U(6.);
  Lattice::initializeHInteraction(interaction, parameters);

  const int origin = RClusterDmn::parameter_type::origin_index();
  for (int r = 0; r < RClusterDmn::dmn_size(); ++r)
    for (int band_1 = 0; band_1 < BandDmn::dmn_size(); ++band_1)
      for (int spin_1 = 0; spin_1 < SpinDmn::dmn_size(); ++spin_1)
        for (int band_2 = 0; band_2 < BandDmn::dmn_size(); ++band_2)
          for (int spin_2 = 0; spin_2 < SpinDmn::dmn_size(); ++spin_2) {
            const bool is_f_interaction = r == origin && band_1 == Lattice::LOCALIZED_BAND &&
                                          band_2 == Lattice::LOCALIZED_BAND && spin_1 != spin_2;
            EXPECT_DOUBLE_EQ(is_f_interaction ? 6. : 0.,
                             interaction(band_1, spin_1, band_2, spin_2, r));
          }
}
