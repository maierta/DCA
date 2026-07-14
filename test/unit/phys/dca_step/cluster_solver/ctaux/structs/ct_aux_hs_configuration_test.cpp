// Copyright (C) 2025 UT-Battelle, LLC
// All rights reserved.
//
// See LICENSE for terms of usage.
// See CITATION.md for citation guidelines, if DCA++ is used for scientific publications.
//
// Unit tests for CT_AUX_HS_configuration, focusing on get_first_non_interacting_spin_index.

#include "dca/phys/dca_step/cluster_solver/ctaux/structs/ct_aux_hs_configuration.hpp"

#include <vector>

#include "dca/testing/gtest_h_w_warning_blocking.h"

#include "dca/math/random/std_random_wrapper.hpp"
#include "dca/phys/dca_step/cluster_solver/ctaux/structs/cv.hpp"

using Scalar = double;

#include "test/mock_mcconfig.hpp"
namespace dca {
namespace config {
using McOptions = MockMcOptions<Scalar>;
}  // namespace config
}  // namespace dca

#include "test/unit/phys/dca_step/cluster_solver/test_setup.hpp"

constexpr char input_name[] =
    DCA_SOURCE_DIR "/test/unit/phys/dca_step/cluster_solver/ctaux/structs/input.json";

using CtauxHsConfigurationTest =
    dca::testing::G0Setup<Scalar, dca::testing::LatticeBilayer, dca::ClusterSolverId::CT_AUX,
                          input_name>;
using Parameters = CtauxHsConfigurationTest::Parameters;

TEST_F(CtauxHsConfigurationTest, EmptyConfig) {
  dca::phys::solver::ctaux::CV<Parameters>::get_H_interaction() = data_->H_interactions;

  std::vector<double> random(200);
  for (auto& x : random)
    x = static_cast<double>(std::rand()) / RAND_MAX;
  Parameters::random_number_generator rng(random);
  dca::phys::solver::ctaux::CT_AUX_HS_configuration<Parameters> config(parameters_, rng);

  // Empty configuration: method returns 0 (same as size() since both are 0).
  EXPECT_EQ(0, config.get_first_non_interacting_spin_index(dca::phys::e_UP));
  EXPECT_EQ(0, config.get_first_non_interacting_spin_index(dca::phys::e_DN));
}

TEST_F(CtauxHsConfigurationTest, AllAnnihilatableReturnsSize) {
  dca::phys::solver::ctaux::CV<Parameters>::get_H_interaction() = data_->H_interactions;

  std::vector<double> random(200);
  for (auto& x : random)
    x = static_cast<double>(std::rand()) / RAND_MAX;
  Parameters::random_number_generator rng(random);
  dca::phys::solver::ctaux::CT_AUX_HS_configuration<Parameters> config(parameters_, rng);
  config.initialize();

  // After initialize(), all vertices are interacting (annihilatable).
  // Method should return size() as the "not found" sentinel.
  const int up_size = config.get(dca::phys::e_UP).size();
  const int dn_size = config.get(dca::phys::e_DN).size();

  EXPECT_EQ(up_size, config.get_first_non_interacting_spin_index(dca::phys::e_UP));
  EXPECT_EQ(dn_size, config.get_first_non_interacting_spin_index(dca::phys::e_DN));
}

TEST_F(CtauxHsConfigurationTest, FirstEntryNonAnnihilatable) {
  dca::phys::solver::ctaux::CV<Parameters>::get_H_interaction() = data_->H_interactions;

  std::vector<double> random(200);
  for (auto& x : random)
    x = static_cast<double>(std::rand()) / RAND_MAX;
  Parameters::random_number_generator rng(random);
  dca::phys::solver::ctaux::CT_AUX_HS_configuration<Parameters> config(parameters_, rng);
  config.initialize();

  // Mark the vertex referenced by configuration_e_UP[0] as non-annihilatable.
  const auto& up = config.get(dca::phys::e_UP);
  if (!up.empty()) {
    int target_idx = up[0].get_configuration_index();
    config[target_idx].set_annihilatable(false);

    // The first entry is now non-annihilatable, so the method should return 0.
    EXPECT_EQ(0, config.get_first_non_interacting_spin_index(dca::phys::e_UP));
    EXPECT_FALSE(config[target_idx].is_annihilatable());
  }

  // Mark the vertex referenced by configuration_e_DN[0] as non-annihilatable.
  const auto& dn = config.get(dca::phys::e_DN);
  if (!dn.empty()) {
    int target_idx = dn[0].get_configuration_index();
    config[target_idx].set_annihilatable(false);

    EXPECT_EQ(0, config.get_first_non_interacting_spin_index(dca::phys::e_DN));
    EXPECT_FALSE(config[target_idx].is_annihilatable());
  }
}

TEST_F(CtauxHsConfigurationTest, InsertNoninteractingVertex) {
  dca::phys::solver::ctaux::CV<Parameters>::get_H_interaction() = data_->H_interactions;

  std::vector<double> random(200);
  for (auto& x : random)
    x = static_cast<double>(std::rand()) / RAND_MAX;
  Parameters::random_number_generator rng(random);
  dca::phys::solver::ctaux::CT_AUX_HS_configuration<Parameters> config(parameters_, rng);
  config.initialize();

  int old_size = config.size();

  // Insert a non-interacting vertex (mark_annihilatable = false).
  config.insert_random_noninteracting_vertex(false);

  EXPECT_EQ(old_size + 1, config.size());
  EXPECT_FALSE(config[old_size].is_annihilatable());

  // Verify e_UP: the first non-annihilatable index should point to a non-annihilatable vertex,
  // and all entries before it must be annihilatable.
  const auto& up = config.get(dca::phys::e_UP);
  int first_up = config.get_first_non_interacting_spin_index(dca::phys::e_UP);
  if (first_up < static_cast<int>(up.size())) {
    int config_idx = up[first_up].get_configuration_index();
    EXPECT_FALSE(config[config_idx].is_annihilatable());
    for (int i = 0; i < first_up; ++i) {
      EXPECT_TRUE(config[up[i].get_configuration_index()].is_annihilatable());
    }
  }

  // Same check for e_DN.
  const auto& dn = config.get(dca::phys::e_DN);
  int first_dn = config.get_first_non_interacting_spin_index(dca::phys::e_DN);
  if (first_dn < static_cast<int>(dn.size())) {
    int config_idx = dn[first_dn].get_configuration_index();
    EXPECT_FALSE(config[config_idx].is_annihilatable());
    for (int i = 0; i < first_dn; ++i) {
      EXPECT_TRUE(config[dn[i].get_configuration_index()].is_annihilatable());
    }
  }
}
