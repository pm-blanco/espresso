/*
 * Copyright (C) 2020-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE "P3M utility functions"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "p3m/common.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <vector>

BOOST_AUTO_TEST_CASE(calc_meshift_false) {
  std::array<std::vector<int>, 3> const ref = {
      {std::vector<int>{0}, std::vector<int>{0, 1, -2, -1},
       std::vector<int>{0, 1, 2, 3, -3, -2, -1}}};

  auto const mesh = Utils::Vector3i{{1, 4, 7}};
  auto const val = calc_p3m_mesh_shift(mesh, false);

  for (std::size_t i = 0u; i < 3u; ++i) {
    for (std::size_t j = 0u; j < ref[i].size(); ++j) {
      BOOST_CHECK_EQUAL(val[i][j], ref[i][j]);
    }
  }
}

BOOST_AUTO_TEST_CASE(calc_meshift_true) {
  std::array<std::vector<int>, 3> const ref = {
      {std::vector<int>{0}, std::vector<int>{0, 1, 0, -1},
       std::vector<int>{0, 1, 2, 0, -3, -2, -1}}};

  auto const mesh = Utils::Vector3i{{1, 4, 7}};
  auto const val = calc_p3m_mesh_shift(mesh, true);

  for (std::size_t i = 0u; i < 3u; ++i) {
    for (std::size_t j = 0u; j < ref[i].size(); ++j) {
      BOOST_CHECK_EQUAL(val[i][j], ref[i][j]);
    }
  }
}
