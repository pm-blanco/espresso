/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

#define BOOST_TEST_MODULE "Lees-Edwards boundary conditions"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "lees_edwards/LeesEdwardsBC.hpp"
#include "lees_edwards/lees_edwards.hpp"

#include <utils/Vector.hpp>

#include <boost/range/algorithm/equal.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

using namespace LeesEdwards;

auto constexpr eps = std::numeric_limits<double>::epsilon();
auto constexpr tol = 100. * eps;

BOOST_AUTO_TEST_CASE(test_shear_direction) {
  LeesEdwardsBC le;
  le.shear_direction = 1;
  BoxGeometry box;
  box.set_lees_edwards_bc(le);
  auto const expected_direction = Utils::Vector3d{{0, 1, 0}};
  BOOST_CHECK_SMALL((shear_direction(box) - expected_direction).norm(), eps);
}

BOOST_AUTO_TEST_CASE(test_update_offset) {
  Particle p;
  p.image_box() = {2, 4, -1};
  p.lees_edwards_offset() = 1.5;

  LeesEdwardsBC le;
  le.shear_direction = 1;
  le.shear_plane_normal = 2;
  le.shear_velocity = 3.5;
  BoxGeometry box;
  box.set_lees_edwards_bc(le);
  UpdateOffset(box, 3.5)(p);
  auto const expected_offset = 1.5 - le.shear_velocity * 0.5 * 3.5 *
                                         (p.image_box()[le.shear_plane_normal]);
  BOOST_CHECK_CLOSE(p.lees_edwards_offset(), expected_offset, tol);
}

BOOST_AUTO_TEST_CASE(test_push) {
  auto const old_offset = -1.2;
  auto const shear_l = 6.;
  auto const shear_normal_l = 4.5;
  auto const dt = 0.8;
  auto const old_pos = Utils::Vector3d{{3., shear_normal_l * 1.1, 10.}};
  auto const old_vel = Utils::Vector3d{{-1.2, 2., 4.1}};

  Particle p;

  p.pos() = old_pos;
  p.v() = old_vel;

  p.image_box() = {2, 4, -1};
  p.lees_edwards_offset() = old_offset;

  LeesEdwardsBC le;
  le.shear_direction = 2;
  le.shear_plane_normal = 1;
  le.pos_offset = 2.5;
  le.shear_velocity = -3.1;

  BoxGeometry box;
  box.set_type(BoxType::LEES_EDWARDS);
  box.set_length({5., shear_normal_l, shear_l});
  box.set_lees_edwards_bc(le);

  // Test transition in one direction
  Push(box, dt)(p);
  auto expected_pos = old_pos - shear_direction(box) * le.pos_offset;
  auto expected_vel = old_vel - shear_direction(box) * le.shear_velocity;
  auto expected_offset =
      old_offset + le.pos_offset -
      le.shear_velocity * 0.5 * dt * p.image_box()[le.shear_plane_normal];
  BOOST_CHECK_SMALL((p.pos() - expected_pos).norm(), eps);
  BOOST_CHECK_SMALL((p.v() - expected_vel).norm(), eps);
  BOOST_CHECK_CLOSE(p.lees_edwards_offset(), expected_offset, tol);

  // Test transition in the other direction
  p.pos()[le.shear_plane_normal] = -1;
  Push(box, dt)(p);
  expected_pos = {old_pos[0], -1., old_pos[2]};
  expected_vel = old_vel;
  expected_offset = old_offset - le.shear_velocity * dt *
                                     p.image_box()[le.shear_plane_normal];
  BOOST_CHECK_SMALL((p.pos() - expected_pos).norm(), eps);
  BOOST_CHECK_SMALL((p.v() - expected_vel).norm(), eps);
  BOOST_CHECK_CLOSE(p.lees_edwards_offset(), expected_offset, tol);
}

BOOST_AUTO_TEST_CASE(protocol_off) {
  auto off = Off();
  BOOST_CHECK_EQUAL(get_pos_offset(17.3, off), 0.0);
  BOOST_CHECK_EQUAL(get_shear_velocity(17.3, off), 0.0);
}

BOOST_AUTO_TEST_CASE(protocol_lin) {
  auto const t0 = 1.2;
  auto const x0 = -2.1;
  auto const v = 2.6;
  auto linear = LinearShear(x0, v, t0);
  BOOST_CHECK_CLOSE(get_pos_offset(3.3, linear), x0 + v * (3.3 - t0), tol);
  BOOST_CHECK_CLOSE(get_shear_velocity(17.3, linear), v, tol);
}

BOOST_AUTO_TEST_CASE(protocol_osc) {
  auto const t0 = 1.2;
  auto const a = 3.1;
  auto const o = 2.1;
  auto osc = OscillatoryShear(a, o, t0);
  BOOST_CHECK_CLOSE(get_pos_offset(3.3, osc), a * sin(o * (3.3 - t0)), tol);
  BOOST_CHECK_CLOSE(get_shear_velocity(3.3, osc), a * o * cos(o * (3.3 - t0)),
                    tol);
}
