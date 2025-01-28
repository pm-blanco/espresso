/*
 * Copyright (C) 2022 The ESPResSo project
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

#pragma once

#include <boost/functional/hash.hpp>
#include <boost/variant.hpp>

#include <array>
#include <cstddef>

namespace BondBreakage {

// Delete Actions
struct DeleteBond {
  int particle_id;
  int bond_partner_id;
  int bond_type;
  std::size_t hash_value() const {
    std::size_t seed = 3875;
    boost::hash_combine(seed, particle_id);
    boost::hash_combine(seed, bond_partner_id);
    boost::hash_combine(seed, bond_type);
    return seed;
  }
  bool operator==(DeleteBond const &) const = default;
  bool operator!=(DeleteBond const &) const = default;
};

struct DeleteAngleBond {
  int particle_id;
  std::array<int, 2> bond_partner_id;
  int bond_type;
  std::size_t hash_value() const {
    std::size_t seed = 3876;
    boost::hash_combine(seed, particle_id);
    boost::hash_combine(seed, bond_partner_id);
    boost::hash_combine(seed, bond_type);
    return seed;
  }
  bool operator==(DeleteAngleBond const &) const = default;
  bool operator!=(DeleteAngleBond const &) const = default;
};

struct DeleteAllBonds {
  int particle_id_1;
  int particle_id_2;
  std::size_t hash_value() const {
    std::size_t seed = 75;
    boost::hash_combine(seed, particle_id_1);
    boost::hash_combine(seed, particle_id_2);
    return seed;
  }
  bool operator==(DeleteAllBonds const &) const = default;
  bool operator!=(DeleteAllBonds const &) const = default;
};

} // namespace BondBreakage

// Hash support for std::unordered_set
namespace boost {
template <> struct hash<BondBreakage::DeleteBond> {
  std::size_t operator()(BondBreakage::DeleteBond const &t) const noexcept {
    return t.hash_value();
  }
};
template <> struct hash<BondBreakage::DeleteAngleBond> {
  std::size_t
  operator()(BondBreakage::DeleteAngleBond const &t) const noexcept {
    return t.hash_value();
  };
};
template <> struct hash<BondBreakage::DeleteAllBonds> {
  std::size_t operator()(BondBreakage::DeleteAllBonds const &t) const noexcept {
    return t.hash_value();
  }
};
} // namespace boost
