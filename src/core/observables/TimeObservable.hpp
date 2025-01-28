/*
 * Copyright (C) 2025 The ESPResSo project
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
#ifndef OBSERVABLES_TIMEOBSERVABLE_HPP
#define OBSERVABLES_TIMEOBSERVABLE_HPP

#include "Observable.hpp"

#include <utils/math/make_lin_space.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Observables {

/** Observable tracking time evolution of contacts between `ids` and
 * `target_ids` within a given `contact_threshold` */
class TimeObservable : virtual public Observable {
public:
  double contact_threshold;
  std::vector<int> const target_ids;
  TimeObservable(double contact_threshold, std::vector<int> const &target_ids)
      : contact_threshold(contact_threshold), target_ids(target_ids) {}
};

} // Namespace Observables
#endif
