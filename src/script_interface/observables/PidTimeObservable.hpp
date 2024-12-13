/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PIDTIMEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PIDTIMEOBSERVABLE_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/observables/Observable.hpp"

#include "core/observables/DensityProfile.hpp"
#include "core/observables/FluxDensityProfile.hpp"
#include "core/observables/ForceDensityProfile.hpp"
#include "core/observables/PidTimeObservable.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace ScriptInterface {
namespace Observables {

template <typename CoreObs>
class PidTimeObservable
    : public AutoParameters<PidTimeObservable<CoreObs>, Observable> {
  using Base = AutoParameters<PidTimeObservable<CoreObs>, Observable>;

public:
  using Base::Base;
  PidTimeObservable() {
    this->add_parameters(
        {{"ids", AutoParameter::read_only,
          [this]() { return pid_time_observable()->ids(); }},
          {"target_ids", AutoParameter::read_only,
          [this]() { return pid_time_observable()->target_ids; }},
         {"contact_threshold", AutoParameter::read_only,
          [this]() { return pid_time_observable()->contact_threshold; }}});
  }

  void do_construct(VariantMap const &params) override {
    ObjectHandle::context()->parallel_try_catch([&]() {
      m_observable =
          make_shared_from_args<CoreObs, std::vector<int>, std::vector<int>, double>(
              params, "ids", "target_ids", "contact_threshold");
    });
  }

  
  std::shared_ptr<::Observables::PidTimeObservable>
  pid_time_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
