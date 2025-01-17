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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_CONTACTTIME_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_CONTACTTIME_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/observables/Observable.hpp"

#include "core/observables/LBVelocityProfile.hpp"
#include "core/observables/ContactTimes.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Observables {

template <typename CoreObs>
class ContactTimes
    : public AutoParameters<ContactTimes<CoreObs>, Observable> {
  using Base = AutoParameters<ContactTimes<CoreObs>, Observable>;

public:
  using Base::Base;
  
  ContactTimes() {
    this->add_parameters(
        {{"ids", AutoParameter::read_only,
          [this]() { return time_observable()->ids(); }},
          {"target_ids", AutoParameter::read_only,
          [this]() { return time_observable()->target_ids; }},
          {"contact_threshold", AutoParameter::read_only,
          [this]() { return time_observable()->contact_threshold; }}});
  }

  void do_construct(VariantMap const &params) override {
    ObjectHandle::context()->parallel_try_catch([&]() {
      m_observable =
          make_shared_from_args<CoreObs, std::vector<int>, std::vector<int>, double>(
              params, "ids", "target_ids", "contact_threshold");
    });
  }

  Variant do_call_method(const std::string &method, VariantMap const &parameters) override {
  if (method == "clean_contact_times") {
    time_observable()->clean_contact_times();
    return {};
  }
  if (method == "shape_last_contact_time") {
    auto const shape = time_observable()->shape_last_contact_time();
    return std::vector<int>{shape.begin(), shape.end()};
  }
  return Base::do_call_method(method, parameters);  // Call base class for unsupported methods
}


  std::shared_ptr<::Observables::ContactTimes> time_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

}  // namespace Observables
}  // namespace ScriptInterface



#endif
