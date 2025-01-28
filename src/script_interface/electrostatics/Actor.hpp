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

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "script_interface/Context.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/get_value.hpp"
#include "script_interface/system/Leaf.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Coulomb {

/**
 * @brief Common interface for electrostatic actors.
 * Several methods are defined in initialize.cpp since they
 * depend on symbols only available in @ref coulomb.hpp, which cannot be
 * included in this header file for separation of concerns reasons.
 */
template <class SIClass, class CoreClass>
class Actor : public AutoParameters<Actor<SIClass, CoreClass>, System::Leaf> {
protected:
  using SIActorClass = SIClass;
  using CoreActorClass = CoreClass;
  using ObjectClass = Actor<SIClass, CoreClass>;
  using AutoParameters<ObjectClass, System::Leaf>::context;
  using AutoParameters<ObjectClass, System::Leaf>::add_parameters;
  using AutoParameters<ObjectClass, System::Leaf>::do_set_parameter;
  using System::Leaf::get_system;
  using System::Leaf::m_system;
  std::shared_ptr<CoreActorClass> m_actor;

  void on_bind_system(::System::System &) override {
    m_actor->bind_system(m_system.lock());
  }

private:
  friend SIClass;
  Actor();

public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  std::shared_ptr<CoreActorClass> actor() { return m_actor; }
  std::shared_ptr<CoreActorClass const> actor() const { return m_actor; }

protected:
  void set_charge_neutrality_tolerance(VariantMap const &params) {
    auto const key_chk = std::string("check_neutrality");
    auto const key_tol = std::string("charge_neutrality_tolerance");
    if (params.count(key_tol)) {
      do_set_parameter(key_tol, params.at(key_tol));
    }
    do_set_parameter(key_chk, params.at(key_chk));
  }
};

} // namespace Coulomb
} // namespace ScriptInterface

#endif // ELECTROSTATICS
