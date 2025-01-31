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

#pragma once

#include "Constraint.hpp"

#include "core/BoxGeometry.hpp"
#include "core/cell_system/CellStructure.hpp"
#include "core/constraints/Constraint.hpp"
#include "core/constraints/ShapeBasedConstraint.hpp"
#include "core/system/System.hpp"

#include "script_interface/shapes/Shape.hpp"

#include <memory>
#include <string>

namespace ScriptInterface {
namespace Constraints {

class ShapeBasedConstraint : public Constraint {
  std::unique_ptr<VariantMap> m_params;
  std::weak_ptr<::System::System const> m_system;

public:
  ShapeBasedConstraint()
      : m_constraint(std::make_shared<::Constraints::ShapeBasedConstraint>()),
        m_shape(nullptr) {
    auto system = ::System::get_system().shared_from_this();
    m_constraint->bind_system(system);
    m_system = system;
    add_parameters({{"only_positive", m_constraint->only_positive()},
                    {"penetrable", m_constraint->penetrable()},
                    {"particle_type",
                     [this](Variant const &value) {
                       m_constraint->set_type(get_value<int>(value));
                     },
                     [this]() { return m_constraint->type(); }},
                    {"shape",
                     [this](Variant const &value) {
                       m_shape =
                           get_value<std::shared_ptr<Shapes::Shape>>(value);
                       if (m_shape) {
                         m_constraint->set_shape(m_shape->shape());
                       }
                     },
                     [this]() { return m_shape; }},
                    {"particle_velocity", m_constraint->velocity()}});
  }

  Variant do_call_method(std::string const &name, VariantMap const &) override {
    if (name == "total_force") {
      return shape_based_constraint()->total_force();
    }
    if (name == "min_dist") {
      auto const system = m_system.lock();
      assert(system);
      return shape_based_constraint()->min_dist(
          *system->box_geo, system->cell_structure->local_particles());
    }
    if (name == "total_normal_force") {
      return shape_based_constraint()->total_normal_force();
    }

    return none;
  }

  std::shared_ptr<::Constraints::Constraint> constraint() override {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }
  std::shared_ptr<const ::Constraints::Constraint> constraint() const override {
    return std::static_pointer_cast<::Constraints::Constraint>(m_constraint);
  }
  std::shared_ptr<::Constraints::ShapeBasedConstraint>
  shape_based_constraint() const {
    return m_constraint;
  }

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
    for (auto const &kv : *m_params) {
      do_set_parameter(kv.first, kv.second);
    }
  }

  void
  bind_system(std::shared_ptr<::System::System const> const &system) override {
    assert(m_params);
    assert(system);
    m_system = system;
    shape_based_constraint()->bind_system(system);
    for (auto const &kv : *m_params) {
      do_set_parameter(kv.first, kv.second);
    }
    m_params.reset();
  }

private:
  /* The actual constraint */
  std::shared_ptr<::Constraints::ShapeBasedConstraint> m_constraint;

  /* Keep a reference to the shape */
  std::shared_ptr<Shapes::Shape> m_shape;
};

} // namespace Constraints
} // namespace ScriptInterface
