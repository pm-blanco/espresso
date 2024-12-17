/*
 * Copyright (C) 2021-2024 The ESPResSo project
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

#include "walberla_bridge/Architecture.hpp"

#include "lb_kernels.hpp"

#include "generated_kernels/DynamicUBBDoublePrecisionCUDA.h"
#include "generated_kernels/DynamicUBBSinglePrecisionCUDA.h"
#include "generated_kernels/FieldAccessorsDoublePrecisionCUDA.cuh"
#include "generated_kernels/FieldAccessorsSinglePrecisionCUDA.cuh"
#include "generated_kernels/InitialPDFsSetterDoublePrecisionCUDA.h"
#include "generated_kernels/InitialPDFsSetterSinglePrecisionCUDA.h"
#include "generated_kernels/PackInfoPdfDoublePrecisionCUDA.h"
#include "generated_kernels/PackInfoPdfSinglePrecisionCUDA.h"
#include "generated_kernels/PackInfoVecDoublePrecisionCUDA.h"
#include "generated_kernels/PackInfoVecSinglePrecisionCUDA.h"
#include "generated_kernels/StreamSweepDoublePrecisionCUDA.h"
#include "generated_kernels/StreamSweepSinglePrecisionCUDA.h"

#include "generated_kernels/CollideSweepDoublePrecisionLeesEdwardsCUDA.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalizedCUDA.h"
#include "generated_kernels/CollideSweepSinglePrecisionLeesEdwardsCUDA.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalizedCUDA.h"

namespace walberla {
namespace detail {

using lbmpy::Arch;

template <> struct KernelTrait<double, Arch::GPU> {
  using CollisionModelThermalized =
      pystencils::CollideSweepDoublePrecisionThermalizedCUDA;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepDoublePrecisionLeesEdwardsCUDA;
  using StreamSweep = pystencils::StreamSweepDoublePrecisionCUDA;
  using InitialPDFsSetter = pystencils::InitialPDFsSetterDoublePrecisionCUDA;
  using PackInfoPdf = pystencils::PackInfoPdfDoublePrecisionCUDA;
  using PackInfoVec = pystencils::PackInfoVecDoublePrecisionCUDA;
};

template <> struct KernelTrait<float, Arch::GPU> {
  using CollisionModelThermalized =
      pystencils::CollideSweepSinglePrecisionThermalizedCUDA;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepSinglePrecisionLeesEdwardsCUDA;
  using StreamSweep = pystencils::StreamSweepSinglePrecisionCUDA;
  using InitialPDFsSetter = pystencils::InitialPDFsSetterSinglePrecisionCUDA;
  using PackInfoPdf = pystencils::PackInfoPdfSinglePrecisionCUDA;
  using PackInfoVec = pystencils::PackInfoVecSinglePrecisionCUDA;
};

template <> struct BoundaryHandlingTrait<double, Arch::GPU> {
  using DynamicUBB = lbm::DynamicUBBDoublePrecisionCUDA;
};

template <> struct BoundaryHandlingTrait<float, Arch::GPU> {
  using DynamicUBB = lbm::DynamicUBBSinglePrecisionCUDA;
};

} // namespace detail
} // namespace walberla
