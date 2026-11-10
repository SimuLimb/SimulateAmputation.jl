module SimulateAmputation

using Comodo
using Comodo.GeometryBasics
using Comodo.LinearAlgebra
using Comodo.Rotations
using Comodo.Statistics
using Comodo.BSplineKit
using FileIO
using Geogram
using FEBio 
using Printf

include("functions.jl")

# Export imported modules for later possible use
export Comodo
export FileIO
export Geogram
export FEBio
export Printf

# Export functions
export SimulateAmputationdir, importGeometry, smoothboundary, taperfunction
export bezierclose, set_z_level_boundary!, getboundaryset, cutEnds!,cut_bones
export close_top_bottom, volume_tet_ideal
export mapfunction_D, bezierMap, map_D,new_shape_residuum
export isotropictraction, traction_shrink
export close_bones_bezier
export mapskin, par2shape
end # module Severus
