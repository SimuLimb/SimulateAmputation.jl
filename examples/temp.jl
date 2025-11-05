using Severus
using Severus.Comodo
using Severus.Comodo.GLMakie
using Severus.Comodo.GeometryBasics

using Severus.Comodo.Statistics
using Severus.Comodo.LinearAlgebra
using Severus.Geogram
using FileIO

V = [Point{3,Float64}(0.0, 0.0, 0.0), Point{3,Float64}(1.0, 0.0, 0.0), Point{3,Float64}(0.5, 1.0, 0.0)]
F = [TriangleFace{Int}(1,2,3)]

R, P = circumcircle(F,V)