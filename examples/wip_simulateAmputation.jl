using Severus
using Severus.Comodo
using Severus.Comodo.GLMakie
using Severus.Comodo.GeometryBasics

using Severus.Comodo.Statistics
using Severus.Comodo.LinearAlgebra
using Severus.Geogram
using FileIO

# 
Z_cut_offset_top = 200.0
Z_cut_level_skin = 160.0
Z_cut_level_tibia = Z_cut_level_skin+50
Z_cut_level_fibula = Z_cut_level_skin+10

Z_thickness_distal = 20 

pointSpacing = 4.0 
pathInputGeom = joinpath(severusdir(),"assets","stl")
fileName_set = ("FMA7163_right_leg_isolated_remesh_25k.stl",
"FMA24474.stl",
"FMA24486.stl",
"FMA24477.stl",
"FMA24480.stl")

nameSet = ("skin","femur","patella","tibia","fibula")
indSkin = 1
indFemur = 2
indPatella = 3
indTibia = 4
indFibula = 5

nameSetCut_top = ("skin","femur")
# nameSetCut_bottom = ("skin","femur","tibia","fibula")

n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))# Cutting plane normal
snapTolerance = 1e-6

# FMAID 	English name
# 7163  	Skin
# 16586 	Right hip bone
# 24474 	Right femur
# 24486     Right patella
# 24477 	Right tibia
# 24480 	Right fibula

FF = []
VV = []
for i in 1:1:length(fileName_set)
    M = load(joinpath(pathInputGeom,fileName_set[i]))
    F = tofaces(faces(M))
    V = topoints(coordinates(M))
    V = [Point{3,Float64}(v) for v in V]
    
    F,V,_,_ = mergevertices(F,V; pointSpacing=pointSpacing)

    # Remesh evenly using desired point spacing

    # Determine number of points to request to closely match desired point spacing
    np = spacing2numvertices(F,V,pointSpacing) 

    # Remesh with desired number of points 
    F,V = ggremesh(F,V; nb_pts=np)

    push!(FF,F)
    push!(VV,V)
end

FF_ori = deepcopy(FF)
VV_ori = deepcopy(VV)

# Cutting surfaces 
v_mid_patella = mean(VV[indPatella])

p = Point{3, Float64}(v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]+Z_cut_offset_top) # Point on cutting plane

for modelName in nameSetCut_top
    indexNow = findall(nameSet.==modelName)[1]
    Fn,Vn,_ = trisurfslice(FF[indexNow],VV[indexNow],n,p; output_type=:below, snapTolerance = snapTolerance) 
    FF[indexNow] = Fn
    VV[indexNow] = Vn
end

p = Point{3,Float64}(v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]-Z_cut_level_skin) # Point on cutting plane
Fn,Vn,_ = trisurfslice(FF[indSkin],VV[indSkin],n,p; output_type=:above, snapTolerance = snapTolerance) 
FF[indSkin] = Fn
VV[indSkin] = Vn

p = Point{3,Float64}(v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]-Z_cut_level_tibia) # Point on cutting plane
Fn,Vn,_ = trisurfslice(FF[indTibia],VV[indTibia],n,p; output_type=:above, snapTolerance = snapTolerance) 
FF[indTibia] = Fn
VV[indTibia] = Vn

p = Point{3,Float64}(v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]-Z_cut_level_fibula) # Point on cutting plane
Fn,Vn,_ = trisurfslice(FF[indFibula],VV[indFibula],n,p; output_type=:above, snapTolerance = snapTolerance) 
FF[indFibula] = Fn
VV[indFibula] = Vn


V_tibia = VV[indTibia]
Eb_tibia = boundaryedges(FF[indTibia])
indCurve_tibia = edges2curve(Eb_tibia)
indCurve_tibia = indCurve_tibia[1:end-1]
v_mid_low_tibia = mean(V_tibia[indCurve_tibia])

p_distal = Point{3,Float64}(v_mid_low_tibia[1],v_mid_low_tibia[2],v_mid_low_tibia[3]-Z_thickness_distal)

zlevel = v_mid_patella[3]

Fs = FF[1]
Vs = VV[1]

Eb = boundaryedges(Fs)
Eb_low = Vector{eltype(Eb)}()
for e in Eb         
    if Vs[e[1]][3] < zlevel        
      push!(Eb_low,e)  
    end
end
indCurve = edges2curve(Eb_low)
indCurve = indCurve[1:end-1] 

# Vs[indCurve] = evenly_sample(Vs[indCurve],length(indCurve); rtol = 1e-8, niter = 1)

# Vs = smoothmesh_hc(Eb,Vs, 25, 0.1, 0.5)

Vs = smoothmesh_laplacian(Eb,Vs, 25, 0.5)


ZF = [mean(Vs[f])[3] for f in Fs]
logicSmooth = ZF.<(minimum(ZF)+20)
constrained_points = elements2indices(boundaryedges(Fs[logicSmooth]))
Vs = smoothmesh_hc(Fs[logicSmooth],Vs, 50, 0.1, 0.5; constrained_points=constrained_points)

# FF[1] = Fs[logicSmooth]
VV[1] = Vs

NV = vertexnormal(Fs,Vs; weighting=:size) # May be inefficient to do this to all vertices (/faces)
NV_indCurve = NV[indCurve] # Vertex normals for the curve
NE_indCurve = Vector{Point{3,Float64}}(undef,length(indCurve))
NB_base = Vector{Point{3,Float64}}(undef,length(indCurve))
m = length(indCurve)
for i in eachindex(indCurve)
    ii = indCurve[i]
    jj = indCurve[mod1(i+1,m)]
    NE_indCurve[i] = normalizevector(Vs[jj] - Vs[ii])

    NB_base[i] = cross(NE_indCurve[i],NV_indCurve[i])
end

numEdgeSmoothSteps = 25

λ = 0.5
for _ in 1:1:numEdgeSmoothSteps
    for i in eachindex(indCurve)
        i_prev = mod1(i-1,m)
        i_next = mod1(i+1,m)
        NB_base[i] = ((1.0-λ) .* NB_base[i]) .+  λ.*(NB_base[i_prev] .+ NB_base[i_next])./2.0
        # NB_base[i] = 0.5.*NB_base[i] .+ 0.25.*NB_base[i_prev] .+ 0.25.*NB_base[i_next]
        NB_base[i] = normalizevector(NB_base[i])
    end
end






numPointBezier = 250
v1 = 150 # Departure velocity of Hermite equivalent spline
v2 = 75 # Arrival velocity of Hermite equivalent spline
VV_B = []
Vb = Vector{Point{3,Float64}}()
for i in eachindex(indCurve)
    nBase = NB_base[i]
    # nEnd = NV_indCurve[i]
    nEnd = normalizevector(Point{3,Float64}(Vs[indCurve[i]][1]-p_distal[1],Vs[indCurve[i]][2]-p_distal[2],0.0))

    pb = [Vs[indCurve[i]], Vs[indCurve[i]] + v1/3*nBase, p_distal + v2/3*nEnd,p_distal]
    vb = nbezier(pb,numPointBezier)
    push!(VV_B,vb)
    append!(Vb,vb)
end

Fb = grid2surf(Vb,length(indCurve);periodicity=(false,true),face_type=:quad)

Fbq = quad2tri(Fb,Vb)

F_skin,V_skin,C_skin = joingeom(FF[1],VV[1],Fbq,Vb)

F_skin,V_skin,_ = mergevertices(F_skin,V_skin)

# Determine number of points to request to closely match desired point spacing
np = spacing2numvertices(F_skin,V_skin,pointSpacing) 

# Remesh with desired number of points 
F_skin,V_skin = ggremesh(F_skin,V_skin; nb_pts=np)

## Visualization

fig = Figure(size=(1200,800))

# ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Imported mesh")
ax1 = LScene(fig[1,1])
for i in 2:1:length(fileName_set)
    hp1 = poly!(ax1,GeometryBasics.Mesh(VV[i],FF[i]), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
end

# hp2 = wireframe!(ax1,GeometryBasics.Mesh(Vs,Eb_low), linewidth=5,color=:red)


# hp2 = lines!(ax1,Vs[indCurve],color=:black,linewidth=3)
# scatter!(ax1,Vs[indCurve], color = :black, markersize = 25)

# # hp1 = poly!(ax1,GeometryBasics.Mesh(Vn,Fn), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=0.25)

# hp3 = dirplot(ax1,Vs[indCurve],NV_indCurve.*mean(edgelengths(Fs,Vs)); color=:red,linewidth=3,scaleval=1.0,style=:from)
# hp4 = dirplot(ax1,Vs[indCurve],NE_indCurve.*mean(edgelengths(Fs,Vs)); color=:green,linewidth=3,scaleval=1.0,style=:from)
# hp5 = dirplot(ax1,Vs[indCurve],NB_base.*mean(edgelengths(Fs,Vs))*3; color=:cyan,linewidth=3,scaleval=1.0,style=:from)

# hp6 = lines!(ax1,V_tibia[indCurve_tibia],color=:yellow,linewidth=6)

# for i in eachindex(indCurve)
#     lines!(ax1,VV_B[i],color=:black,linewidth=1)
# end

# hp1 = poly!(ax1,GeometryBasics.Mesh(Vb,Fb), color=:blue, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=0.25)
# hp7 = mesh!(ax1,GeometryBasics.Mesh(Vb,Fb), color=:white, shading = FastShading, transparency=true)
# hp7 = mesh!(ax1,GeometryBasics.Mesh(V_skin,F_skin), color=:white, shading = FastShading, transparency=true)
hp7 = poly!(ax1,GeometryBasics.Mesh(V_skin,F_skin), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=0.25)

fig

# save("/home/kevin/Desktop/leg.stl",GeometryBasics.Mesh(V_skin,F_skin))