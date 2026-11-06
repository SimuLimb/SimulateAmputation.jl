using SimulateAmputation
using SimulateAmputation.Comodo
using SimulateAmputation.Comodo.GLMakie
using SimulateAmputation.Comodo.GeometryBasics

using SimulateAmputation.Comodo.Statistics
using SimulateAmputation.Comodo.LinearAlgebra
using SimulateAmputation.Geogram
using FileIO

GLMakie.closeall()
# 
# Z_cut_offset_top = 200.0
# Z_cut_level_skin = 160.0
# Z_cut_level_tibia = Z_cut_level_skin+50
# Z_cut_level_fibula = Z_cut_level_skin+10

# Z_thickness_distal = 20 
Z_cut_offset_top = 200.0
Z_cut_level_skin = 100.0
Z_cut_level_tibia = 150.0 #surgical standards
# Z_cut_level_tibia = Z_cut_level_skin+90
# Z_cut_level_fibula = Z_cut_level_skin+10
Z_cut_level_fibula = Z_cut_level_tibia -10 #surgical standards
Z_thickness_distal = 20
# Z_taper_offset = 80

pointSpacing = 4.0 
pathInputGeom = joinpath(SimulateAmputationdir(),"assets","stl")
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
#top surface cut it, ssnap it to the top, ggremeh it
v_mid_patella = mean(VV[indPatella])

p = Point{3, Float64}(v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]+Z_cut_offset_top) # Point on cutting plane

for modelName in nameSetCut_top
    indexNow = findall(nameSet.==modelName)[1]
    Fn,Vn,_ = trisurfslice(FF[indexNow],VV[indexNow],n,p; output_type=:below, snapTolerance = snapTolerance) 
    np = spacing2numvertices(Fn,Vn,pointSpacing)
    Fn,Vn = ggremesh(Fn,Vn; nb_pts = np)
    set_z_level_boundary!(Fn,Vn,p[3])
    FF[indexNow] = Fn
    VV[indexNow] = Vn
end
# cut bottom surface, smoothen it 
nameSetCut_down = ("skin","tibia","fibula")
Z_cut_offset_down =[Z_cut_level_skin,Z_cut_level_tibia,Z_cut_level_fibula] 
for model in eachindex(nameSetCut_down)
    Z_cut_down = Z_cut_offset_down[model]
    modelName = nameSetCut_down[model]
    indexNow = findall(nameSet .== modelName)[1]    
    p = Point{3,Float64}(v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]-Z_cut_down) # Point on cutting plane   
    Fn,Vn,_ = trisurfslice(FF[indexNow],VV[indexNow],n,p; output_type=:above, snapTolerance = snapTolerance) 
    FF[indexNow] = Fn
    VV[indexNow] = Vn
end 

# obtain the boundary sets of tibia, fibula 

V_tibia = VV[indTibia]
Eb_tibia = boundaryedges(FF[indTibia])
indCurve_tibia = edges2curve(Eb_tibia)
indCurve_tibia = indCurve_tibia[1:end-1]
V_fibula = VV[indFibula]
Eb_fibula = boundaryedges(FF[indFibula])
indCurve_fibula = edges2curve(Eb_fibula)
indCurve_fibula = indCurve_fibula[1:end-1]



v_mid_low_tibia = mean(V_tibia[indCurve_tibia])

p_distal = Point{3,Float64}(v_mid_low_tibia[1],v_mid_low_tibia[2],v_mid_low_tibia[3]-Z_thickness_distal)

zlevel = v_mid_patella[3]
Fs = FF[1]
Vs = VV[1]
#skin boundary edges isolate the top and bottom can be cleaned further 
Eb = boundaryedges(Fs)
Eb_low = Vector{eltype(Eb)}()
for e in Eb         
    if Vs[e[1]][3] < zlevel        
      push!(Eb_low,e)  
    end
end
indCurve = edges2curve(Eb_low)
indCurve = indCurve[1:end-1] 

#nsmoothen the skin, tibia and fibula 

nameSetSmooth = ("skin", "tibia", "fibula")
Z_level_smooth = [Z_thickness_distal, nothing, nothing]
for model in eachindex(nameSetSmooth)
    Z_thickness_distal = Z_level_smooth[model]
    modelName = nameSetSmooth[model]
    indexNow = findall(nameSet .== modelName)[1]    
    F_boundary = FF[indexNow]
    V_boundary = VV[indexNow]
    Fs = FF[1]
    Vs = VV[1]
    Eb = boundaryedges(F_boundary)
    Vs_new = smoothboundary(VV[indexNow], FF[indexNow],Eb; Z_thickness_distal=Z_thickness_distal)
    VV[indexNow] = Vs_new
end 
#bezeir curve
VV[1] = Vs
bezier = ("skin", "tibia", "fibula")
velocity = ((100,100),(25,20),(25,20))
z_distance = [p_distal, nothing, nothing]
boundary_curve = [indCurve,indCurve_tibia,indCurve_fibula]
numPointBezier = 250
for model in eachindex(bezier)
    v_vec = velocity[model]
    z_dist = z_distance[model]
    ind_curve= boundary_curve[model]
    modelName = bezier[model]
    indexNow = findall(nameSet .== modelName)[1]    
    # println(indexNow)
    Fn,Vn = bezierclose(FF[indexNow],VV[indexNow],v_vec[1],v_vec[2],ind_curve; numPointBezier = numPointBezier,p_distal=z_dist, numEdgeSmoothSteps = 25, λ = 0.5, pointSpacing=pointSpacing)
    FF[indexNow] = Fn
    VV[indexNow] = Vn    
end 
F_skin = FF[indSkin]
V_skin = VV[indSkin]

#fixing the zlvel 
set_z_level_boundary!(F_skin,V_skin,v_mid_patella[3]+Z_cut_offset_top)

#closing the top surface 
#trimesh the skin and then the femur 
V_femur = VV[indFemur]
F_femur = FF[indFemur]
boundary_set_skin, z_level_skin = getboundaryset(F_skin,V_skin)
boundary_set_femur, z_level_femur = getboundaryset(F_femur,V_femur)

VT = (V_skin[boundary_set_skin[1]], V_femur[boundary_set_femur[1]]) # Curves

R = ([1,2],) 
# R_femur = ([2],) 
P = (pointSpacing,pointSpacing)  # Point spacings
Ft,Vt,Ct = regiontrimesh(VT,R,P)
# Ffe,Vfe,Cfe = regiontrimesh(VT,R_femur,P)
# F_skin,V_skin,Cskin = joingeom(F_skin,V_skin,Ft,Vt)
# F_skin,V_skin = mergevertices(F_skin,V_skin)

# F_femur,V_femur,Cfemur = joingeom(F_femur,V_femur,Ffe,Vfe)
VV[indSkin] = V_skin
FF[indSkin] = F_skin
# VV[indFemur] = V_femur
# FF[indFemur] = F_femur
#ggremesh? 

## Visualization

fig = Figure(size=(1200,800))

ax1 = LScene(fig[1,1])
for i in 1:1:length(fileName_set)
    hp1 = poly!(ax1,GeometryBasics.Mesh(VV[i],FF[i]), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
end
hp2 = poly!(ax1,GeometryBasics.Mesh(Vt,Ft), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
display(GLMakie.Screen(),fig)

#Tetgen 
#tetgen only for bones 

stringOpt = "paAqY"

function volume_region_length(F,V)
    vol_re = mean(V)
    D = edgelengths(F,V)
    vol = mean(D)^3 / (6.0*sqrt(2.0))
    return vol,vol_re
end 

bone_set = [indFemur,indPatella,indTibia,indFibula]
region_vol = []
V_regions = []

for index in bone_set 
    Fn = FF[index]
    Vn = VV[index]
    vol,vol_re = volume_region_length(Fn,Vn)
    push!(region_vol,vol)
    push!(V_regions,vol_re)
end

region_skin = mean(VV[indSkin]) .-30
D = edgelengths(FF[indSkin],VV[indSkin])
vol_skin = mean(D)^3 / (6.0*sqrt(2.0))
push!(V_regions,region_skin)
push!(region_vol,vol_skin)
region_top = mean(Vt) .-30 
D = edgelengths(Ft,Vt)
vol_top = mean(D)^3 / (6.0*sqrt(2.0))
push!(V_regions,region_top)
push!(region_vol,vol_top)


Fw,Vw,Cw= joingeom(FF[indFemur],VV[indFemur],FF[indPatella],VV[indPatella],FF[indTibia],VV[indTibia],FF[indFibula],VV[indFibula],FF[indSkin],VV[indSkin],Ft,Vt)
Fw,Vw = mergevertices(Fw,Vw)
E,V,CE,Fb,Cb= tetgenmesh(Fw,Vw; facetmarkerlist=Cw, V_regions=V_regions,region_vol=region_vol, stringOpt)


#Visualization
cmap = cgrad(:Spectral, 5, categorical = true)

F = element2faces(E) # Triangular faces
CE_F = repeat(CE,inner=4) #why is it repeating

Fbs,Vbs = separate_vertices(Fb,V)
Cbs_V = simplex2vertexdata(Fbs,Cb)

Fs,Vs = separate_vertices(F,V)
CE_Vs = simplex2vertexdata(Fs,CE_F)
M = GeometryBasics.Mesh(Vs,Fs)

strokewidth = 1 

indNodesTop = unique(reduce(vcat,Fb[Cb.==6]))
indNodesBones = unique(reduce(vcat,Fb[[in(c,[1,2,3,4]) for c in Cb]]))
fig3 = Figure(size=(800,800))

ax1 = Axis3(fig3[1, 1][1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary surfaces")
hp1 = mesh!(ax1,GeometryBasics.Mesh(Vbs,Fbs), color=Cbs_V, shading = FastShading, transparency=true, overdraw=false,colorrange = (1,6),colormap=cmap)

display(GLMakie.Screen(),fig3)

