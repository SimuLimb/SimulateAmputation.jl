using Severus
using Severus.Comodo
using Severus.Comodo.GeometryBasics
using Severus.Comodo.GLMakie
using Severus.Comodo.LinearAlgebra
using Severus.Comodo.Statistics
using Severus.FEBio
using Severus.FEBio.XML
using Severus.Printf 


#Full leg 
Z_cut_offset_top = 200.0
Z_cut_offset_bottom = 250.0

pointSpacing = 5.0 
pathInputGeom = joinpath(severusdir(),"assets","stl")

fileName_set = ("FMA7163_right_leg_isolated_remesh_25k.stl", # Skin
"FMA24474.stl", # Right femur
"FMA24486.stl", # Right patella
"FMA24477.stl", # Right tibia
"FMA24480.stl") # Right fibula
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

# importing the geometry form the folder 
FF, VV = importGeometry(pathInputGeom, fileName_set; pointSpacing = pointSpacing)
FF_ori = deepcopy(FF)
VV_ori = deepcopy(VV)

# Cutting surfaces 
# the reference point is at patella 
#top surface cut it, ssnap it to the top, ggremeh it
v_mid_patella = mean(VV[indPatella])
zCut_level_top = v_mid_patella[3]+Z_cut_offset_top
zCut_level_bottom = v_mid_patella[3]-Z_cut_offset_bottom

cutEnds!(FF, VV, [zCut_level_top, zCut_level_bottom],n;snapTolerance,pointSpacing)
#generating the top and bottom boundary 
bone_set_bound = ("femur","tibia","fibula") # bone set used for closing the top and bottom 
# # obtain the boundary sets of tibia, fibula 
    # # ## Visualization
Ft,Vt,Ft_b,Vt_b = close_top_bottom(FF,VV,nameSet,bone_set_bound;pointSpacing=pointSpacing)
# # #Tetgen 
# # #tetgen only for bones 

stringOpt = "paAqY"

region_vol = [] #the interior volume for the residum 
V_regions = [] #interior points inside for the geometry to be meshed 
#instrior point inside the skin manually adjusted 
indFace_skin = Int(round(length(FF[indSkin])*0.47))
point_interior = faceinteriorpoint(FF[indSkin],VV[indSkin],indFace_skin ; w = 0.5)
point_interior = Point{3,Float64}(point_interior[1],point_interior[2],point_interior[3]+22.0) 
#ideal tet hed volume 
vol_skin = 15.0*volume_tet_ideal(FF[indSkin],VV[indSkin])
#interior points that need to be meshed 
push!(V_regions,point_interior)
push!(region_vol,vol_skin)
#bone set is not meshed, they are empty regions in the mesh
bone_set = ("femur","patella","tibia","fibula")
for modelName in bone_set 
    indexNow = findall(nameSet.==modelName)[1]
    F = FF[indexNow]
    V = VV[indexNow]
    indFace = Int(round(length(F)/2))
    interior_point =  faceinteriorpoint(F,V,indFace;w = 0.5 )   
    vol = volume_tet_ideal(F,V)
    push!(region_vol,vol)
    push!(V_regions,interior_point)
end  
V_holes = V_regions[2:5]
V_regions=V_regions[1]
region_vol=region_vol[1]

#the whole geometry is meshed
Fw,Vw,Cw= joingeom(FF[indSkin],VV[indSkin],FF[indFemur],VV[indFemur],FF[indPatella],VV[indPatella],FF[indTibia],VV[indTibia],FF[indFibula],VV[indFibula],Ft,Vt,Ft_b,Vt_b)
Fw,Vw = mergevertices(Fw,Vw)
fig1 = Figure(size=(800,800))
color = [:white, 3 ,5, 7, 9]
ax1 = LScene(fig1[1,1],show_axis =false)
for i in 1:1:length(fileName_set)
    hp1 = poly!(ax1,GeometryBasics.Mesh(VV[i],FF[i]), color= color[i], colormap = :viridis, colorrange = (1, 9), shading = FastShading, transparency=true)
end

# hp3 = poly!(ax1,GeometryBasics.Mesh(V_shrunk,F_shrunk), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
display(GLMakie.Screen(),fig1)

# # # #tetgen
   
E,V,CE,Fb,Cb= tetgenmesh(Fw,Vw; facetmarkerlist=Cw, V_regions=[V_regions], region_vol=[region_vol], V_holes = V_holes,  stringOpt)
#atrophy part 
κ = 1.0
ratio = 100.0 #ratio = κ/μ

J_desired = 1/(2^2)
exe_path ="C:/Program Files/FEBioStudio/bin/febio4"
#######
UT,VT,FF,VV,F_shrunk,V_shrunk,U_mag,U_vec = shrinky(exe_path,FF,VV,E,V,CE,Fb,Cb,κ,ratio,J_desired)
# Visualization
fig2 = Figure(size=(800,800))

ax1 = LScene(fig2[1,1])
# hp3 = poly!(ax1,GeometryBasics.Mesh(VT[end],Fb[Cb.==1]), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
hp3 = poly!(ax1,GeometryBasics.Mesh(V_shrunk,F_shrunk), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
display(GLMakie.Screen(),fig2)
# # normalplot(ax1,Ft,Vt)


 ######################mapping the coordinates##############################
V_skin = VV[indSkin]
F_skin = FF[indSkin]
#actual cut of the limb can be user generated 
Z_taper_start = v_mid_patella[3]
Z_taper_end = v_mid_patella[3]-127.0 
#normalizing in the z direction returning the D and the correspoinding indices 
D,ind = mapfunction_D(V_skin,Z_taper_start,Z_taper_end)
D_map = D[ind]
#linear
# T_bez_in = [0.0, 0.5, 1.0]
# T_bez_map_in = [0.0, 0.5, 1.0]
T_bez_in1     = [0.0, 0.2, 0.5, 0.7]
T_bez_map_in1 = [0.0, 0.0, 1.0, 1.0]
T_bez1, T_map_bez1 = bezierMap(T_bez_in1,T_bez_map_in1; n=100)

T_bez_in2     = [0.7, 0.9, 1.0]
T_bez_map_in2 = [1.0, 1.0, 0.5]
T_bez2, T_map_bez2 = bezierMap(T_bez_in2,T_bez_map_in2; n=100)

T_bez = T_bez1
append!(T_bez,T_bez2[2:end])
T_map_bez = T_map_bez1
append!(T_map_bez,T_map_bez2[2:end])

#control points in between
# T_bez_in = [0.0, 0.25, 0.25, 1.0]
# T_bez_map_in = [0.0, 0.0, 0.25, 1.0]
#less than linear
#  T_bez_in = [0.0, 0.1, 0.1, 1.0]
#  T_bez_map_in = [0.0, 0.0, 0.75, 1.0]
#correspoinding spline 
# T_bez, T_map_bez = bezierMap(T_bez_in,T_bez_map_in; n=200)
#mapping it to the skin 
D_map = map_D(T_bez, T_map_bez, D_map; spline_order=4)
Fs,Vs = new_shape_residum(V_skin,F_skin,D_map,ind,U_vec,v_mid_patella)

# # normalplot(ax1,Ft,Vt)
 display(GLMakie.Screen(),fig)
#cut bones tibia and fibula 
#cut bones
Z_cut_level_tibia = 150.0 #surgical standards
Z_cut_level_fibula = Z_cut_level_tibia -10 #surgical standards
nameSetCut_down = ("tibia","fibula")
Z_cut_offset_down =[Z_cut_level_tibia,Z_cut_level_fibula] 
output_type=:above
FF,VV = cut_bones(FF,VV,nameSet,nameSetCut_down,Z_cut_offset_down,v_mid_patella,n;output_type,snapTolerance)
#bezeir close for bones 
bezier_bones = ("tibia", "fibula")
velocity = ((25,20),(25,20))
z_distance = [nothing, nothing]
# boundary_curve = [indCurve_tibia,indCurve_fibula]
numPointBezier = 250
indCurve_tibia,_ = getboundaryset(FF[indTibia],VV[indTibia])
V_tibia = VV[indTibia]
F_tibia = VV[indTibia]

v_mid_low_tibia = mean(V_tibia[indCurve_tibia[1]])
v_mid_low_tibia = mean(V_tibia[indCurve_tibia[1]])

FF,VV = close_bones_bezier(FF,VV,nameSet,bezier_bones,velocity,z_distance;numPointBezier,pointSpacing)

#close residum skin
boundary_set_s1,_ = getboundaryset(Fs,Vs)
# indCurve_tibia,_ = getboundaryset(FF[indTibia],VV[indTibia])
# V_tibia=VV[indTibia]

numPointBezier = 250
Z_thickness_distal = 30
p_distal = Point{3,Float64}(v_mid_low_tibia[1],v_mid_low_tibia[2],v_mid_low_tibia[3]-Z_thickness_distal)
Fs2,Vs2 = bezierclose(Fs,Vs,150,75,boundary_set_s1[2]; numPointBezier = numPointBezier,p_distal=p_distal, numEdgeSmoothSteps = 25, λ = 0.5, pointSpacing=pointSpacing)
V_bez = Vs2
F_bez = Fs2

#visualization
fig = Figure(size=(800,800))
# ax1 = LScene(fig[1,1])
ax1 = Axis3(fig[1, 1], clip =false, viewmode = :free, aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Test")

for i in 2:1:length(fileName_set)
    hp1 = poly!(ax1,GeometryBasics.Mesh(VV[i],FF[i]), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)
end
#  hp4 = poly!(ax1,GeometryBasics.Mesh(Vs,Fs), color=:white, shading = FastShading, transparency=true)
hp3 = poly!(ax1,GeometryBasics.Mesh(V_bez,F_bez), color=:white, shading = FastShading, transparency=false,strokecolor=:black,strokewidth=1)

# hp3 = poly!(ax1,GeometryBasics.Mesh(V_shrunk,F_shrunk), color=:blue, shading = FastShading, transparency=true)
 display(GLMakie.Screen(),fig)