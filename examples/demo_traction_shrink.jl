using SimulateAmputation
using SimulateAmputation.Comodo
using SimulateAmputation.Comodo.GeometryBasics
using SimulateAmputation.Comodo.GLMakie
using SimulateAmputation.Comodo.LinearAlgebra
using SimulateAmputation.Comodo.Statistics
using SimulateAmputation.Comodo.Rotations
using SimulateAmputation.FEBio
using SimulateAmputation.FEBio.XML
using SimulateAmputation.Printf 
using SimulateAmputation.Geogram
using FileIO

GLMakie.closeall()

# Control parameters 
testCase = 2
if testCase == 1
    pathInputGeom = joinpath(SimulateAmputationdir(),"assets","stl")
    fileName_set = ("FMA7163_right_leg_isolated_remesh_25k.stl", # Skin
    "FMA24474.stl", # Right femur
    "FMA24486.stl", # Right patella
    "FMA24477.stl", # Right tibia
    "FMA24480.stl") # Right fibula
elseif testCase == 2
    pathInputGeom = joinpath(SimulateAmputationdir(),"assets","stl","VHF")
    fileName_set = ("VHF_right_skin_leg_5p0_mm.obj", # Skin
    "VHF_right_femur.obj", # Right femur
    "VHF_right_patella.obj", # Right patella
    "VHF_right_tibia.obj", # Right tibia
    "VHF_right_fibula.obj") # Right fibula
end

nameSet = ["skin","femur","patella","tibia","fibula"]
indSkin = 1
indFemur = 2
indPatella = 3
indTibia = 4
indFibula = 5

pointSpacing = 6.0 
pointSpacingExport = 6.0

# Cutting parameters for full leg (used for shrinkage simulation)
Z_cut_offset_top = 200.0
Z_cut_offset_bottom = 250.0

Z_cut_level_tibia = 150.0 #surgical standards
Z_cut_level_fibula = Z_cut_level_tibia -10 #surgical standards
Z_cut_offset_skin =  maximum([Z_cut_level_tibia, Z_cut_level_fibula]) #127

Z_taper_offset_top = 30.0
Z_taper_offset_bottom = Z_cut_offset_skin+15

Z_thickness_distal = 30

numPointBezier = 250
v1_bezier_distal = 75.0
v2_bezier_distal = 75.0

# Material parameters and bc's for shrinkage
κ = 1.0
κ_ratio = 100.0 # κ_ratio = κ/μ
J_desired = 0.2 # 1/(2^2)
exe_path ="febio4"

nCut = normalizevector(Vec{3, Float64}(0.0, 0.0, 1.0)) # Cutting plane normal

# importing the geometry form the folder 
FF, VV = importGeometry(pathInputGeom, fileName_set; pointSpacing = pointSpacing)
# VV[5] .-= Point{3,Float64}(1.0, 0.0, 1.0) # shift fibula 

# Centre on patella
v_mid_patella = surface_centroid(FF[indPatella], VV[indPatella])
v_mid_femur = surface_centroid(FF[indFemur], VV[indFemur])
for (i,V) in enumerate(VV)
    VV[i] = V.-v_mid_patella
end

if testCase == 2
    # Reorient based on tibia direction 
    svd_reorient = pointsvd(VV[indTibia])
    R = deepcopy(svd_reorient.Vt)
    nz = R[:,1]
    Q = rotation_between([0.0, 0.0, 1.0], nz)
    for (i,V) in enumerate(VV)
        VV[i] = [Point{3, Float64}(Q*v) for v ∈ V] # Rotate the coordinates
    end

    Q = RotXYZ(0.0, 0.0, (3/180)*pi)
    for (i,V) in enumerate(VV)
        VV[i] = [Point{3, Float64}(Q*v) for v ∈ V] # Rotate the coordinates
    end
end

# nx = v_mid_patella-v_mid_femur
# nx = Point{3,Float64}(nx[1], nx[2], 0.0)
# nx = normalizevector(nx)
# Q = rotation_between([1.0, 0.0, 0.0], nx)
# for (i,V) in enumerate(VV)
#     VV[i] = [Point{3, Float64}(Q*v) for v ∈ V] # Rotate the coordinates
# end

# Ensure Patella points forward

FF_ori = deepcopy(FF)
VV_ori = deepcopy(VV)


# Visualisation of surface 
cMap_alpha = cgrad(:viridis, length(FF), categorical = true, alpha=0.2)
cMap = cgrad(:viridis, length(FF), categorical = true, alpha=1.0)
fig = Figure(size=(800,800))
ax1 = AxisGeom(fig[1, 1], title = "Surface meshes")
# hp1 = meshplot!(ax1, Fw, Vw; color=RGBAf(0.95,0.95,0.95,0.5), transparency=true, strokewidth=0.0)

for i in eachindex(FF)
    if i==1
        hp1 = meshplot!(ax1, FF[i], VV[i]; color=cMap_alpha[i], strokewidth=0.0, transparency = true)            
    else
        hp1 = meshplot!(ax1, FF[i], VV[i]; color=cMap[i], strokewidth=0.0, transparency = false)            
    end    
end

display(GLMakie.Screen(),fig)

# -------------- Cutting surfaces 
# The reference point for cutting is the centre of the patella. The top surface 
# is cut, snapped to the cut level, and ggremesh is used for remeshing. 
v_mid_patella = surface_centroid(FF[indPatella], VV[indPatella])
zCut_level_top = v_mid_patella[3]+Z_cut_offset_top
zCut_level_bottom = v_mid_patella[3]-Z_cut_offset_bottom
Z_taper_start = v_mid_patella[3] + Z_taper_offset_top
Z_taper_end = v_mid_patella[3] - Z_taper_offset_bottom
Z_cut_level_skin = v_mid_patella[3] - Z_cut_offset_skin

cutEnds!(FF, VV, [zCut_level_top, zCut_level_bottom], nCut; snapTolerance=1e-6, pointSpacing=pointSpacing)

# Closing over the top and bottom
bone_set_bound = ("femur","tibia","fibula") # bone set used for closing the top and bottom 
Ft,Vt,Ft_b,Vt_b = close_top_bottom(FF,VV,nameSet,bone_set_bound; pointSpacing=pointSpacing)

# -------------- Tetrahedral meshing using tetgen 
# Create combined and closed surface set
Fw,Vw,Cw= joingeom(FF[indSkin],VV[indSkin],FF[indFemur],VV[indFemur],FF[indPatella],VV[indPatella],FF[indTibia],VV[indTibia],FF[indFibula],VV[indFibula],Ft,Vt,Ft_b,Vt_b)
Fw,Vw = mergevertices(Fw,Vw) # Merge nodes so they are new_shape_residum

# Visualisation of surface 
cMap_alpha = cgrad(:viridis, maximum(Cw), categorical = true, alpha=0.2)
cMap = cgrad(:viridis, maximum(Cw), categorical = true, alpha=1.0)
fig = Figure(size=(800,800))
ax1 = AxisGeom(fig[1, 1], title = "Surface meshes")
# hp1 = meshplot!(ax1, Fw, Vw; color=RGBAf(0.95,0.95,0.95,0.5), transparency=true, strokewidth=0.0)

for c in 1:1:maximum(Cw)
    if c==1
        hp1 = meshplot!(ax1, Fw[Cw.==c], Vw; color=cMap_alpha[c], strokewidth=0.0, transparency = true)            
    else
        hp1 = meshplot!(ax1, Fw[Cw.==c], Vw; color=cMap[c], strokewidth=0.0, transparency = false)            
    end    
end

display(GLMakie.Screen(),fig)

# Define hole and region point 
V_hole = faceinteriorpoint(FF[indPatella],VV[indPatella],1;w = 0.5 ) # Patella forms a hole
V_region = faceinteriorpoint(Fw,Vw,Int(round(length(FF[indSkin])/2)); w = 0.5 ) # This is a point inside the tissue to be meshed
region_vol = 5.0*volume_tet_ideal(Fw,Vw) # Tetrahedral volume approximate in tissue away from boundary 

# Visualisation of surface 
fig = Figure(size=(800,800))
ax1 = AxisGeom(fig[1, 1], title = "Limb surfaces meshes ready for Tetgen")
hp1 = meshplot!(ax1, Fw, Vw; color=RGBAf(0.95,0.95,0.95,0.5), transparency=true, strokewidth=0.0)
hp2 = scatter!(ax1, V_region, color=:red, markersize=15)
Legend(fig[1,2], [hp1, hp2], ["Surface", "Interior point"])
display(GLMakie.Screen(),fig)

stringOpt = "paAqYQ"
E, V, CE, Fb, Cb= tetgenmesh(Fw,Vw; facetmarkerlist=Cw, V_regions=[V_region], region_vol=[region_vol], V_holes = [V_hole],  stringOpt)

# -------------- Use FEBio to simulate atrophy by reducing element volumes 
UT,VT,FF,VV,F_shrunk,V_shrunk,U_mag,U_vec, JT_E, DD_J = traction_shrink(exe_path, FF, VV, E, V, CE, Fb, Cb, κ, κ_ratio, J_desired)

# Visualisation of FEBio shrinking results 
numInc = length(UT)
incRange = 0:1:numInc-1
stepStart = 0#incRange[end]

J_mean = mean(JT_E[stepStart+1])
VE = simplexcenter(E,V)
B = [v[1]>-90.0 for v in VE]

F = element2faces(E[B])
# indB = boundaryfaceindices(F)
J_F = repeat(JT_E[stepStart+1][B], inner=4)
# Fb = F[indB]

Fs,Vs = separate_vertices(F,VT[stepStart+1])
J_Vs = simplex2vertexdata(Fs,J_F,Vs)

vol0 = sum(tetvolume(E,VT[1]))

## Visualise shrinkage 

fig = Figure(size=(800,800))
ax = AxisGeom(fig[1, 1], title = "Step: $stepStart. J = " * @sprintf("%.3f", J_mean))
hp = meshplot!(ax, Fs, Vs; strokewidth=0.0, color=J_Vs , colormap = :viridis, colorrange=(J_desired-1e-3,1.0))

Colorbar(fig[1, 2], hp.plots[1], label = "Volume ratio []", ticks = J_desired:0.1:1.0, labelsize=35, ticklabelsize=25) 

hSlider = Slider(fig[2, 1], range = incRange, startvalue = stepStart,linewidth=30)
on(hSlider.value) do stepIndex         
    J_F = repeat(JT_E[stepIndex+1][B], inner=4)
    Fs, Vs = separate_vertices(F, VT[stepIndex+1])
    J_Vs = simplex2vertexdata(Fs,J_F)

    vol = sum(tetvolume(E,VT[stepIndex+1]))

    hp[1] = GeometryBasics.Mesh(Vs, Fs)
    hp.color = J_Vs    

    J_tot = vol/vol0
    # hp.color = UT_mag[stepIndex+1]
    ax.title = "Step: $stepIndex. J = " * @sprintf("%.3f", J_tot)
end

screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "FEBio example")

## Visualise J distribution graph

xlim=(0,1.01)
ylim=(0,50)
fig = Figure(size=(800,800), fontsize = 25)
ax = Axis(fig[1,1], xlabel="J", xticks = xlim[1]:0.1:xlim[2], yticks = ylim[1]:5:ylim[2], title="Desired J=$J_desired")
nBins = 500
hist!(ax, JT_E[end], bins = nBins, color = :values, strokewidth = 0.0, strokecolor = :black, normalization = :pdf)
lines!(ax, [Point{2,Float64}(J_desired,ylim[1]), Point{2,Float64}(J_desired, ylim[2])], color=:red, linewidth=2, linestyle=:dash)
xlims!(ax, xlim[1], xlim[2])
ylims!(ax, ylim[1], ylim[2])
screen = display(GLMakie.Screen(), fig)
GLMakie.set_title!(screen, "Histogram")