"""
    severusdir()

# Description 

This function simply returns the string for the Severus path. This is helpful for instance to load items, such as meshes, from the `assets`` folder. 
"""
function severusdir()
    joinpath(@__DIR__, "..")
end

function importGeometry(pathInputGeom, fileName_set; pointSpacing = nothing)
    FF = Vector{Vector{TriangleFace{Int}}}()
    VV = Vector{Vector{Point{3,Float64}}}()
    for i in 1:1:length(fileName_set)
        M = load(joinpath(pathInputGeom,fileName_set[i]))
        F = tofaces(faces(M))
        V = topoints(coordinates(M))
        V = [Point{3,Float64}(v) for v in V]
        F,V,_,_ = mergevertices(F,V)

        # Remesh evenly using desired point spacing
        if !isnothing(pointSpacing)
            np = spacing2numvertices(F,V,pointSpacing) # Determine number of points to request to closely match desired point spacing
            F,V = ggremesh(F,V; nb_pts=np) # Remesh with desired number of points 
        end
        # Append current face and vertex set
        push!(FF,F)
        push!(VV,V)
    end
    return FF, VV
end

"""
Function to smooth the boundary 

"""

function smoothboundary(Vs,Fs,Eb; Z_thickness_distal=nothing)
    if isnothing(Z_thickness_distal)
        indCurve_boundary = edges2curve(Eb)
        constrained_points =nothing 
        Vs_smooth = smoothmesh_hc(Fs[indCurve_boundary],Vs, 50, 0.1, 0.5; constrained_points=constrained_points)
    else
        Vs = smoothmesh_laplacian(Eb,Vs, 25, 0.5)   
        ZF = [mean(Vs[f])[3] for f in Fs]
        logicSmooth = ZF.<(minimum(ZF)+Z_thickness_distal)
        constrained_points = elements2indices(boundaryedges(Fs[logicSmooth]))
        Vs_smooth = smoothmesh_hc(Fs[logicSmooth],Vs, 50, 0.1, 0.5; constrained_points=constrained_points)
    end       
    return Vs_smooth
end 

"""
function for adding taper to the skin

"""

function taperfunction(Vs,V_tibia,indCurve_tibia,taperFraction,Z_taper_offset,taperThreshold)
    Zs = [v[3] for v in Vs]
    min_Z = minimum(Zs)
    z_level_taper = min_Z + Z_taper_offset
    # Vt_end = V_tibia[indCurve_tibia]
    D= zeros(length(Vs))
    Dp = zeros(length(Vs))
    indMin1 = Int.(zeros(length(Vs)))
    R = fill(Point{3, Float64}(0.0, 0.0, 0.0), length(Vs))
    for (i,v) in enumerate(Vs)
        if v[3] < z_level_taper
            D[i] = ((v[3] - z_level_taper) / (min_Z-z_level_taper))^2        
            d = Vector{Float64}(undef,length(indCurve_tibia))
            for (i_d,j) in enumerate(indCurve_tibia)
                d[i_d] = sqrt((v[1]-V_tibia[j][1])^2 + (v[2]-V_tibia[j][2])^2)
            end
            dp, index = findmin(d) 
            if dp <= taperThreshold
                Dp[i] = taperThreshold
                indMin1[i] = Int(index)
            else 
                Dp[i] = dp
                indMin1[i] = Int(index)
            end
        
        R[i] = V_tibia[indCurve_tibia[Int(index)]] .- v
        #   push!(R1,r)
        end
    end
    Dp ./= maximum(Dp)
    
    R = R.*D.*Dp.*taperFraction
    #process skin taper
    Vs1 = Vs
    Vs_new = [Point{3,Float64}(v[1]+r[1],v[2]+r[2],v[3]) for (v,r) in zip(Vs1,R)]
    return Vs_new 
end 

"""
Function for bezierclose
"""
function bezierclose(Fs, Vs_new, v1, v2, indCurve; numPointBezier = numPointBezier, p_distal = nothing, numEdgeSmoothSteps = 25, λ = 0.5, pointSpacing=nothing, remesh=true)
    
    if isnothing(pointSpacing)
        pointSpacing = mean(edgelengths(Fs,Vs_new))
    end

    #Bezier
    NV = vertexnormal(Fs,Vs_new; weighting=:size) # May be inefficient to do this to all vertices (/faces)
    NV_indCurve = NV[indCurve] # Vertex normals for the curve
    NE_indCurve = Vector{Point{3,Float64}}(undef,length(indCurve))
    NB_base = Vector{Point{3,Float64}}(undef,length(indCurve))
    m = length(indCurve)
    for i in eachindex(indCurve)
        ii = indCurve[i]
        jj = indCurve[mod1(i+1,m)]
        NE_indCurve[i] = normalizevector(Vs_new[jj] - Vs_new[ii])

        NB_base[i] = cross(NE_indCurve[i],NV_indCurve[i])
    end

    for _ in 1:1:numEdgeSmoothSteps
        for i in eachindex(indCurve)
            i_prev = mod1(i-1,m)
            i_next = mod1(i+1,m)
            NB_base[i] = ((1.0-λ) .* NB_base[i]) .+  λ.*(NB_base[i_prev] .+ NB_base[i_next])./2.0    
            NB_base[i] = normalizevector(NB_base[i])
        end
    end

    VV_B = Vector{Vector{Point{3,Float64}}}()
    Vb = Vector{Point{3,Float64}}()
    
    if isnothing(p_distal)     
        vz = [v[3] for v in Vs_new]
        len = abs(maximum(vz) - minimum(vz))
        v_low = mean(Vs_new[indCurve])
        p_distal = Point{3,Float64}(v_low[1], v_low[2], v_low[3] - len*0.03)        
    end
    
    if !isa(v1, Vector)
        v1 = fill(v1, length(indCurve))
    end

    if !isa(v2, Vector)
        v2 = fill(v2, length(indCurve))
    end

    for i in eachindex(indCurve)
        nBase = NB_base[i]
        nEnd = normalizevector(Point{3,Float64}(Vs_new[indCurve[i]][1]-p_distal[1],Vs_new[indCurve[i]][2]-p_distal[2],0.0))

        pb = [Vs_new[indCurve[i]], Vs_new[indCurve[i]] + v1[i]/3*nBase, p_distal + v2[i]/3*nEnd,p_distal]
        vb = nbezier(pb,numPointBezier)
        push!(VV_B,vb)
        append!(Vb,vb)
    end   

    Fb = grid2surf(Vb,length(indCurve);periodicity=(false,true),face_type=:quad)   
    Fbq = quad2tri(Fb,Vb)   
    F_n,V_n,_ = joingeom(Fs,Vs_new,Fbq,Vb)
    F_n,V_n,_ = mergevertices(F_n,V_n)

    if remesh == true
        # Determine number of points to request to closely match desired point spacing
        np = spacing2numvertices(F_n,V_n,pointSpacing) 
        # Remesh with desired number of points 
        Fn,Vn = ggremesh(F_n,V_n; nb_pts=np)
        return Fn,Vn
    else
        return F_n,V_n
    end
end 

####bones bezier close function 
function close_bones_bezier(FF,VV,nameSet,bezier_bones,velocity,z_distance, pointSpacing; numPointBezier=numPointBezier)
    for model in eachindex(bezier_bones)
         modelName = bezier_bones[model]             
        indexNow = findall(nameSet .== modelName)[1]    
        v_vec = velocity[model]
        z_dist = z_distance[model]
        ind_curve,_= getboundaryset(FF[indexNow],VV[indexNow])
        
        zMax = maximum([v[3] for v in VV[indexNow]])
        Fn,Vn = bezierclose(FF[indexNow],VV[indexNow],v_vec[1],v_vec[2],ind_curve[1]; numPointBezier = numPointBezier,p_distal=z_dist, numEdgeSmoothSteps = 25, λ = 0.5, pointSpacing=pointSpacing)
        edgeSet = boundaryedges(Fn)
        indicesCurve = edges2curve(edgeSet; remove_last = true)
        for i in indicesCurve
            Vn[i] = Point{3,Float64}(Vn[i][1],Vn[i][2],zMax)
        end
       
        FF[indexNow] = Fn
        VV[indexNow] = Vn    
    end 
    return FF,VV
end 

#snap to z z_level

function set_z_level_boundary!(Fn,Vn,z_level;indCurve=nothing)
    #this works only for a single boundary 
    if isnothing(indCurve)
        Eb = boundaryedges(Fn)
        indCurve = edges2curve(Eb)
        indCurve = indCurve[1:end-1]
    end 
    for i in indCurve
        Vn[i] = Point{3,Float64}(Vn[i][1],Vn[i][2],z_level)
    end
    # return Vn 
end 

function getboundaryset(F,V)
    Eb = boundaryedges(F)
    Gb = meshgroup(Eb)
    numGroups = maximum(Gb)
    z_level = zeros(numGroups)
     boundary_set = Vector{Vector{Int}}()
    z = [v[3] for v in V]
    for g in 1:1:numGroups
        ind = edges2curve(Eb[Gb.==g])
        ind = ind[1:end-1]  
        push!(boundary_set,ind)     
        z_level[g] = mean(z[ind])
    end
    indSort = sortperm(z_level, rev=true)
    boundary_set = boundary_set[indSort]
    z_level = z_level[indSort]

    for (i,ind) in enumerate(boundary_set)
        n = edgecrossproduct(V[ind])
        if dot(n, Vec{3, Float64}(0.0, 0.0, 1.0))<0.0
            boundary_set[i] = reverse(ind)
        end
    end

    return boundary_set, z_level
end

function cutEnds!(FF, VV, zCut_levels,n;snapTolerance,pointSpacing)
    for i in eachindex(FF)
        F = FF[i]
        V = VV[i] 
        # Process cuts 
        for (i_cut,z_cut_level) in enumerate(zCut_levels)
            if i_cut == 1
                output_type =:below
            else
                output_type =:above
            end
            if maximum(x-> x[3],V)>z_cut_level && minimum(x-> x[3],V)<z_cut_level # Shape is cutable
                p = Point{3,Float64}(0.0, 0.0, z_cut_level) # Point on cutting plane
                F,V,_ = trisurfslice(F,V,n,p; output_type=output_type, snapTolerance = snapTolerance) 
                np = spacing2numvertices(F,V,pointSpacing)
                F,V = ggremesh(F,V; nb_pts = np)
                boundary_set,z_level = getboundaryset(F,V)
                # println(z_level)
                if length(z_level) == 1
                    set_z_level_boundary!(F,V,z_cut_level;indCurve=boundary_set[1])
                elseif (output_type==:above ) && z_level[2] < abs(mean(v->v[3],V)) && z_level[2] > z_cut_level
                     set_z_level_boundary!(F,V,zCut_levels[1];indCurve=boundary_set[1])
                     set_z_level_boundary!(F,V,zCut_levels[2];indCurve=boundary_set[2])
                end             
                  FF[i] = F
                  VV[i] = V
            end
            # FF[i] = F
            # VV[i] = V 
            # println(i,"###",maximum(v->v[3],V),"###",minimum(v->v[3],V),"###",z_cut_level)
        end 
        
    end
end
#### alternative cut_bones
function cut_bones(FF,VV,nameSet,nameSetCut_down,Z_cut_offset_down,v_mid_patella,n;output_type,snapTolerance)
    for model in eachindex(nameSetCut_down)
        Z_cut_down = Z_cut_offset_down[model]
        modelName = nameSetCut_down[model]
        indexNow = findall(nameSet .== modelName)[1]    
        p = Point{3,Float64}(v_mid_patella[1],v_mid_patella[2],v_mid_patella[3]-Z_cut_down) # Point on cutting plane   
        Fn,Vn,_ = trisurfslice(FF[indexNow],VV[indexNow],n,p; output_type=output_type, snapTolerance = snapTolerance) 
        FF[indexNow] = Fn
        VV[indexNow] = Vn
    end 
    return FF,VV
end     

#
function close_top_bottom(FF,VV,nameSet,bone_set;pointSpacing=pointSpacing)
    indCurve_boneset = []

    for modelName in bone_set 
        indexNow = findall(nameSet.==modelName)[1]
        F = FF[indexNow]
        V = VV[indexNow]
        i_c,_ = getboundaryset(F,V)
        push!(indCurve_boneset,i_c[1])
    end 

    # Close top 
    F_skin = FF[findall(nameSet.=="skin")[1]]
    V_skin = VV[findall(nameSet.=="skin")[1]]
    V_femur = VV[findall(nameSet.=="femur")[1]]
    boundary_set_skin,_ = getboundaryset(F_skin,V_skin)
    VT = (V_skin[boundary_set_skin[1]], V_femur[indCurve_boneset[1]]) # Curves
    R = ([1,2],) 
    P = (pointSpacing,pointSpacing)  # Point spacings
    Ft,Vt,_ = regiontrimesh(VT,R,P)
    # Close bottom   
    V_tibia =  VV[findall(nameSet.=="tibia")[1]]
    V_fibula =  VV[findall(nameSet.=="fibula")[1]]
    VT_b = (V_skin[boundary_set_skin[2]],V_tibia[indCurve_boneset[2]],V_fibula[indCurve_boneset[3]],) # Curves
    R_b = ([1,2,3],)
    P_b = (pointSpacing,pointSpacing,pointSpacing)
    Ft_b,Vt_b,_ = regiontrimesh(VT_b,R_b,P_b) 
    return Ft,Vt,Ft_b,Vt_b
end

#
function volume_tet_ideal(F,V)     
    # indFace = Int(round(length(F)/2))
    # vol_re =  faceinteriorpoint(F,V,indFace;w = 0.5 )
         D = edgelengths(F,V)
        vol = mean(D)^3 / (6.0*sqrt(2.0))
    return vol
end 
#mapping functions 
function mapfunction_D(V_skin,Z_taper_start,Z_taper_end)
    w = Z_taper_start-Z_taper_end
    D = zeros(length(V_skin))
    for (i,v) in enumerate(V_skin)       
        if v[3] <= Z_taper_start && v[3] >= Z_taper_end 
            D[i] = (Z_taper_start - v[3]) / w
        elseif v[3] < Z_taper_end         
            D[i] = 1.0
        end
    end 
    return D 
end 

### geting the spline points 
function bezierMap(t,t_map; n=100)
    P = [Point{2, Float64}(t[i],t_map[i]) for i in eachindex(t)] # Define control points
    V = nbezier(P,n) # Get Bezier fit points
    T_bez = [p[1] for p in V]
    T_map_bez = [p[2] for p in V]
    return T_bez, T_map_bez
end
# map it to the skin 
function map_D(T,T_map, D; spline_order=4)                                 
    S = BSplineKit.interpolate(T, T_map, BSplineOrder(spline_order)) # Create interpolator
    return S.(D)
end

#new shape  
function new_shape_residuum(V_skin,F_skin,D_map,U_vec, v_mid_patella, Z_level_cut)
    V2 = deepcopy(V_skin)
    for i in eachindex(V_skin)
        V2[i] += Point{3,Float64}(U_vec[i] .* D_map[i])
    end

    F2 = deepcopy(F_skin)
    
    p = Point{3,Float64}(v_mid_patella[1],v_mid_patella[2], Z_level_cut) # Point on cutting plane 
    n = normalizevector(Vec{3, Float64}(0.0,0.0,1.0))# Cutting plane normal  
    F2, V2,_ = trisurfslice(F2, V2, n, p; output_type=:above, snapTolerance =  1e-6)     
    boundary_set_s1,_ = getboundaryset(F2,V2)
    set_z_level_boundary!(F2,V2, Z_level_cut; indCurve = boundary_set_s1[2])
    return F2, V2
end 

#FEBio
function isotropictraction(E,ν,J)
    # Lamé parameters
    μ = E./(2*(1+ν))
    λ = (ν*E)/((1+ν)*(1-2*ν))

    # Hydrostatic stress
    σₕ = (μ/J) * (J^(2.0/3.0)-1.0) + (λ/J) * log(J)

    # Required contraction given J and material parameters 
    T = -σₕ*J^(1.0/3.0)
    return T
end

function shrinky(exe_path,FF,VV,E,V,CE,Fb,Cb,κ,ratio,J_desired; endCapLabels=[6,7])
    ######
    # Set FEBio exec path or name
    # const FEBIO_EXEC = "febio4" # FEBio executable
     FEBIO_EXEC = exe_path # Path to FEBio executable

    ###### 
    # Control parameters 
    μ = κ/ratio 
    E_youngs = 9.0*κ*μ/(3.0*κ + μ)
    ν = (3.0*κ - 2.0 *μ)/(2.0*(3.0*κ + μ))
    T0 = isotropictraction(E_youngs,ν,J_desired)
    
    # Node sets definitions for end caps and bones 
    boolFacesEndCaps = [in(c,endCapLabels) for c in Cb]
    indNodesEndCaps = unique(reduce(vcat,Fb[boolFacesEndCaps]))
    indNodesBones = unique(reduce(vcat,Fb[[in(c,[2,3,4,5]) for c in Cb]]))

    ######
    # Define file names
    saveDir = joinpath(febiojl_dir(),"assets","temp") # Main directory to save FEBio input and output files
    if !isdir(saveDir)
        mkdir(saveDir)      
    end

    filename_FEB = joinpath(saveDir,"febioInputFile_01.feb")   # The main FEBio input file
    filename_xplt = joinpath(saveDir,"febioInputFile_01.xplt") # The XPLT file for viewing results in FEBioStudio
    filename_log = joinpath(saveDir,"febioInputFile_01_LOG.txt") # The log file featuring the full FEBio terminal output stream
    filename_disp = "febioInputFile_01_DISP.txt" # A log file for results saved in same directory as .feb file  e.g. nodal displacements
    filename_stress = "febioInputFile_01_STRESS.txt"
    filename_J = "febioInputFile_01_J.txt"

    ######
    # Define febio input file XML
    doc,febio_spec_node = feb_doc_initialize()

    aen(febio_spec_node,"Module"; type = "solid") # Define Module node: <Module type="solid"/>

    time_steps = 10 
    step_size = 1.0/time_steps

    control_node = aen(febio_spec_node,"Control") # Define Control node: <Control>
        aen(control_node,"analysis","STATIC")                
        aen(control_node,"time_steps",time_steps)
        aen(control_node,"step_size",step_size)
        aen(control_node,"plot_zero_state",1)
        aen(control_node,"plot_range",@sprintf("%.2f, %.2f",0,-1))
        aen(control_node,"plot_level","PLOT_MAJOR_ITRS")
        aen(control_node,"plot_stride",1)
        aen(control_node,"output_level","OUTPUT_MAJOR_ITRS")
        aen(control_node,"adaptor_re_solve",1)

    time_stepper_node = aen(control_node,"time_stepper"; type = "default")
        aen(time_stepper_node,"max_retries",5)
        aen(time_stepper_node,"opt_iter",10)
        aen(time_stepper_node,"dtmin",1e-3)
        aen(time_stepper_node,"dtmax",step_size)
        aen(time_stepper_node,"aggressiveness",0)
        aen(time_stepper_node,"cutback",5e-1)
        aen(time_stepper_node,"dtforce",0)

    solver_node = aen(control_node,"solver"; type = "solid")
        aen(solver_node,"symmetric_stiffness",1)
        aen(solver_node,"equation_scheme",1)
        aen(solver_node,"equation_order","default")
        aen(solver_node,"optimize_bw",0)
        aen(solver_node,"lstol",9e-1)
        aen(solver_node,"lsmin",1e-2)
        aen(solver_node,"lsiter",5)
        aen(solver_node,"max_refs",70)
        aen(solver_node,"check_zero_diagonal",0)
        aen(solver_node,"zero_diagonal_tol",0)
        aen(solver_node,"force_partition",0)
        aen(solver_node,"reform_each_time_step",1)
        aen(solver_node,"reform_augment",0)
        aen(solver_node,"diverge_reform",1)
        aen(solver_node,"min_residual",1e-20)
        aen(solver_node,"max_residual",0)
        aen(solver_node,"dtol",1e-3)
        aen(solver_node,"etol",1e-2)
        aen(solver_node,"rtol",0)
        aen(solver_node,"rhoi",0)
        aen(solver_node,"alpha",1)
        aen(solver_node,"beta",2.5e-01)
        aen(solver_node,"gamma",5e-01)
        aen(solver_node,"logSolve",0)
        aen(solver_node,"arc_length",0)
        aen(solver_node,"arc_length_scale",0)
    qn_method_node = aen(solver_node,"qn_method"; type = "BFGS")
        aen(qn_method_node,"max_ups",0)
        aen(qn_method_node,"max_buffer_size",0)
        aen(qn_method_node,"cycle_buffer",0)
        aen(qn_method_node,"cmax",0)

    Globals_node   = aen(febio_spec_node,"Globals")

    Constants_node = aen(Globals_node,"Constants")
        aen(Constants_node,"R",8.3140000e-06)
        aen(Constants_node,"T",298)
        aen(Constants_node,"F",9.6485000e-05)

    Material_node = aen(febio_spec_node,"Material")

    material_node = aen(Material_node,"material"; id = 1, name="Material1", type="solid mixture")    
            solid_node = aen(material_node,"solid"; type="neo-Hookean")                    
                            aen(solid_node,"E",E_youngs)
                            aen(solid_node,"v",ν)

            solid_node = aen(material_node,"solid"; type="prescribed isotropic active contraction")
                            aen(solid_node,"T0",@sprintf("%.16e",T0); lc=1)
                            
    Mesh_node = aen(febio_spec_node,"Mesh")

    # Nodes
    Nodes_node = aen(Mesh_node,"Nodes"; name="nodeSet_all")
        for (i,v) in enumerate(V)        
            aen(Nodes_node,"node", join([@sprintf("%.16e",x) for x ∈ v],','); id = i)     
        end
        
    # Elements
    Elements_node = aen(Mesh_node,"Elements"; name="Part1", type="tet4")
        for (i,e) in enumerate(E)        
            aen(Elements_node,"elem", join([@sprintf("%i",i) for i ∈ e],','); id = i)
            # aen(Elements_node,"elem",join(map(string, e), ','); id = i)
        end
        
    # Node sets
    bcSupportList = "bcSupportList"
    bcSupportList_z = "bcSupportList_z"
    aen(Mesh_node,"NodeSet",join([@sprintf("%i", i) for i ∈ indNodesBones],','); name=bcSupportList)
    aen(Mesh_node,"NodeSet",join([@sprintf("%i", i) for i ∈ indNodesEndCaps],','); name=bcSupportList_z)

    MeshDomains_node = aen(febio_spec_node, "MeshDomains")
        aen(MeshDomains_node,"SolidDomain"; mat = "Material1", name="Part1")

    Boundary_node = aen(febio_spec_node, "Boundary")

    bc_node = aen(Boundary_node,"bc"; name="zero_displacement", node_set=bcSupportList, type="zero displacement")
        aen(bc_node,"x_dof",1)
        aen(bc_node,"y_dof",1)
        aen(bc_node,"z_dof",1)

    bc_node = aen(Boundary_node,"bc"; name="zero_displacement_z", node_set=bcSupportList_z, type="zero displacement")
        aen(bc_node,"x_dof",0)
        aen(bc_node,"y_dof",0)
        aen(bc_node,"z_dof",1)
        
    LoadData_node = aen(febio_spec_node,"LoadData")

    load_controller_node = aen(LoadData_node,"load_controller"; id=1, name="LC_1", type="loadcurve")
        aen(load_controller_node,"interpolate","LINEAR")
        
    points_node = aen(load_controller_node,"points")
        aen(points_node,"pt",@sprintf("%.2f, %.2f",0.0,0.0))
        aen(points_node,"pt",@sprintf("%.2f, %.2f",1.0,1.0))


    Output_node = aen(febio_spec_node,"Output")

    plotfile_node = aen(Output_node,"plotfile"; type="febio")
        aen(plotfile_node,"var"; type="displacement")
        aen(plotfile_node,"var"; type="stress")
        aen(plotfile_node,"var"; type="relative volume")
        aen(plotfile_node,"var"; type="reaction forces")
        aen(plotfile_node,"compression",@sprintf("%i",0))

    logfile_node = aen(Output_node,"logfile"; file=filename_log)
        aen(logfile_node,"node_data"; data="ux;uy;uz", delim=",", file=filename_disp)
        aen(logfile_node,"element_data"; data="s1;s2;s3", delim=",", file=filename_stress)
        aen(logfile_node,"element_data", join([@sprintf("%i",i) for i ∈ 1:length(E)], ", "); data="J", delim=",", file=filename_J)

    #######
    # Write FEB file
    XML.write(filename_FEB, doc)

    #######
    # Run FEBio
    run_febio(filename_FEB,FEBIO_EXEC)

    #######
    # Import results
    DD_disp = read_logfile(joinpath(saveDir,filename_disp))
    DD_stress = read_logfile(joinpath(saveDir,filename_stress))
    DD_J = read_logfile(joinpath(saveDir,filename_J))
    numInc = length(DD_disp)
    incRange = 0:1:numInc-1
    J_end = DD_J[numInc-1].data

    # Create time varying vectors
    UT = fill(V,numInc) 
    VT = fill(V,numInc)
    J_E = fill(zeros(length(E)),numInc)
    UT_mag = fill(zeros(length(V)),numInc)
    ut_mag_max = zeros(numInc)
    @inbounds for i in 0:1:numInc-1    
        UT[i+1] = [Point{3,Float64}(u) for u in DD_disp[i].data]
        VT[i+1] += UT[i+1]
        UT_mag[i+1] = norm.(UT[i+1])
        J_E[i+1] = [J[1] for J in DD_J[i].data]
        ut_mag_max[i+1] = maximum(UT_mag[i+1]) 
    end

    min_p = minp([minp(V) for V in VT])
    max_p = maxp([maxp(V) for V in VT])
    V_skin = VV[1]
    F_skin = FF[1]
    V_shrunk = VT[end]
    U_vec = UT[end]
    U_mag = UT_mag[end]

    # indNodesSkin = unique(reduce(vcat,Fb[Cb.==1]))
    indUsed = elements2indices(Fb[Cb.==1])
    U_vec = U_vec[indUsed]
    # F_shrunk = Fb[Cb.==1]
    # V_shrunk = V_shrunk[indNodesSkin]
    U_mag = U_mag[indUsed]
    # U_vec = U_vec[indNodesSkin]    
    F_shrunk,V_shrunk = remove_unused_vertices(Fb[Cb.==1],VT[end])

    return UT,VT,FF,VV,F_shrunk,V_shrunk,U_mag,U_vec, J_E, DD_J
end 
