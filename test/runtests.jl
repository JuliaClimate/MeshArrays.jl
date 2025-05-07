using Test, Documenter, Suppressor, MeshArrays, CairoMakie
import DataDeps, JLD2, Shapefile, GeoJSON, Proj

MeshArrays.GRID_LL360_download()
MeshArrays.GRID_LLC90_download()
MeshArrays.GRID_LLC270_download()
MeshArrays.GRID_CS32_download()

p=dirname(pathof(MeshArrays))
include(joinpath(p,"../examples/Demos.jl"))

@testset "MeshArrays tests:" begin
    for nTopo=1:4
        if nTopo==1; grTopo="CubeSphere"; nFaces=6; N=200;
        elseif nTopo==2; grTopo="LatLonCap"; nFaces=5; N=200;
        elseif nTopo==3; grTopo="PeriodicChannel"; nFaces=1; N=400;
        elseif nTopo==4; grTopo="PeriodicDomain"; nFaces=1; N=400;
        end;
        Npt=nFaces*N*N
        γ,Γ=Grids_simple.GridOfOnes(grTopo,nFaces,N;option="full")
        @test γ.class == grTopo
        Rini= 0.; Rend= 0.;
        (Rini,Rend,DXCsm,DYCsm)=demo2(Γ);
        @test isa(Rend,MeshArray)
        @test sum(isfinite.(Rend)) == Npt
        Sini=sqrt(sum(Rini*Rini)/(Npt-1.0))
        Send=sqrt(sum(Rend*Rend)/(Npt-1.0))
        #println([Sini Send])
        @test isapprox(Sini,1.000; atol=1e-2)
        @test isapprox(Send,0.093; atol=1e-2)
        (dRdx,dRdy)=gradient(Rend,Γ)
        exchange(Rend,2)
        exchange(dRdx,dRdy,1)        
        (dRdx_e,dRdy_e)=exchange(dRdx,dRdy,1)        
        @test isa(dRdx_e,MeshArray_wh)
    end
end

@testset "Vertical Dimension:" begin
    γ=GridSpec(ID=:onedegree)
    Γ=GridLoad(γ;option="full")
    θ=Float64.(Γ.hFacC)
    nk=length(Γ.RC)
    #[θ[:,k]=cosd.((nk-k)/nk*Γ.YC) for k in 1:nk]
    #why fail? [θ[:,k]=(nk-k) .+ cosd.(0.5*Γ.YC[:]) for k in 1:nk]
    [θ[1,k]=0.01*(nk-k) .+ cosd.(Γ.YC[1]) for k in 1:nk]
    θ[findall(Γ.hFacC.==0.0)].=NaN
    d=isosurface(θ,1.1,Γ)

    @test isapprox(d[1][180,90],-2204.8384919029777)
end

@testset "Regional Integration:" begin
    G,M,files=Integration.example()
    @suppress show(M)

    allones=1.0 .+0*G.hFacC
    vol0=sum(G.RAC*G.DRF*G.hFacC)

    vol=Integration.volumes(M,G)
    @test isapprox(sum(vol),vol0)

    G,M,files=Integration.example(option=:streamlined_loop)
    vol=[b(allones) for b in M.h_sum]
    @test isapprox(sum(vol),vol0)

    rgns=Integration.define_regions(option=:basins,grid=G)
    rgns=Integration.define_regions(option=:dlat_10,grid=G)
    rgns=Integration.define_regions(option=(30,10),grid=G)
    @test isa(rgns,NamedTuple)

    G,M,files=Integration.example()
    files=fill("?",3)
    rd0(F,var,tim,tmp)=tim*ones(tmp)
    H=Integration.loops(M,files=files,rd=rd0)
    @test isa(H,Array)

    G,M,files=Integration.example(option=:streamlined_loop)
    files=fill("?",3)
    H=Integration.streamlined_loop(M,files=files,rd=rd0)
    @test isa(H,Array)
end

@testset "Transport computations:" begin
    #Load grid and transport / vector field
    γ=GridSpec(ID=:LLC90)
    @suppress show(γ)
    Γ=GridLoad(ID=:LLC90,option=:full)
    @suppress show(Γ.XC)

    Tx=γ.read(MeshArrays.GRID_LLC90*"TrspX.bin",MeshArray(γ,Float32))
    Ty=γ.read(MeshArrays.GRID_LLC90*"TrspY.bin",MeshArray(γ,Float32))
    plot(Γ.XC)

    hFacC=GridLoadVar("hFacC",γ)
    μ=land_mask(hFacC[:,1])

    lons=[-68 -63]; lats=[-54 -66]; name="Drake Passage"
    Trsct=Transect(name,lons,lats,Γ,segment=:long,format=:NamedTuple)
    Trsct=Transect(name,lons,lats,Γ)

    mask=demo.extended_basin(demo.ocean_basins(),:Pac)
    edge=edge_path("Pacific Ocean Edge",mask,Γ)

    #Various vector operations
    hFacW=GridLoadVar("hFacW",γ)
    hFacS=GridLoadVar("hFacS",γ)
    RAZ=GridLoadVar("RAZ",γ)

    U=0*hFacW; V=0*hFacS;
    UVtoTransport(U,V,Γ)
    UVtoUEVN(U[:,1],V[:,1],Γ)
    curl(U[:,1],V[:,1], merge(Γ,(hFacW=hFacW,hFacS=hFacS,RAZ=RAZ)) )
    dD=zeros(γ)
    MeshArrays.UVtoSpeed!(U[:,1],V[:,1],Γ,dD)
    
    #Meridional transport integral
    uv=Dict("U"=>Tx,"V"=>Ty,"dimensions"=>["x","y"])
    L=-85.0:5.0:85.0; LC=LatitudeCircles(L,Γ,format=:gridpath)
    T=Array{Float64,1}(undef,length(LC))
    [T[i]=1e-6*ThroughFlow(uv,LC[i],Γ) for i=1:length(LC)]
    plot(LC)
    plot(LC[1])

    x=zeros(γ)
    fill!(x,1.0)
    y=fill(-1.0,γ)
    extrema(y)
    @test minimum(y)<minimum(x)

    y*ones(3,2)
    ones(γ)
    ones(y)
    zeros(y)

    GM_PsiX=read(randn(90,1170,50),Γ.hFacW)
    GM_PsiY=read(randn(90,1170,50),Γ.hFacS)
    bolusU, bolusV, bolusW=MeshArrays.calc_bolus(GM_PsiX,GM_PsiY, Γ)
    
    read(rand(90*1170),γ)
    read(rand(90*1170*2),γ)
    read(rand(90,1170,2,2),γ)

    #See: OceanTransports/helper_functions.jl
    #u,v,uC,vC=rotate_uv(uv,Γ);

    #Transport convergence   
    TrspCon=μ.*convergence(Tx,Ty)
    #Scalar potential
    TrspPot=ScalarPotential(TrspCon)
    #Divergent transport component
    (TxD,TyD)=gradient(TrspPot,Γ)
    TxD=TxD.*Γ.DXC
    TyD=TyD.*Γ.DYC
    #Rotational transport component
    TxR = Tx-TxD
    TyR = Ty-TyD
    #vector Potential
    TrspPsi=VectorPotential(TxR,TyR, merge(Γ,(hFacC=hFacC,hFacW=hFacW,hFacS=hFacS)) );

    tst=extrema(write(TrspPsi))

    @test prod(isapprox.(tst,(-8.08f7,2.19f8)))

    #Ekman transport
    de4=demo4()
    @test isapprox(de4.Tr[120],1.317456)

end

@testset "gcmfaces type:" begin
    for ID in (:CS32, :LLC270)
        γ=GridSpec(ID=ID)
        MeshArrays.gcmfaces(γ)
        MeshArrays.gcmfaces(γ,Float32)
        tmp=MeshArrays.gcmfaces(γ,Float32,3)
        tmp[3,1,2]=1.0
        view(tmp,:,1,1:2)
        MeshArrays.gcmfaces(γ,tmp.f)
        MeshArrays.fsize(tmp)
        MeshArrays.fsize(tmp,2)
        MeshArrays.fsize(tmp.f)
        MeshArrays.fsize(tmp.f,2)
        size(tmp)
        size(tmp,3)
        tmp1=similar(tmp)
        2 .*tmp1
        findall(tmp.>0)
        MeshArrays.nFacesEtc(tmp)
        @suppress show(tmp)

        x=tmp[1:10,1,1:2]; y=x[2]; x[3]=1.0
        view(x,1:3,:,1)
        MeshArrays.gcmsubset(γ,x.f,x.fSize,x.aSize,x.i,x.iSize)
        MeshArrays.fsize(x.f)
        MeshArrays.fsize(x.f,2)
        size(x)
        size(x,3)

        MeshArrays.fijind(tmp,10)
        MeshArrays.fsize(tmp)
        @test isa(tmp,MeshArrays.gcmfaces)

        tmp=MeshArray(γ)
        tmp1=findall(tmp.>0)
        tmp[tmp1]
        tmp[tmp1].=1.0
        size(tmp1)
        tmp1[2]
        view(tmp1,:)
        @suppress show(tmp1)
        similar(tmp1)
        #tmp[tmp1]

        @suppress show(tmp)
        MeshArrays.getindexetc(tmp,2)
        MeshArray(γ,tmp.f,meta=tmp.meta)
        MeshArray(γ,meta=tmp.meta)
        MeshArray(γ,Float32,3,4)
        MeshArray(γ,Float32,tmp.fSize,tmp.fIndex,2,3)

        tmp1=MeshArray(γ,Float32,3)
        MeshArray(γ,tmp1.f,meta=tmp1.meta)
        MeshArrays.getindexetc(tmp1,2,1)
    end
end

@testset "UnitGrid:" begin
    C=MeshArray(randn(20,10))
    D=MeshArray(randn(20,10,3))

    (Γ,γ)=Grids_simple.UnitGrid( (80,90) , (20,30) ; option="full")
    @test isa(γ,gcmgrid)

    tmp=Grids_simple.UnitGrid(γ)
    @test isa(tmp,NamedTuple)

    #various read/write functions
    read(write(Γ.XC),γ)
    tmp=tempname()
    write(tmp,Γ.XC)
    MeshArrays.read_tiles(tmp,Γ.XC)
    MeshArrays.write_tiles(Γ.XC)
    MeshArrays.write_tiles(tmp,Γ.XC)    
    @test isfile(tmp)

    xy=Grids_simple.xy_OISST()
    xy=Grids_simple.xy_Oscar()

    xy=Grids_simple.xy_IAP()
    gr=Grids_simple.grid_factors(xy)
    dep=[10 100 1000]; msk=ones(gr[:XC].fSize[1]...,3)
    gr=Grids_simple.grid_add_z(gr,dep,msk)

    @test haskey(gr,:hFacC)
end

@testset "Plotting:" begin
    γ=GridSpec(ID=:LLC90)
    Γ=GridLoad(γ;option="light")
    D=Γ.Depth
    λ=interpolation_setup(Γ=Γ,
        lon=[i for i=-170.:20.0:170., j=-80.:20.0:80.], 
        lat=[j for i=-170.:20.0:170., j=-80.:20.0:80.])
    λ=interpolation_setup()

    basins=demo.ocean_basins()
    AtlExt=demo.extended_basin(basins,:Atl)
    sections,path_sec=demo.ocean_sections(Γ)
    my_section=demo.one_section(Γ,[127 127],[-25 -68])

    fig=MeshArrays.plot_examples(:smoothing_demo,D,D)
    (fig1,fig2,fig3)=MeshArrays.plot_examples(:interpolation_demo,Γ)

    MeshArrays.plot_examples(:meriodional_overturning,Γ,rand(179,50))
    MeshArrays.plot_examples(:northward_transport,rand(179))

    MeshArrays.plot_examples(:gradient_EN,λ,D,D)
    MeshArrays.plot_examples(:gradient_xy,λ,D,D)

    ## more methods

    scatter(Γ.XC,Γ.YC,color=:black)
    heatmap(D,interpolation=λ)
    scatter!(current_axis(),Γ.XC,Γ.YC,color=:red)

    heatmap(D) #will display tile by tile
    heatmap(D,interpolation=λ,title="ocean depth") #same but w title

    lon0=-160
    proj=Proj.Transformation(MA_preset=2,lon0=lon0)
    Dint=reshape(Interpolate(D,λ.f,λ.i,λ.j,λ.w),size(λ.lon))

    Interpolate(D,λ)
    InterpolationFactors(Γ,30.0,30.0)

    ###

    fil=demo.download_polygons("countries.geojson")
    pol=MeshArrays.read_polygons(fil)

    MeshArraysMakieExt = Base.get_extension(MeshArrays, :MeshArraysMakieExt)
    dest="+proj=eqearth +lon_0=$(lon0) +lat_1=0.0 +x_0=0.0 +y_0=0.0 +ellps=GRS80"
    MeshArraysMakieExt.split(pol,dest)
    MeshArraysMakieExt.split(pol,Observable(dest))
    MeshArraysMakieExt.split(Observable(pol),Observable(dest))
    MeshArraysMakieExt.split(Observable(pol),dest)

    ###

    fil=demo.download_polygons("ne_110m_admin_0_countries.shp")
    pol=MeshArrays.read_polygons(fil)

    f = Figure()
    ax = f[1, 1] = Axis(f, aspect = DataAspect(), title = "Ocean Depth (m)")
	pr_ax=MeshArrays.ProjAxis(ax; proj=proj,lon0=lon0)
	surf = surface!(pr_ax,λ.lon,λ.lat,0*λ.lat; color=Dint, 
			colorrange=(0.0,6000.0), colormap=:berlin, shading = NoShading)
	lines!(pr_ax; polygons=pol,color=:black,linewidth=0.5)
	MeshArrays.grid_lines!(pr_ax;color=:lightgreen,linewidth=0.5)
	f

	meta=(colorrange=(0.0,6000.0),cmap=:BrBG_10,ttl="Ocean Depth (m)",lon0=lon0)
	data=(lon=λ.lon,lat=λ.lat,var=Dint,meta=meta) #,polygons=pol)
    plot_examples(:projmap,data,lon0,proj)
    plot_examples(:simple_heatmap,data)

end

@testset "nanmath" begin
    x=[NaN 1 2]
    nansum(x)
    nansum(x,1)
    nanmax(x,2)
    nanmin(x,2)

    nanmean(NaN,1)
    nanmean(1,NaN)
    nanmean(2,1)
    nanmean(NaN,NaN)
end

@testset "plotting" begin
    lon,lat,earth_img=demo.get_basemap()
    plot_examples(:basemap,lon,lat,earth_img)

    pol_file=demo.download_polygons("ne_110m_admin_0_countries.shp")
    pol=MeshArrays.read_polygons(pol_file)

    lon0=-160
    proj=Proj.Transformation(MA_preset=2,lon0=lon0)
    plot_examples(:baseproj,proj,lon0,pol=pol)
end

@testset "doctests" begin
    doctest(MeshArrays; manual = false)
end
