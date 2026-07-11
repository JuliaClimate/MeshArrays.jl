@testset "Transport computations:" begin
    #Load grid and transport / vector field
    γ=GridSpec(ID=:LLC90)
    @suppress show(γ)
    Γ=GridLoad(γ,option=:full)
    @suppress show(Γ.XC)

    path=MeshArrays.Dataset("GRID_LLC90")
    Tx=γ.read(joinpath(path,"TrspX.bin"),MeshArray(γ,Float32))
    Ty=γ.read(joinpath(path,"TrspY.bin"),MeshArray(γ,Float32))
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

    LC_nt=LatitudeCircles(L,Γ,format=:NamedTuple)
    @test isa(LC_nt,Array) && isa(LC_nt[1],NamedTuple)
    @test LC_nt[1].lat == L[1]
    @test isapprox(ThroughFlow(uv,LC_nt[1],Γ), ThroughFlow(uv,LC[1],Γ))

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
    @test isapprox(de4.Tr[120],1.318311; atol=1e-5)

end
