using Test, Documenter
using MeshArrays

MeshArrays.GRID_LL360_download()
MeshArrays.GRID_LLC90_download()
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
        γ,Γ=MeshArrays.GridOfOnes(grTopo,nFaces,N;option="full")
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
        @test isa(dRdx_e,MeshArray)
    end
end

@testset "Vertical Dimension:" begin
    γ=GridSpec("PeriodicChannel",MeshArrays.GRID_LL360)
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

@testset "Transport computations:" begin
    #Load grid and transport / vector field
    γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
    Tx=γ.read(MeshArrays.GRID_LLC90*"TrspX.bin",MeshArray(γ,Float32))
    Ty=γ.read(MeshArrays.GRID_LLC90*"TrspY.bin",MeshArray(γ,Float32))
    Γ=GridLoad(γ;option="full")
    
    μ =Γ.hFacC[:,1]
    μ[findall(μ.>0.0)].=1.0
    μ[findall(μ.==0.0)].=NaN

    lons=[-68 -63]; lats=[-54 -66]; name="Drake Passage"
    Trsct=Transect(name,lons,lats,Γ)

    #Various vector operations
    U=0*Γ.hFacW; V=0*Γ.hFacS;
    UVtoTransport(U,V,Γ)
    UVtoUEVN(U[:,1],V[:,1],Γ)
    curl(U[:,1],V[:,1],Γ)
    
    #Meridional transport integral
    uv=Dict("U"=>Tx,"V"=>Ty,"dimensions"=>["x","y"])
    L=-85.0:5.0:85.0; LC=LatitudeCircles(L,Γ)
    T=Array{Float64,1}(undef,length(LC))
    [T[i]=1e-6*ThroughFlow(uv,LC[i],Γ) for i=1:length(LC)]

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
    TrspPsi=VectorPotential(TxR,TyR,Γ);

    tst=extrema(write(TrspPsi))

    @test prod(isapprox.(tst,(-8.08f7,2.19f8)))
end

@testset "gcmfaces type:" begin
    γ=GridSpec("CubeSphere",MeshArrays.GRID_CS32)

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
    show(tmp)

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
    show(tmp1)
    similar(tmp1)
    #tmp[tmp1]

    show(tmp)
    MeshArrays.getindexetc(tmp,2)
    MeshArray(γ,tmp.f,meta=tmp.meta)
    MeshArray(γ,meta=tmp.meta)
    MeshArray(γ,Float32,3,4)
    MeshArray(γ,Float32,tmp.fSize,tmp.fIndex,2,3)

    tmp1=MeshArray(γ,Float32,3)
    MeshArray(γ,tmp1.f,meta=tmp1.meta)
    MeshArrays.getindexetc(tmp1,2,1)
end

@testset "UnitGrid:" begin
    (Γ,γ)=UnitGrid( (80,90) , (20,30) ; option="full")
    @test isa(γ,gcmgrid)

    tmp=UnitGrid(γ)
    @test isa(tmp,NamedTuple)

    #various read/write functions
    read(write(Γ.XC),γ)
    tmp=tempname()
    write(tmp,Γ.XC)
    MeshArrays.read_tiles(tmp,Γ.XC)
    MeshArrays.write_tiles(Γ.XC)
    MeshArrays.write_tiles(tmp,Γ.XC)    
    @test isfile(tmp)
end

@testset "doctests" begin
    doctest(MeshArrays; manual = false)
end
