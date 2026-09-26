@testset "Vertical Dimension:" begin
    γ=GridSpec(ID=:onedegree)
    Γ=GridLoad(γ;option="full")
    θ=Float64.(Γ.hFacC)
    nk=length(Γ.RC)
    [θ[:,k]=0.01*(nk-k) .+ cosd.(Γ.YC[:]) for k in 1:nk]
    θ[findall(Γ.hFacC.==0.0)].=NaN
    d=isosurface(θ,1.1,Γ)
    @test isapprox(d[1][180,90],-2204.8384919029777)

    mc=MeshArrays.coldlayer(θ,1.1,Γ)
    mh=MeshArrays.hotlayer(θ,1.1,Γ);
    ml=MeshArrays.layerfraction(θ,1.1,1.2,Γ);
    @test isapprox(sum(mc),835263.7654388269)
    @test isapprox(sum(mh),693042.234561173)
    @test isapprox(sum(ml),242277.4251312976)

end
