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
