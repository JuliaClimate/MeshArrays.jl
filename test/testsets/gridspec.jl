@testset "GridSpec:" begin
    γ = GridSpec(ID=:Oscar)
    γ = GridSpec(ID=:IAP)
    GridLoad(GridSpec(ID=:OISST))
    GridLoad(GridSpec("ones"))
end
