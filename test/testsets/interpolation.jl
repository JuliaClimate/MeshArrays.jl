@testset "Interpolation" begin
    λ=interpolation_setup()
    @test isa(λ,NamedTuple)
    @test all(isfinite.(λ.f))
    # Test LLC90 interpolation coefficients
    γ=GridSpec(ID=:LLC90)
    λ=interpolation_setup(γ)
    @test isa(λ,NamedTuple)
    @test all(isfinite.(λ.f))
    # Test LLC270 interpolation coefficients
    γ=GridSpec(ID=:LLC270)
    λ=interpolation_setup(γ)
    @test isa(λ,NamedTuple)
    @test all(isfinite.(λ.f))
    # Test outside LLC90/LLC270
    γ=GridSpec(ID=:CS32)
    λ=interpolation_setup(γ)
    @test isa(λ,NamedTuple)
    @test all(isfinite.(λ.f))

    Γ=GridLoad(γ;option="light")
    λ=interpolation_setup(Γ=Γ,
        lon=[i for i=-170.:20.0:170., j=-80.:20.0:80.],
        lat=[j for i=-170.:20.0:170., j=-80.:20.0:80.])
    @test isa(λ,NamedTuple)
    @test all(isfinite.(λ.f))
end

