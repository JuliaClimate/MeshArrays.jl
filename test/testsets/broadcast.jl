@testset "broadcast .= equivalence" begin
    γ = GridSpec(ID=:LLC90)
    Γ = GridLoad(γ, option="full")
    u = MeshArray(γ, Float32, 10)
    for I in eachindex(u.f)
        u.f[I] = rand(Float32, size(u.f[I])...)
    end

    # Test 1: A .= B  (plain copy)
    u2 = similar(u, allocate=true)
    u2 .= u
    @test all(I -> u2.f[I] == u.f[I], eachindex(u.f))

    # Test 2: A .= scalar .* B
    u3 = similar(u, allocate=true)
    u3 .= 2.0f0 .* u
    @test all(I -> u3.f[I] == 2.0f0 .* u.f[I], eachindex(u.f))

    # Test 3: A .= B .* C (the Drifters.jl pattern, whole array)
    u4 = similar(u, allocate=true)
    ref = similar(u, allocate=true)
    for I in eachindex(u.f)
        ref.f[I] = u.f[I] .* 2.0f0
    end
    u4 .= u .* 2.0f0
    @test all(I -> u4.f[I] == ref.f[I], eachindex(u.f))
end