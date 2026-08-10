γ = GridSpec(ID=:LLC90)
Γ = GridLoad(γ, option="full")
u = MeshArray(γ, Float32, 10)
for I in eachindex(u.f)
    u.f[I] = rand(Float32, size(u.f[I])...)
end

@testset "broadcast .= equivalence" begin
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

@testset "broadcast slice .= equivalence" begin
    # A[:,k] .= B[:,k] .* C  — the Drifters.jl pattern
    u_eq  = similar(u, allocate=true)
    u_deq = similar(u, allocate=true)

    for k in 1:10
        u_eq[:,k]  = u[:,k] .* Γ.DXC   # plain =
        u_deq[:,k] .= u[:,k] .* Γ.DXC  # .= should match
    end

    @test all(I -> u_eq.f[I] == u_deq.f[I], eachindex(u_eq.f))
end