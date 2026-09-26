"""
    isosurface(θ, T, Γ)

Return the uppermost in-water depth at which `θ` crosses the threshold `T`.

Crossing depths are linearly interpolated between adjacent wet vertical
cell centers, using `Γ.RC`. A crossing is considered only where both
adjacent cells have `Γ.hFacC > 0`; a wet/dry transition is interpreted
as the sea floor and does not constitute an isosurface.

The returned field is `NaN` where no in-water crossing is found,
including land points, below-bottom points, columns entirely above or
below `T`, and columns with only one wet vertical level.

The function scans downward and retains only the first valid crossing,
i.e. the uppermost crossing. It supports non-monotonic vertical profiles,
which may contain multiple crossings.

Dry-cell values in `θ` are ignored. Values of `θ` must be finite in wet
cells.
"""
function isosurface(θ, T, Γ)
    # A 2-D MeshArray, initialized to NaN.
    # NaN means: no in-water isotherm has yet been found.
    d = zeros(θ[:, 1])
    d .= NaN
    nr = size(θ, 2)

    for j in 1:size(d, 1)
        # k = nr is intentionally excluded: there is no k+1 center
        # with which to interpolate an interior θ=T crossing.
        for k in 1:nr-1
            θk  = θ[j, k]
            θk1 = θ[j, k+1]
            wetk  = Γ.hFacC[j, k]   .> 0
            wetk1 = Γ.hFacC[j, k+1] .> 0

            # A crossing can only be interpolated between two wet centers.
            # In particular, wet / dry is the sea floor, not an isotherm.
            pairwet = wetk .& wetk1
            # Optional safety guard. It is harmless if θ is finite everywhere,
            # and prevents accidental interpolation through missing wet data.
            validθ = isfinite.(θk) .& isfinite.(θk1)
            # θ goes from >T at k to <=T at k+1.
            downcross = (θk .> T) .& (θk1 .<= T)
            # θ goes from <=T at k to >T at k+1.
            upcross = (θk .<= T) .& (θk1 .> T)
            # Since k increases downward, the first assignment is the
            # uppermost valid crossing in the column.
            i = findall(
                isnan.(d[j]) .&
                pairwet .&
                validθ .&
                (downcross .| upcross),
            )

            if !isempty(i)
                # Fraction of the interval from RC[k] to RC[k+1].
                # This formula works for either crossing direction.
                a = (θk[i] .- T) ./ (θk[i] .- θk1[i])
                d[j][i] .= (1 .- a) .* Γ.RC[k] .+
                            a .* Γ.RC[k+1]
            end
        end
    end

    return d
end

"""
    layerfraction(θ, Tlo, Thi, Γ)

Return the fraction of each wet tracer cell satisfying

    Tlo ≤ θ < Thi.

The result is computed as

    hotlayer(θ, Tlo, Γ) - hotlayer(θ, Thi, Γ),

and thus inherits the vertical interpolation, wet/dry, surface, and
sea-floor conventions of [`hotlayer`](@ref).

`Tlo` must not exceed `Thi`. The returned fraction refers to the nominal
tracer cell; multiply by `Γ.hFacC` and cell volume for volume-weighted
integrals.

Special cases include:

    coldlayer(θ, T, Γ) == layerfraction(θ, -Inf, T, Γ)
    hotlayer(θ, T, Γ)  == layerfraction(θ, T, Inf, Γ)

up to floating-point roundoff.

See also: [`hotlayer`](@ref), [`coldlayer`](@ref), [`isosurface`](@ref).
"""
function layerfraction(θ, Tlo, Thi, Γ)
    Tlo <= Thi || throw(ArgumentError("require Tlo ≤ Thi"))

    m = hotlayer(θ, Tlo, Γ) .- hotlayer(θ, Thi, Γ)

    return clamp.(m, 0.0, 1.0)
end

"""
    hotlayer(θ, T, Γ)

Return the fraction of each wet tracer cell for which `θ ≥ T`.

The returned `MeshArray` has values between zero and one. It accounts for
linear variation of `θ` between adjacent wet cell centers, and therefore
allows fractional contributions from cells crossed by the threshold.
Non-monotonic profiles and multiple crossings are supported.

At the upper ocean boundary and at a wet/dry transition at the sea floor,
the relevant half-cell is classified using the adjacent wet cell-center
value. Dry cells, identified by `Γ.hFacC == 0`, are returned as zero.

The result is a fraction of the nominal tracer cell; `Γ.hFacC` is used
only to distinguish wet from dry cells. For a physical-volume integral,
multiply the result by `Γ.hFacC` and the cell volume.

See also: [`coldlayer`](@ref), [`layerfraction`](@ref), [`isosurface`](@ref).
"""
function hotlayer(θ, T, Γ)
    nr = size(θ, 2)
    m = zeros(θ)
    for j in 1:size(m, 1)

        # -- Upper half of the surface cell -------------------------------
        # This is a contribution from z = 0 to RC[1].
        # hFacC prevents a land point with an arbitrary θ fill value from contributing.
        wet1 = Γ.hFacC[j, 1] .> 0
        hot1 = θ[j, 1] .>= T
        m[j, 1] .+= 0.5 .* (wet1 .& hot1)

        # -- Interior interfaces ------------------------------------------
        for k in 1:nr-1
            θk  = θ[j, k]
            θk1 = θ[j, k+1]
            wetk  = Γ.hFacC[j, k]   .> 0
            wetk1 = Γ.hFacC[j, k+1] .> 0
            hotk  = θk  .>= T
            hotk1 = θk1 .>= T
            # Only interpolate across interfaces for which both cells are wet. 
            # A wet/dry interface is a sea-floor boundary, not a θ=T crossing.
            pairwet = wetk .& wetk1

            hk  = Γ.DRF[k]
            hk1 = Γ.DRF[k+1]
            zF  = Γ.RF[k+1]   # interface between levels k and k+1

            # Case 1: hot / hot
            i = findall(pairwet .& hotk .& hotk1)
            m[j, k][i]   .+= 0.5
            m[j, k+1][i] .+= 0.5

            # Case 2: hot / cold
            i = findall(pairwet .& hotk .& (.!hotk1))
            if !isempty(i)
                a = (θk[i] .- T) ./ (θk[i] .- θk1[i])
                zc = (1 .- a) .* Γ.RC[k] .+ a .* Γ.RC[k+1]
                # Lower half of cell k
                m[j, k][i] .+= clamp.(
                    (Γ.RC[k] .- zc) ./ hk,
                    0.0, 0.5,
                )
                # Upper half of cell k+1
                m[j, k+1][i] .+= clamp.(
                    (zF .- zc) ./ hk1,
                    0.0, 0.5,
                )
            end

            # Case 3: cold / hot
            i = findall(pairwet .& (.!hotk) .& hotk1)
            if !isempty(i)
                a = (θk[i] .- T) ./ (θk[i] .- θk1[i])
                zc = (1 .- a) .* Γ.RC[k] .+ a .* Γ.RC[k+1]
                # Lower half of cell k
                m[j, k][i] .+= clamp.(
                    (zc .- zF) ./ hk,
                    0.0, 0.5,
                )
                # Upper half of cell k+1
                m[j, k+1][i] .+= clamp.(
                    (zc .- Γ.RC[k+1]) ./ hk1,
                    0.0, 0.5,
                )
            end

            # -- Bottom half at a wet/dry transition ----------------------
            # If k is the deepest wet cell in a column, then its lower
            # half is bounded by the sea floor rather than another θ value.
            i = findall(wetk .& (.!wetk1) .& hotk)
            m[j, k][i] .+= 0.5
        end

        # -- Lower half of level nr ---------------------------------------
        # Needed only for columns that are wet through the final model
        # level. Shallower columns were handled by the wet/dry case above.
        wetn = Γ.hFacC[j, nr] .> 0
        hotn = θ[j, nr] .>= T
        m[j, nr] .+= 0.5 .* (wetn .& hotn)
    end

    return m
end

"""
    coldlayer(θ, T, Γ)

Return the fraction of each wet tracer cell for which `θ < T`.

This is the wet-cell complement of [`hotlayer`](@ref), which uses the
convention `θ ≥ T`. It therefore has identical vertical interpolation,
surface, sea-floor, and wet/dry conventions.

Dry cells, identified by `Γ.hFacC == 0`, are returned as zero. The result
is a fraction of the nominal tracer cell; include `Γ.hFacC` separately
when computing physical-volume integrals.

See also: [`hotlayer`](@ref), [`layerfraction`](@ref), [`isosurface`](@ref).
"""
function coldlayer(θ, T, Γ)
    m = 1.0 .* (Γ.hFacC .> 0)
    m .-= hotlayer(θ, T, Γ)
    return m
end
