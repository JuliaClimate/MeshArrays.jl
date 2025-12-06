
"""
    Tiles(γ::gcmgrid,ni::Int,nj::Int)

Define sudomain `tiles` of size `ni,nj`. Each tile is defined by a `Dict` where
`tile,face,i,j` correspond to tile ID, face ID, index ranges.

```jldoctest; output = false
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.Dataset("GRID_LLC90"))
τ=Tiles(γ,30,30)

isa(τ[1],NamedTuple)

# output

true
```
"""
function Tiles(γ::gcmgrid,ni::Int,nj::Int)
    nt=Int(prod(γ.ioSize)/ni/nj)
    τ=Array{Dict,1}(undef,nt)
    #
    cnt=0
    for iF=1:γ.nFaces
        for jj=Int.(1:γ.fSize[iF][2]/nj)
            for ii=Int.(1:γ.fSize[iF][1]/ni)
                cnt=cnt+1
                i=(1:ni).+ni*(ii-1)
                j=(1:nj).+nj*(jj-1)
                τ[cnt]=Dict("tile" => cnt, "face" => iF, "i" => i, "j" => j)
            end
        end
    end
    #
    return Dict_to_NamedTuple.(τ)
end

"""
    Tiles(τ::Array{Dict},x::MeshArray)

Return an `Array` of tiles which cover `x` according to tile partition `τ`.

```jldoctest; output = false
using MeshArrays
γ=GridSpec("LatLonCap",MeshArrays.Dataset("GRID_LLC90"))
d=γ.read(γ.path*"Depth.data",MeshArray(γ,γ.ioPrec))
τ=Tiles(γ,30,30)
td=Tiles(τ,d)

D=similar(d)
Tiles!(τ,td,D)

isa(td[1],Array)

# output

true
```
"""
function Tiles(τ::Array,x::MeshArray)
    nt=length(τ)
    tx=Array{typeof(x[1]),1}(undef,nt)
    dn=size(x[1],1)-x.fSize[1][1]
    for ii=1:nt
        f=τ[ii].face
        i0=minimum(τ[ii].i)
        i1=maximum(τ[ii].i)
        j0=minimum(τ[ii].j)
        j1=maximum(τ[ii].j)
        tx[ii]=view(x[f],i0:i1+dn,j0:j1+dn)
    end
    return tx
end

"""
    Tiles!(τ::Array,tx::Array,x::MeshArrays)

Map tiles in `tx` according to tile partition `τ` into `x`.
"""
function Tiles!(T::Array,tx::Array,x::MeshArray)
    nt=length(T)
    dn=size(x[1],1)-x.fSize[1][1]
    for ii=1:nt
        f=T[ii].face
        i0=minimum(T[ii].i)
        i1=maximum(T[ii].i)
        j0=minimum(T[ii].j)
        j1=maximum(T[ii].j)
        x[f][i0:i1+dn,j0:j1+dn].=tx[ii]
    end
end
