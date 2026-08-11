
## This file contains the exchange and exch_UV functions
# along with grid-specific methods (exch_T_N.jl, etc.)

## User Front Ends
#
# exchange!  — in-place variants; caller pre-allocates output MeshArray_wh
#
# Allocate a suitable output buffer with:
#   buf = exchange_alloc(x)          # for scalar field x
#   buf_u, buf_v = exchange_alloc(u, v)  # for vector field u, v

"""
    exchange_alloc(fld::AbstractMeshArray, N=1)

Allocate a `MeshArray_wh` output buffer compatible with `exchange!(buf, fld)`.
"""
function exchange_alloc(fld::AbstractMeshArray, N::Integer=1)
    nf = fld.grid.nFaces
    # fld.fSize is always a 1D array of (nx,ny) tuples, one per face
    fs = fld.fSize
    nf == 5 ? fs = vcat(fs, fs[3]) : nothing
    out = MeshArray_wh(similar(fld; m=fld.meta), N)
    if length(size(fld)) == 1
        for i in 1:nf; out.MA.f[i] = fill(0.0, fs[i] .+ 2N); end
    else
        for k in 1:size(fld)[2], i in 1:nf
            out.MA.f[i, k] = fill(0.0, fs[i] .+ 2N)
        end
    end
    out
end

"""
    exchange_alloc(u::AbstractMeshArray, v::AbstractMeshArray, N=1)

Allocate two `MeshArray_wh` output buffers compatible with `exchange!(bu, bv, u, v)`.
"""
function exchange_alloc(u::AbstractMeshArray, v::AbstractMeshArray, N::Integer=1)
    exchange_alloc(u, N), exchange_alloc(v, N)
end

"""
    exchange!(y::MeshArray_wh, x::AbstractMeshArray)

In-place exchange: fill pre-allocated `y` from `x` without allocating new arrays.
`y` must have been created with `exchange_alloc(x)` (or equivalent).
"""
function exchange!(y::MeshArray_wh, x::AbstractMeshArray)
    N = y.HS
    if length(size(x)) == 1
        exchange_main!(y, x, N)
    else
        for k in 1:size(x)[2]
            exchange_main!(MeshArray_wh(y.MA[:, k], N), x[:, k], N)
        end
    end
    y
end

"""
    exchange!(yu::MeshArray_wh, yv::MeshArray_wh, u::AbstractMeshArray, v::AbstractMeshArray)

In-place exchange for a vector field pair.
`yu`, `yv` must have been created with `exchange_alloc(u, v)`.
"""
function exchange!(yu::MeshArray_wh, yv::MeshArray_wh,
                   u::AbstractMeshArray, v::AbstractMeshArray)
    N = yu.HS
    if length(size(u)) == 1
        exchange_main!(yu, yv, u, v, N)
    else
        for k in 1:size(u)[2]
            exchange_main!(MeshArray_wh(yu.MA[:, k], N),
                           MeshArray_wh(yv.MA[:, k], N),
                           u[:, k], v[:, k], N)
        end
    end
    yu, yv
end

## in-place exchange_main! — dispatch over grid types

function exchange_main!(y::MeshArray_wh, x::AbstractMeshArray, N::Integer)
    if x.grid.class == "LatLonCap" || x.grid.class == "CubeSphere"
        exch_T_N_cs!(y, x, N)
    elseif x.grid.class == "PeriodicChannel"
        exch_T_N_PeriodicChannel!(y, x, N)
    elseif x.grid.class == "PeriodicDomain"
        exch_T_N_PeriodicDomain!(y, x, N)
    else
        error("unknown grid.class case")
    end
    y
end

function exchange_main!(yu::MeshArray_wh, yv::MeshArray_wh,
                        u::AbstractMeshArray, v::AbstractMeshArray, N::Integer)
    if u.grid.class == "LatLonCap" || u.grid.class == "CubeSphere"
        exch_UV_N_cs!(yu, yv, u, v, N)
    elseif u.grid.class == "PeriodicChannel"
        exch_UV_N_PeriodicChannel!(yu, yv, u, v, N)
    elseif u.grid.class == "PeriodicDomain"
        exch_UV_N_PeriodicDomain!(yu, yv, u, v, N)
    else
        error("unknown grid.class case")
    end
    yu, yv
end

##

"""
    exchange(fld::AbstractMeshArray)

Exchange / transfer data between neighboring arrays. Other methods are

    exchange(fld::AbstractMeshArray,N::Integer)
    exchange(u::AbstractMeshArray,v::AbstractMeshArray)
    exchange(u::AbstractMeshArray,v::AbstractMeshArray,N::Integer)
"""
function exchange(x::AbstractMeshArray)
	y=MeshArray_wh(similar(x),1)
	if length(size(x))==1
    tmp=exchange_main(x).MA
		[y.MA.f[k]=tmp.f[k] for k in 1:length(x)] 
	else
		for k in 1:size(x)[2]
      tmp=exchange_main(x[:,k]).MA
			[y.MA.f[kk,k]=tmp.f[kk] for kk in 1:size(x,1)] 
		end
	end
	y
end

function exchange(u::AbstractMeshArray,v::AbstractMeshArray)
	uu=MeshArray_wh(similar(u),1)
  vv=MeshArray_wh(similar(v),1)
	if length(size(v))==1
    tmp_u,tmp_v=exchange_main(u,v)
		[uu.MA.f[k]=tmp_u.MA.f[k] for k in 1:length(u)] 
		[vv.MA.f[k]=tmp_v.MA.f[k] for k in 1:length(v)] 
	else
		for k in 1:size(v)[2]
      tmp_u,tmp_v=exchange_main(u[:,k],v[:,k])
			[uu.MA.f[kk,k]=tmp_u.MA.f[kk] for kk in 1:size(u,1)] 
			[vv.MA.f[kk,k]=tmp_v.MA.f[kk] for kk in 1:size(v,1)] 
		end
	end
	uu,vv
end

function exchange_main(fld::AbstractMeshArray)
  MeshArray_wh(exch_T_N(fld,1),1)
end

function exchange_main(fld::AbstractMeshArray,N::Integer)
  MeshArray_wh(exch_T_N(fld,N),N)
end

function exchange_main(u::AbstractMeshArray,v::AbstractMeshArray)
  MeshArray_wh.(exch_UV_N(u,v,1),1)
end

function exchange_main(u::AbstractMeshArray,v::AbstractMeshArray,N::Integer)
  MeshArray_wh.(exch_UV_N(u,v,N),N)
end

## dispatch over grid types

#note: the "CubeSphere" implementation covers both cs and llc

function exch_T_N(fld,N)

if fld.grid.class=="LatLonCap"
  FLD=exch_T_N_cs(fld,N)
elseif fld.grid.class=="CubeSphere"
  FLD=exch_T_N_cs(fld,N)
elseif fld.grid.class=="PeriodicChannel"
  FLD=exch_T_N_PeriodicChannel(fld,N)
elseif fld.grid.class=="PeriodicDomain"
  FLD=exch_T_N_PeriodicDomain(fld,N)
else
  error("unknown grid.class case")
end

return FLD

end

function exch_UV_N(u,v,N)

if u.grid.class=="LatLonCap"
  (uex,vex)=exch_UV_N_cs(u,v,N)
elseif u.grid.class=="CubeSphere"
  (uex,vex)=exch_UV_N_cs(u,v,N)
elseif u.grid.class=="PeriodicChannel"
  (uex,vex)=exch_UV_N_PeriodicChannel(u,v,N)
elseif u.grid.class=="PeriodicDomain"
  (uex,vex)=exch_UV_N_PeriodicDomain(u,v,N)
else
  error("unknown grid.class case")
end

return uex,vex

end

function exch_UV(u,v)

if u.grid.class=="LatLonCap"
  (uex,vex)=exch_UV_cs(u,v)
elseif u.grid.class=="CubeSphere"
  (uex,vex)=exch_UV_cs(u,v)
elseif u.grid.class=="PeriodicChannel"
  (uex,vex)=exch_UV_PeriodicChannel(u,v)
elseif u.grid.class=="PeriodicDomain"
  (uex,vex)=exch_UV_PeriodicDomain(u,v)
else
  error("unknown grid.class case")
end

return uex,vex

end

## Grid-specific implementations:

include("PeriodicDomain.jl")
include("PeriodicChannel.jl")
include("CubeSphere.jl")
