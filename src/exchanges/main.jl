
## This file contains the exchange and exch_UV functions
# along with grid-specific methods (exch_T_N.jl, etc.)

## User Front Ends

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
    tmp=exchange_main(x[:,k]).MA
		for k in 1:size(x)[2]
			[y.MA.f[kk,k]=tmp.f[kk] for kk in 1:size(x,1)] 
		end
	end
	y
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
