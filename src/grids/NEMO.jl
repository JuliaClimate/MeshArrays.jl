module NEMO_GRID

import MeshArrays: GridSpec, MeshArray_wh, exchange
import MeshArrays: gcmgrid, varmeta, defaultmeta
import MeshArrays: InnerArray, OuterArray, AbstractMeshArray, thisversion, gcmarray

struct nemoarray{T, N, AT} <: AbstractMeshArray{T, N}
   grid::gcmgrid
   meta::varmeta
   f::OuterArray{AT,N}
   fSize::OuterArray{NTuple{2, Int}}
   fIndex::OuterArray{Int,1}
   version::String
end

nemoarray(stuff...;kwargs...)=gcmarray(stuff...;kwargs...)

one_nemoarray(g::gcmgrid) = begin
	T=g.ioPrec
	f=[zeros(g.fSize[1])]
	fIndex=[1]
	nemoarray{T,1,InnerArray{T,2}}(g,defaultmeta,f,g.fSize,[1],thisversion)
end

gcmarray_to_nemorarray(A::gcmarray) = 
	nemoarray{eltype(A),ndims(A),InnerArray{eltype(A),2}}(A.grid,defaultmeta,copy(A.f),copy(A.fSize),copy(A.fIndex),thisversion)
nemorarray_to_gcmarray(A::nemoarray) = 
	gcmarray{eltype(A),ndims(A),InnerArray{eltype(A),2}}(A.grid,defaultmeta,copy(A.f),copy(A.fSize),copy(A.fIndex),thisversion)

#this approach works in 3D	
function Base.similar(A::nemoarray;m::varmeta=defaultmeta)
    if ndims(A)==1
        B=gcmarray(similar(A.grid),eltype(A),copy(A.fSize),copy(A.fIndex); meta=m)
    else
        B=gcmarray(similar(A.grid),eltype(A),copy(A.fSize),copy(A.fIndex),size(A)[2:end]...; meta=m)
    end
    return gcmarray_to_nemorarray(B)
end

#Base.getindex(A::nemoarray,args...;kwargs...)=
#	Base.getindex(nemorarray_to_gcmarray(A),args...;kwargs...)

Base.dataids(A::nemoarray) = (Base.dataids(A.f)..., Base.dataids(A.fSize)..., Base.dataids(A.fIndex)...)

##

class="PeriodicChannel"
ioSize=(1440, 1020)
ioPrec=Float64

GridSpec_NEMO(path) = gcmgrid(path,class,1,[ioSize],ioSize,ioPrec,read,write)

##

variable_list_2d=(
	(:hdept, :Depth, :T),
	(:e1t, :DXW, :T),
	(:e1u, :DXC, :U),
	(:e1v, :DXG, :V),
	(:e1f, :DXS, :Zeta),
	(:e2t, :DYS, :T),
	(:e2u, :DYG, :U),
	(:e2v, :DYC, :V),
	(:e2f, :DYW, :Zeta),
	(:gphit, :YC, :T),
	(:gphiu, :YW, :U),
	(:gphiv, :YS, :V),
	(:gphif, :YG, :Zeta),
	(:glamt, :XC, :T),
	(:glamu, :XW, :U),
	(:glamv, :XS, :V),
	(:glamf, :XG, :Zeta),
)
#derived variables	
variable_list_derived=(
	(:e1t_X_e2t, :RAC, :T),
	(:e1u_X_e2u, :RAW, :U),
	(:e1v_X_e2v, :RAS, :V),
	(:e1f_X_e2f, :RAZ, :Zeta),
)
#3D arrays:	
variable_list_3d=(
	(:tmask, :hFacC, :T),
	(:umask, :hFacW, :U),
	(:vmask, :hFacS, :V),
)
#1D profiles:
variable_list_1d=(
	(:gdept_0,:RC, :RC),	
	(:gdepw_0,:RF, :RF),	
	(:e3t_0,:DRF, :RC),	
	(:e3w_0,:DRC, :RF),	
)

variable_NTA(variable_list=variable_list_2d)=
	[(NEMO=x[1],MITgcm=x[2],location=x[3]) for x in variable_list]

variable_in_NEMO(v,vl=variable_list_2d)=
	vl[findall([i[2]==v for i in vl])[1]][1]
#DataFrame version:	
#	filter(p->p.:MITgcm==v,variable_df(variable_list))[1,:NEMO]

function convert_one_grid_variable(grid_data,df_line; 
				is_3d=false, verbose=false)	
	nam_in=df_line.NEMO
	nam_out=df_line.MITgcm
	loc=df_line.location
	verbose ? println([nam_in nam_out loc]) : nothing
	read_one(grid_data[nam_in],loc,is_3d)
end

nomissing64(x) = [(ismissing(a) ? 0 : Float64(a)) for a in x]
nomissing64(x::UnitRange) = Float64.(x)
nonan64(x) = [(isnan(a) ? 0 : Float64(a)) for a in x]

function read_one(fld,loc,is_3d)
	#for T/U points : jj=1:1020
	#for V/Zeta points : add a line of zeros at South on V + 1:1019
	ii=2:1441
	jj=1:1020

	list_loc=[:T :Zeta :U :V :RC :RF]
	list_di =[ 0 -1    -1  0  0   0  ]
	list_dj =[ 0 -1     0 -1  0   0  ]
	
	k=findall(list_loc.==loc)[1]
	di=list_di[k]
	dj=list_dj[k]

	jj_m_1=jj[2:end].-1
	
	tmp=if is_3d
		if dj<0
			tmp=fld[ii.+di,jj_m_1,:]
			nr=size(tmp,3)
			cat(zeros(length(ii),1,nr),tmp, dims=2)
		else
			fld[ii.+di,jj.+dj,:]
		end
	else
		if dj<0
			tmp=fld[ii.+di,jj_m_1,1]
			cat(zeros(length(ii),1),tmp, dims=2)
		else
			fld[ii.+di,jj.+dj,1]
		end
	end
	nomissing64(tmp)
end

"""
    GridLoad_NEMO(grid_data=Any[]; verbose=false)
    
Read in the grid variables from `NEMO` grid file 
and convert them into a `MeshArray` grid.

```
using MeshArrays
γ=MeshArrays.GridSpec_NEMO(joinpath("data","mesh_mask_ORCA025.nc"))
using NCDatasets; grid_data=Dataset(γ.path)
Γ=MeshArrays.GridLoad_NEMO(grid_data)

using CairoMakie
heatmap(Γ.nanmask*Γ.RAC)
```
"""
function GridLoad_NEMO(grid_data=Any[]; verbose=false)
	G=read_nc_grid(grid_data,verbose=verbose)
	G=grid_to_MeshArrays(G)
	add_angle_CS_SN(G)
end

load(grid_data; verbose=false)=begin
	G=read_nc_grid(grid_data,verbose=verbose)
	G=grid_to_MeshArrays(G)
	add_angle_CS_SN(G)
end

function read_nc_grid(all_grid_data; verbose=false)
	if isa(all_grid_data,Tuple)
		grid_data=all_grid_data[1]
		ds_e3=all_grid_data[2]
	else
		grid_data=all_grid_data
		ds_e3=missing
	end

	grid=Dict()

	for i in variable_NTA(variable_list_2d)
		tmp=convert_one_grid_variable(grid_data,i,verbose=verbose)
		merge!(grid,Dict(i.MITgcm=>tmp))
	end

	merge!(grid,Dict(:RAC=>grid[:DXW].*grid[:DYS]))
	merge!(grid,Dict(:RAW=>grid[:DXC].*grid[:DYG]))
	merge!(grid,Dict(:RAS=>grid[:DXG].*grid[:DYC]))
	merge!(grid,Dict(:RAZ=>grid[:DXS].*grid[:DYW]))

	for i in variable_NTA(variable_list_3d)
		tmp=convert_one_grid_variable(grid_data,i,is_3d=true,verbose=verbose)
		merge!(grid,Dict(i.MITgcm=>tmp))
	end

	add_one_dim_variables!(grid,grid_data)

	ismissing(ds_e3) ? nothing : overwrite_hFac!(grid,ds_e3)

	mask=Float64.(1*grid[:hFacC][:,:,1])
	mask[findall(mask.==0)].=NaN
	mask[findall(mask.>0)].=1.0
	merge!(grid,Dict("nanmask"=>mask))

	verbose ? println(keys(grid)) : nothing
	
	NamedTuple((Symbol(key),value) for (key,value) in grid)
end

function grid_to_MeshArrays(grid)
	ni,nj=size(grid.XC)

    g=GridSpec("PeriodicChannel",ioPrec=Float64)
    g.fSize[1]=(ni,nj)
    g.ioSize.=[ni nj]

	grid2=Dict()
	for i in keys(grid)
		if ndims(grid[i])>1
			tmp=read(grid[i],g)
			tmp=gcmarray_to_nemorarray(tmp)
			merge!(grid2,Dict(i=>tmp))
		else
			merge!(grid2,Dict(i=>grid[i]))
		end
	end

	NamedTuple((Symbol(key),value) for (key,value) in grid2)
end

function add_one_dim_variables!(grid,grid_data; verbose=false)
	for df_line in variable_NTA(variable_list_1d)
		nam_in=df_line.NEMO
		nam_out=df_line.MITgcm
		loc=df_line.location

		fac=(in(nam_out,[:RC,:RF]) ? -1 : 1)
		verbose ? println(typeof(grid_data[nam_in])) : nothing
		tmp=fac*nomissing64.(grid_data[nam_in][:])
		merge!(grid,Dict(nam_out=>tmp))
	end
end

## 3D grid factors needed for transport calculations

function overwrite_hFac!(Γ,ds; verbose=false)
	verbose ? println("overwriting hFac") : nothing
	DRF=Γ[:DRF]
	e3t=read_one(ds["e3t_0_field"],:T,true)
	e3u=read_one(ds["e3u_0_field"],:U,true)
	e3v=read_one(ds["e3v_0_field"],:V,true)
#why is e3t_0_field ~ 2xDRF at last point? 
#and this seems necessary to compute correct top-bottom transports
#	e3t[:,:,end-1].=e3t[:,:,end-1]./2
#	e3u[:,:,end-1].=e3u[:,:,end-1]./2
#	e3v[:,:,end-1].=e3v[:,:,end-1]./2
	[Γ[:hFacC][:,:,k].=e3t[:,:,k]./Float64(DRF[k]) for k in 1:75]
	[Γ[:hFacW][:,:,k].=e3u[:,:,k]./Float64(DRF[k]) for k in 1:75]
	[Γ[:hFacS][:,:,k].=e3v[:,:,k]./Float64(DRF[k]) for k in 1:75]
end

##

#exchange(x::nemoarray) = NEMO_exchange(x)

exchange(x::nemoarray)=begin
	y=MeshArray_wh(similar(x),1)
	if length(size(x))==1
		y.MA.f[1]=NEMO_exchange(x).MA.f[1]
	else
		for k in 1:size(x)[2]
			y.MA.f[k]=NEMO_exchange(x[:,k]).MA.f[1]
		end
	end
	y
end

"""
    NEMO_exchange(x; fac=1.0, U_shift=false)

Add Halos with neighboring values, based on NEMO's folds.
"""
function NEMO_exchange(x; fac=1.0, U_shift=false)
	y=zeros(1442,1022)
	y[2:end-1,2:end-1].=x[1]
	y[1,:].=y[end-1,:]
	y[end,:].=y[2,:]
	if !U_shift
		y[2:end-1,end].=fac*reverse(y[3:end,end-2])
		y[1,end]=y[end-1,end]
		y[end,end]=y[2,end]
	else
		tmp1=circshift(y[3:end,end-2],-1)
		y[2:end-1,end].=fac*reverse(tmp1)
		y[1,end]=y[end-1,end]
		y[end,end]=y[2,end]
	end
	y

	yy=similar(x)
	yy[1]=y
	MeshArray_wh(yy,1)
end

function add_angle_CS_SN(Γ)
	r_earth=sqrt(sum(Γ.RAC))./(4 .*pi)
	ni,nj=Γ.RAC.fSize[1]
	grid=Γ.RAC.grid

	psi=NEMO_exchange(-deg2rad(1)*r_earth*Γ.YC).MA[1]
	uZ=psi[1:ni,1:nj]-psi[1:ni,2:nj+1]
	vZ=psi[2:ni+1,1:nj]-psi[1:ni,1:nj]

	DXG=NEMO_exchange(Γ.DXG).MA[1]
	DYG=NEMO_exchange(Γ.DYG).MA[1]
	uZ=uZ./DYG[1:ni,1:nj]#shift by 1 point?
	vZ=vZ./DXG[1:ni,1:nj]

	uZ=read(uZ,grid)
	vZ=read(vZ,grid)
	uu=NEMO_exchange(uZ,fac=-1).MA[1]
	vv=NEMO_exchange(vZ,fac=-1).MA[1]

	uZc=(uu[1:ni,1:nj]+uu[2:ni+1,1:nj])/2
	vZc=(vv[1:ni,1:nj]+vv[1:ni,2:nj+1])/2

	norm=sqrt.(uZc.*uZc+vZc.*vZc)
	AngleCS =  read(uZc./norm,grid)
	AngleSN =  read(-vZc./norm,grid)

	merge(Γ,(AngleCS=AngleCS,AngleSN=AngleSN))
end

end
