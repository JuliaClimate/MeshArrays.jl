module NEMO_GRID

import MeshArrays: GridSpec, MeshArray_wh

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

	ii=2:1441
	jj=1:1020

	#note on jj:
	#for V/Zeta points : add a line of zeros at South on V + 1:1019
	#for T/U points : 1:1020
	
	nam_in=df_line.NEMO
	nam_out=df_line.MITgcm
	loc=df_line.location
	verbose ? println([nam_in nam_out loc]) : nothing

	list_loc=[:T :Zeta :U :V :RC :RF]
	list_di =[ 0 -1    -1  0  0   0  ]
	list_dj =[ 0 -1     0 -1  0   0  ]
	
	k=findall(list_loc.==loc)[1]
	di=list_di[k]
	dj=list_dj[k]

	jj_m_1=jj[2:end].-1
	
	if is_3d
		if dj<0
			tmp=grid_data[nam_in][ii.+di,jj_m_1,:]
			nr=size(tmp,3)
			cat(zeros(length(ii),1,nr),tmp, dims=2)
		else
			grid_data[nam_in][ii.+di,jj.+dj,:]
		end
	else
		if dj<0
			tmp=grid_data[nam_in][ii.+di,jj_m_1,1]
			cat(zeros(length(ii),1),tmp, dims=2)
		else
			grid_data[nam_in][ii.+di,jj.+dj,1]
		end
	end
end

"""
    load(grid_data; verbose=false)
    
Read in the grid variables from `NEMO` grid file 
and convert them into a `MeshArray` grid.

```
using NCDatasets
grid_data=Dataset("mesh_mask_ORCA025.nc")

using MeshArrays
Γ=NEMO_GRID.load(grid_data)

using CairoMakie,
heatmap(Γ.nanmask*Γ.RAC)

LC=LatitudeCircles(-89.0:89.0,Γ)

using JLD2
lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5];
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5];
λ=MeshArrays.interpolation_setup(Γ=Γ,lon=lon,lat=lat);

lo,la,z=MeshArrays.Interpolate(Γ.Depth,λ)
heatmap(z)

path="data/NEMO_sample/monthly/"
function read_vel(path,month=1)
	m=string(month)
	uT=Dataset(joinpath(path,"ORAS5_uT_1993-"*m*".nc"))["uT"][2:end-1,1:1020,1]
	uT=read(uT,g)
	vT=Dataset(joinpath(path,"ORAS5_vT_1993-"*m*".nc"))["vT"][2:end-1,1:1020,1]
	vT=read(vT,g)
	(uT,vT)
end

G=NEMO_GRID.calc_angle(Γ)
(uT,vT)=read_vel(path,1)
uT_E,vT_N=UVtoUEVN(uT,vT,G)
_,_,z_E=MeshArrays.Interpolate(uT_E,λ);
_,_,z_N=MeshArrays.Interpolate(vT_N,λ);
heatmap(z_E,colorrange=(-1,1).*1e3)
heatmap(z_N,colorrange=(-1,1).*1e3)

MT=zeros(179,12)
for m in 1:12
	uT,vT=read_vel(path,m)
	UV=Dict("U"=>Γ.DYG*uT,"V"=>Γ.DXG*vT,"dimensions"=>["x","y"])
	MT=1e-15*4e6*[ThroughFlow(UV,lc,Γ) for lc in LC]
end
```
"""
load(grid_data; verbose=false)=
	grid_to_MeshArrays(read_nc_grid(grid_data,verbose=verbose))

function read_nc_grid(grid_data; verbose=false)
	grid=Dict()

	for i in variable_NTA(variable_list_2d)
		tmp=convert_one_grid_variable(grid_data,i,verbose=verbose)
		merge!(grid,Dict(i.MITgcm=>tmp))
	end

	merge!(grid,Dict("RAC"=>grid[:DXW].*grid[:DYS]))
	merge!(grid,Dict("RAW"=>grid[:DXC].*grid[:DYG]))
	merge!(grid,Dict("RAS"=>grid[:DXG].*grid[:DYC]))
	merge!(grid,Dict("RAZ"=>grid[:DXS].*grid[:DYW]))

	for i in variable_NTA(variable_list_3d)
		tmp=convert_one_grid_variable(grid_data,i,is_3d=true,verbose=verbose)
		merge!(grid,Dict(i.MITgcm=>tmp))
	end

	add_one_dim_variables!(grid,grid_data)

	mask=Float64.(1*grid[:hFacC][:,:,1])
	mask[findall(mask.==0)].=NaN
	mask[findall(mask.>0)].=1.0
	merge!(grid,Dict("nanmask"=>mask))

	verbose ? println(keys(grid)) : nothing
	
	NamedTuple((Symbol(key),value) for (key,value) in grid)
end

function grid_to_MeshArrays(grid)
	ni,nj=size(grid.XC)

    g=GridSpec("PeriodicChannel")
    g.fSize[1]=(ni,nj)
    g.ioSize.=[ni nj]

	grid2=Dict()
	for i in keys(grid)
		if ndims(grid[i])>1
			tmp=read(grid[i],g)
			merge!(grid2,Dict(i=>tmp))
		else
			merge!(grid2,Dict(i=>grid[i]))
		end
	end

	NamedTuple((Symbol(key),value) for (key,value) in grid2)
end

function add_one_dim_variables!(grid,grid_data)
	for df_line in variable_NTA(variable_list_1d)
		nam_in=df_line.NEMO
		nam_out=df_line.MITgcm
		loc=df_line.location

		fac=(in(nam_out,[:RC,:RF]) ? -1 : 1)
		tmp=fac*grid_data[nam_in][:]
		merge!(grid,Dict(nam_out=>tmp))
	end
end

##

"""
    exchange(x; V_fac=1.0, U_shift=false)

Add Halos with neighboring values, based on NEMO's folds.
"""
function exchange(x; V_fac=1.0, U_shift=false)
	y=zeros(1442,1022)
	y[2:end-1,2:end-1].=x[1]
	y[1,:].=y[end-1,:]
	y[end,:].=y[2,:]
	if !U_shift
		y[2:end-1,end].=V_fac*reverse(y[3:end,end-2])
		y[1,end]=y[end-1,end]
		y[end,end]=y[2,end]
	else
		tmp1=circshift(y[3:end,end-2],-1)
		y[2:end-1,end].=V_fac*reverse(tmp1)
		y[1,end]=y[end-1,end]
		y[end,end]=y[2,end]
	end
	y

	yy=similar(x;m=x.meta);
	yy[1]=y
	MeshArray_wh(yy,1)
end

function calc_angle(Γ)
	r_earth=sqrt(sum(Γ.RAC))./(4 .*pi)
	ni,nj=Γ.RAC.fSize[1]
	grid=Γ.RAC.grid

	psi=exchange(-deg2rad(1)*r_earth*Γ.YC).MA[1]
	uZ=psi[1:ni,1:nj]-psi[1:ni,2:nj+1]
	vZ=psi[2:ni+1,1:nj]-psi[1:ni,1:nj]

	DXG=exchange(Γ.DXG).MA[1]
	DYG=exchange(Γ.DYG).MA[1]
	uZ=uZ./DYG[1:ni,1:nj]#shift by 1 point?
	vZ=vZ./DXG[1:ni,1:nj]

	uZ=read(uZ,grid)
	vZ=read(vZ,grid)
	uu=exchange(uZ).MA[1]
	vv=exchange(vZ).MA[1]

	uZc=(uu[1:ni,1:nj]+uu[2:ni+1,1:nj])/2
	vZc=(vv[1:ni,1:nj]+vv[1:ni,2:nj+1])/2

	norm=sqrt.(uZc.*uZc+vZc.*vZc)
	AngleCS =  read(uZc./norm,grid)
	AngleSN =  read(-vZc./norm,grid)

	merge(Γ,(AngleCS=AngleCS,AngleSN=AngleSN))
end

end
