


module NEMO_GRID

import MeshArrays

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
	(:gphif, :gphif, :Zeta),
	(:glamt, :XC, :T),
	(:glamu, :XW, :U),
	(:glamv, :XS, :V),
	(:glamf, :glamf, :Zeta),
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
				ii=2:1441, jj=2:1021, is_3d=false, verbose=false)

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
	
	if is_3d
		grid_data[nam_in][ii.+di,jj.+dj,:]
	else
		grid_data[nam_in][ii.+di,jj.+dj,1]
	end
end

"""
    read_nc_grid(grid_data; verbose=false)
    
Read in the grid variables from `NEMO`, so they can 
then easily be converted to `MeshArrays`.

```
using NCDatasets
data=Dataset("mesh_mask_ORCA025.nc")

import MeshArrays
grid=MeshArrays.NEMO_GRID.read_nc_grid(data)
Γ=MeshArrays.NEMO_GRID.grid_to_MeshArrays(grid)

using CairoMakie
heatmap(Γ.nanmask*Γ.RAC)
```
"""
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

    g=MeshArrays.GridSpec("PeriodicChannel")
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

end
