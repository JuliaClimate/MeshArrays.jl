
module Polygons

import GeoInterface as GI
using Glob

## functions to derive polygons from XG,YG read from native grid 

function to_sphere(XC,YC)
	x=sin.(pi/2 .-YC*pi/180).*cos.(XC*pi/180);
	y=sin.(pi/2 .-YC*pi/180).*sin.(XC*pi/180);
	z=cos.(pi/2 .-YC*pi/180);
	(x,y,z)
end

function treat_180lon!(x::Matrix{Float64}; level=2)
	if (maximum(x)-minimum(x))>180
		x[x.>0].=x[x.>0].-360
        if level>1
            x[x.>180].=180
            x[x.<-180].=-180
        end
	end
end

function to_LineStrings3D(pol)
	arr2=Array{Any}(undef,32,32)
	b0=Array{Any}(undef,5)
	for ij in eachindex(pol)
		geom=GI.getgeom(pol[ij])
		a=GI.getpoint.(geom)
		b=GI.coordinates.(geom)
		b1=GI.coordinates.(a[1])
		for p in 1:5
#			b0[p,:].=to_sphere(GI.coordinates(a[1][p])...)
			b0[p]=GI.Point(to_sphere(GI.coordinates(a[1][p])...))
		end
#		arr2[ij]=deepcopy(b0)
		arr2[ij]=deepcopy(GI.LineString(b0))
	end
	arr2
end

function to_Polygons(XG,YG,ff=1)
	arr2=Array{Any}(undef,32,32)
	IJ=1:32
	for i in IJ
		for j in IJ
			x=[XG[ff][i,j] XG[ff][i+1,j] XG[ff][i+1,j+1] XG[ff][i,j+1] XG[ff][i,j]]
			treat_180lon!(x)
			y=[YG[ff][i,j] YG[ff][i+1,j] YG[ff][i+1,j+1] YG[ff][i,j+1] YG[ff][i,j]]  
			arr2[i,j]=GI.Polygon([ GI.LinearRing(GI.Point.(zip(vec(x),vec(y)))) ])
#			arr2[i,j]=GI.LineString(GI.Point.(zip(vec(x),vec(y))))
		end
	end
	arr2
end

function to_LineStrings2D(XG,YG,ff=1;do_sphere=true)
	arr2=Array{Any}(undef,32,32)
	IJ=1:32
	for i in IJ
		for j in IJ
			x=[XG[ff][i,j] XG[ff][i+1,j] XG[ff][i+1,j+1] XG[ff][i,j+1] XG[ff][i,j]] 
			treat_180lon!(x)
			y=[YG[ff][i,j] YG[ff][i+1,j] YG[ff][i+1,j+1] YG[ff][i,j+1] YG[ff][i,j]]  
			arr2[i,j]=GI.LineString(GI.Point.(zip(vec(x),vec(y))))
		end
	end
	arr2
end

## demo_grid to read in cs32 grid as test case

"""
    demo_grid(path_grid)


Example that associates each grid cell with a polygon.
Then plot all polygons on a map using `Depth` as color.

```
using MeshArrays, MITgcm, CairoMakie

path_MITgcm=MITgcm.getdata("mitgcmsmallverif")
path_grid=joinpath(path_MITgcm,"MITgcm",
    "verification","tutorial_held_suarez_cs","input")

pols,pols3D=MeshArrays.Polygons.demo_grid(path_grid)
Depth=GridLoadVar("Depth",GridSpec(ID=:CS32))
MeshArrays.plot_examples(:polygons_plot,pols,color=Depth)
```
"""
function demo_grid(path_grid)
    files_grid=glob("grid_cs32.face00?.bin",path_grid)
    list_fields=["XC","YC","DXF","DYF","RAC","XG","YG","DXV","DYU","RAZ",
        "DXC","DYC","RAW","RAS","DXG","DYG","AngleCS","AngleSN"]
    grid=native_grid(files_grid,list_fields)
    (; XG, YG) = grid

    pols=[to_Polygons(XG,YG,ff) for ff in 1:length(grid.XG)]
	pols3D=[to_LineStrings3D(p) for p in pols]

	return pols,pols3D
end

## helper functions for demo_grid

function read_native_grid(files_grid,list_fields,ff=1,va="XG")
    fil=files_grid[ff]
    fid = open(fil)
    (n1,n2,n3)=(33,33,18)
    xx = Array{Float64,1}(undef,n1*n2*n3)
    read!(fid,xx)
    xx = reshape(hton.(xx),(n1,n2,n3))
    close(fid)
    vv=findall(list_fields.==va)[1]
    xx[:,:,vv]
end

function native_grid(files_grid,list_fields)
    XG=Array{Any}(undef,6)
    YG=Array{Any}(undef,6)
    for ff=1:6
        XG[ff]=read_native_grid(files_grid,list_fields,ff,"XG")
        YG[ff]=read_native_grid(files_grid,list_fields,ff,"YG")
    end
    (XG=XG,YG=YG)
end

end
