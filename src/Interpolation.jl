
using NearestNeighbors
import NearestNeighbors: knn

"""
    knn(xgrid,ygrid::MeshArray,x,y::Array{T,1},k::Int)

Find k nearest neighbors to each point in x,y on xgrid,ygrid

```
lon=collect(0.1:0.5:2.1); lat=collect(0.1:0.5:2.1);
(f,i,j,c)=knn(Γ.XC,Γ.YC,lon,lat)
```
"""
function knn(xgrid::MeshArray,ygrid::MeshArray,
        xvec::Array{T,1},yvec::Array{T,1},k=1::Int) where {T}

        #ancillary variables
        γ=xgrid.grid
        a_f=MeshArray(γ,Int); [a_f[ii][:,:].=ii for ii=1:γ.nFaces]
        a_i=MeshArray(γ,Int); [a_i[ii]=collect(1:γ.fSize[ii][1])*ones(Int,1,γ.fSize[ii][2]) for ii=1:γ.nFaces]
        a_j=MeshArray(γ,Int); [a_j[ii]=ones(Int,γ.fSize[ii][1],1)*collect(1:γ.fSize[ii][2])' for ii=1:γ.nFaces]

        #convert to flat Array format
        a_x=write(xgrid)
        a_y=write(ygrid)
        a_f=write(a_f)
        a_i=write(a_i)
        a_j=write(a_j)

        #vector of grid points in Cartesian, 3D, coordinates
        kk=findall(isfinite.(a_x))
        x=sin.(pi/2 .-a_y[kk]*pi/180).*cos.(a_x[kk]*pi/180)
        y=sin.(pi/2 .-a_y[kk]*pi/180).*sin.(a_x[kk]*pi/180)
        z=cos.(pi/2 .-a_y[kk]*pi/180);

        #vector of target points in Cartesian, 3D, coordinates
        xx=sin.(pi/2 .-yvec*pi/180).*cos.(xvec*pi/180);
        yy=sin.(pi/2 .-yvec*pi/180).*sin.(xvec*pi/180);
        zz=cos.(pi/2 .-yvec*pi/180);

        #define tree
        kdtree = KDTree([x y z]')

        #find nearest neighbors
        idxs, _ = knn(kdtree, [xx yy zz]', k, true)
        idxs=[idxs[i][j] for i=1:length(idxs),j=1:k]

        return a_f[kk[idxs]],a_i[kk[idxs]],a_j[kk[idxs]],kk[idxs]
end

knn(xgrid::MeshArray,ygrid::MeshArray,lon::Number,lat::Number) = knn(xgrid::MeshArray,ygrid::MeshArray,[lon],[lat])

"""
    Interpolate(z_in::MeshArray,f,i,j,w)

```
using MeshArrays

γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
Γ=GridLoad(γ; option="full")

lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,vec(lon),vec(lat))
DD=Interpolate(Γ.Depth,f,i,j,w)

using CairoMakie
heatmap(vec(lon[:,1]),vec(lat[1,:]),DD,colorrange=(0.,6000.))
```
"""
function Interpolate(z_in::MeshArray,f,i,j,w)
    z_out=NaN*similar(f[:,1])
    for jj=1:size(f,1)
        if !isnan(sum(w[jj,:]))
            x=[z_in[f[jj,ii]][i[jj,ii],j[jj,ii]] for ii=1:4]
            kk=findall(isfinite.(x))
            ~isempty(kk) ? z_out[jj]=sum(w[jj,kk].*x[kk])/sum(w[jj,kk]) : nothing
        end
    end
    return z_out
end

"""
    Interpolate(z_in::MeshArray,f,i,j,w)

```
using MeshArrays

γ=GridSpec("LatLonCap",MeshArrays.GRID_LLC90)
Γ=GridLoad(γ; option="full")

lon=[i for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
lat=[j for i=-179.5:1.0:179.5, j=-89.5:1.0:89.5]
(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,vec(lon),vec(lat))
L=(lon=lon, lat=lat, f=f, i=i, j=j, w=w)

using CairoMakie
heatmap(Interpolate(Γ.Depth,L)...,colorrange=(0.,6000.))
```

or 

```
heatmap(Γ.Depth,interpolation=L,colorrange=(0.,6000.))
```
"""
function Interpolate(z_in::MeshArray,λ::NamedTuple)
        z=Interpolate(z_in,λ.f,λ.i,λ.j,λ.w)
        z=reshape(z,size(λ.lon))
        return λ.lon[:,1],λ.lat[1,:],z
end

"""
    InterpolationFactors(Γ,lon::Array{T,1},lat::Array{T,1})

Compute interpolation coefficients etc from grid `Γ` to `lon,lat`

```jldoctest; output = false
using MeshArrays
γ=GridSpec("CubeSphere",MeshArrays.GRID_CS32)
Γ=GridLoad(γ; option="full")
lon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)
(f,i,j,w,j_f,j_x,j_y)=InterpolationFactors(Γ,lon,lat)
YC=Interpolate(Γ.YC,f,i,j,w)
extrema(i)==(9,10)

# output

true
```
"""
function InterpolationFactors(Γ,lon::Array{T,1},lat::Array{T,1}) where {T}
#to match gcmfaces test case (`interp=gcmfaces_interp_coeffs(0.1,0.1);`)
#set: iiTile=17; XC0=6.5000; YC0=-0.1994

        #main loop
        i_f=fill(0,length(lon),4)
        i_i=fill(0,length(lon),4)
        i_j=fill(0,length(lon),4)
        i_w=fill(NaN,length(lon),4)
        j_f=fill(0,length(lon),1)
        j_x=fill(0.0,length(lon),1)
        j_y=fill(0.0,length(lon),1)

        #ancillary variables
        (f,i,j,c)=knn(Γ.XC,Γ.YC,lon,lat)

        #1. tiles and τ        
        fs=Γ.XC.fSize
        s=fill(0,2*length(fs))
        [s[collect(1:2) .+ (i-1)*2]=collect(fs[i]) for i in 1:length(fs)]
        ni=gcd(s); nj=gcd(s); γ=Γ.XC.grid
        τ=Tiles(γ,ni,nj); tiles=MeshArray(γ,Int);
        [tiles[τ[ii].face][τ[ii].i,τ[ii].j].=ii for ii in 1:length(τ)]

        #2. t_XC, t_XC, t_f, t_i, t_j        
        t=vec(write(tiles)[c])
        t_list=unique(t)
        t_XC=Tiles(τ,exchange(Γ.XC))
        t_YC=Tiles(τ,exchange(Γ.YC))

        t_f=MeshArray(γ,Int); [t_f[ii][:,:].=ii for ii=1:γ.nFaces]
        t_i=MeshArray(γ,Int); [t_i[ii]=collect(1:γ.fSize[ii][1])*ones(Int,1,γ.fSize[ii][2]) for ii=1:γ.nFaces]
        t_j=MeshArray(γ,Int); [t_j[ii]=ones(Int,γ.fSize[ii][1],1)*collect(1:γ.fSize[ii][2])' for ii=1:γ.nFaces]

        t_f=Tiles(τ,exchange(t_f))
        t_i=Tiles(τ,exchange(t_i))
        t_j=Tiles(τ,exchange(t_j))

        x_q=fill(0.0,1,4)
        y_q=fill(0.0,1,4)
        tmpx=fill(0.0,1,4)
        tmpy=fill(0.0,1,4)
        w=fill(0.0,1,4)
        angsum=fill(0.0,prod(size(t_f[1])))

        #main loop
        for ii=1:length(t_list)
                tt=t_list[ii]

                ff=τ[tt].face
                ii0=minimum(τ[tt].i)+Int(ni/2)
                jj0=minimum(τ[tt].j)+Int(nj/2)
                XC0=Γ.XG.f[ff][ii0,jj0]
                YC0=Γ.YG.f[ff][ii0,jj0]
                #
                (x_grid,y_grid)=StereographicProjection(XC0,YC0,t_XC[tt],t_YC[tt])
                (x_quad,y_quad,i_quad,j_quad)=QuadArrays(x_grid,y_grid)
                #
                x=minimum(τ[tt].i)-0.5 .+collect(-1:ni)*ones(Int,1,nj+2)
                y=minimum(τ[tt].j)-0.5 .+ones(Int,ni+2,1)*collect(-1:nj)'

                for pp in findall(t.==tt)
                        (x_trgt,y_trgt)=StereographicProjection(XC0,YC0,lon[pp],lat[pp])
                        PolygonAngle(x_quad,y_quad,x_trgt,y_trgt,angsum)
                        if sum(angsum.>180.)>0
                                kk=findall(angsum.>180.)[end]
                                i_f[pp,:].=[t_f[tt][i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                i_i[pp,:].=[t_i[tt][i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                i_j[pp,:].=[t_j[tt][i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                x_q[:].=x_quad[kk,:]
                                y_q[:].=y_quad[kk,:]
                                QuadCoeffs(x_q,y_q,x_trgt,y_trgt,w)
                                i_w[pp,:].=w[:]
                                #
                                [tmpx[i]=x[i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                [tmpy[i]=y[i_quad[kk,i]+1,j_quad[kk,i]+1] for i=1:4]
                                #
                                j_f[pp]=ff
                                j_x[pp]=sum(tmpx[:].*i_w[pp,:])
                                j_y[pp]=sum(tmpy[:].*i_w[pp,:])
                        end
                end
        end

        return i_f,i_i,i_j,i_w,j_f,j_x,j_y
end

InterpolationFactors(Γ,lon::Number,lat::Number) = InterpolationFactors(Γ,[lon],[lat])

"""
    interpolation_setup(fil::String)

Read e.g. `interp_coeffs_halfdeg.jld2`
    
```
fil=joinpath(tempdir(),"interp_coeffs_halfdeg.jld2")
λ=MeshArrays.interpolation_setup(fil)
```
"""
function interpolation_setup(fil::String)
        λ = MeshArrays.read_JLD2(fil)
        λ = MeshArrays.Dict_to_NamedTuple(λ)
end

"""
    interpolation_setup(;Γ,lon,lat,path,url)
    
Download or recompute interpolation coefficients.

- `λ=interpolation_setup()` to download "interp_coeffs_halfdeg.jld2" 
- `λ=interpolation_setup(Γ=Γ)` to recompute interpolation to `lon,lat`
"""
function interpolation_setup(;Γ=NamedTuple(),
        lon=[i for i=-179.:2.0:179., j=-89.:2.0:89.],
        lat=[j for i=-179.:2.0:179., j=-89.:2.0:89.])

        if isempty(Γ)
                fil=joinpath(MeshArrays.MA_datadep("interp_halfdeg"),"interp_coeffs_halfdeg.jld2")
        else
		(f,i,j,w)=InterpolationFactors(Γ,vec(lon),vec(lat))
                fil=tempname()*"_interp_coeffs.jld2"
		MeshArrays.write_JLD2(fil; lon=lon, lat=lat, f=f, i=i, j=j, w=w)
        end
        interpolation_setup(fil)
end

##

"""
    StereographicProjection(XC0,YC0,XC,YC)

Apply stereographic projection that puts `XC0,YC0` at `0.0,0.0`
to target point(s) `XC,YC`

```
lon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)
x,y=StereographicProjection(45.,60.,lon,lat)
```
"""
function StereographicProjection(XC0::Number,YC0::Number,XC,YC)
        #compute spherical coordinates:
        phi=XC; theta=90 .-YC;
        phi0=XC0; theta0=90-YC0;

        #compute cartesian coordinates:
        X=sind.(theta).*cosd.(phi);
        Y=sind.(theta).*sind.(phi);
        Z=cosd.(theta);

        x=X; y=Y; z=Z;

        #bring chosen point to the north pole:
        xx=x; yy=y; zz=z;
        x=cosd(phi0).*xx+sind(phi0).*yy;
        y=-sind(phi0).*xx+cosd(phi0).*yy;
        z=zz;

        xx=x; yy=y; zz=z;
        x=cosd(theta0)*xx-sind(theta0)*zz;
        y=yy;
        z=sind(theta0)*xx+cosd(theta0)*zz;

        #stereographic projection from the south pole:
        xx=x./(1 .+z);
        yy=y./(1 .+z);

        #nrm=sqrt(xx.^2+yy.^2);
        #msk=1+0*nrm; msk(nrm>tan(pi/4/2))=NaN;%mask points outside of pi/4 cone

        return xx,yy
end

"""
    PolygonAngle(px::Array,py::Array,angsum::Array)
"""
function PolygonAngle(px::Array,py::Array,angsum::Array)
        M=size(px,1)
        N=size(px,2)
        angsum .= 0.0
        for ii=0:N-1
                i1=mod1(ii+1,N)
                i2=mod1(ii+2,N)
                i4=mod1(ii+4,N)
                for jj=1:M
                #compute sum of sector angles
                v1x=px[jj,i2]-px[jj,i1]
                v1y=py[jj,i2]-py[jj,i1]
                v2x=px[jj,i4]-px[jj,i1]
                v2y=py[jj,i4]-py[jj,i1]

                tmp=( v1x.*v2x+v1y.*v2y )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                g_acos=acos.( min.(max.(tmp,-1.0),1.0) )
                g_sin= ( v1x.*v2y-v1y.*v2x )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                angsum[jj] += rad2deg(g_acos*sign(g_sin))
                end
        end
end

"""
    function PolygonAngle(px::Array,py::Array,x::Array,y::Array,angsum)

Compute sum of interior angles (`angsum`) for polygons or points-to-polygons (when
`px,py,x,y,angsum` is provided as input). `px,py` are `MxN` matrices where each line
specifies one polygon. (optional) `x,y` are position vectors.

```jldoctest; output = false
using MeshArrays
px=[0. 0. 1. 1.]; py=[0. 1. 1. 0.];
x=collect(-1.0:0.25:2.0); y=x;
angsum=fill(0.0,1,length(x))
MeshArrays.PolygonAngle(px,py,x,y,angsum)

isa(angsum,Array)

# output

true
```

"""
function PolygonAngle(px::Array,py::Array,x::Array,y::Array,angsum)
        for ii in 1:length(x)
                PolygonAngle(px,py,x[ii],y[ii],view(angsum,:,ii))
        end        
end

"""
    PolygonAngle(px::Array,py::Array,x::Number,y::Number,angsum)
"""
function PolygonAngle(px::Array,py::Array,x::Number,y::Number,angsum)
        M=size(px,1)
        N=size(px,2)
        angsum .= 0.0
        for ii=0:N-1
                i1=mod1(ii+1,N)
                i2=mod1(ii+2,N)
                for jj=1:M
                #compute sum of sector angles
                v1x=px[jj,i1] -x
                v1y=py[jj,i1] -y
                v2x=px[jj,i2] -x
                v2y=py[jj,i2] -y

                tmp=( v1x.*v2x+v1y.*v2y )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                g_acos=acos.( min.(max.(tmp,-1.0),1.0) )
                g_sin= ( v1x.*v2y-v1y.*v2x )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                angsum[jj] += rad2deg(g_acos*sign(g_sin))
                end
        end
end


"""
    QuadArrays(x_grid,y_grid)

Transform x_grid,y_grid (size ni+2,nj+2) into x_quad,y_quad,i_quad,j_quad
quadrilaterals (size ni+1*nj+1,4) where i_quad,j_quad are point indices
"""
function QuadArrays(x_grid::Array{T,2},y_grid::Array{T,2}) where {T}
        ni,nj=size(x_grid) .-2

        x_quad=Array{Float64,2}(undef,(ni+1)*(nj+1),4)
        y_quad=Array{Float64,2}(undef,(ni+1)*(nj+1),4)
        i_quad=Array{Int64,2}(undef,(ni+1)*(nj+1),4)
        j_quad=Array{Int64,2}(undef,(ni+1)*(nj+1),4)

        didj=[[0 0];[1 0];[1 1];[0 1]]
        for pp=1:4
                di=didj[pp,1]
                dj=didj[pp,2]

                #note the shift in indices due to exchange above
                tmp=x_grid[1+di:ni+1+di,1+dj:nj+1+dj]
                x_quad[:,pp]=vec(tmp)
                tmp=y_grid[1+di:ni+1+di,1+dj:nj+1+dj]
                y_quad[:,pp]=vec(tmp)

                tmp=collect(0+di:ni+di)*ones(1,nj+1)
                i_quad[:,pp]=vec(tmp)
                tmp=ones(ni+1,1)*transpose(collect(0+dj:nj+dj));
                j_quad[:,pp]=vec(tmp)
        end

        return x_quad,y_quad,i_quad,j_quad
end

"""
    QuadCoeffs(px,py,ox=[],oy=[],ow=[])

Compute bilinear interpolation coefficients for `ox,oy` within `px,py`
by remapping these quadrilaterals to the `unit square`.
- `px,py` are `Mx4` matrices where each line specifies one quadrilateral.
- `ox,oy` are `MxP` position matrices
- `ow` (output) are the `MxPx4` bilinear interpolation weights

```
QuadCoeffs([-1., 8., 13., -4.]',[-1., 3., 11., 8.]',0.,6.)
QuadCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)
```
"""
function QuadCoeffs(px,py,ox=[],oy=[],ow=[])
        #test case from https://www.particleincell.com/2012/quad-interpolation/
        #  QuadCoeffs([-1., 8., 13., -4.]',[-1., 3., 11., 8.]',0.,6.)
        #However the case of a perfect parallelogram needs special treatment
        #  QuadCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)
        #Deals with this situtation by falling back to ParaCoeffs
        #  ParaCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)

        #1. solve linear problem (`a,b` vectors from `px,py`)
        #  A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0]; AI = inv(A);
        #  AI=[1 0 0 0;-1 1 0 0;-1 0 0 1; 1 -1 1 -1];
        #  a = AI*px';
        #  b = AI*py';
        #This defines the mapping from logical `l,m` to physical `x,y` as
        #  x=a(1)+a(2)*l+a(3)*m+a(2)*l*m;
        #  y=b(1)+b(2)*l+b(3)*m+b(2)*l*m;

        a=[px[1] -px[1]+px[2] -px[1]+px[4] px[1]-px[2]+px[3]-px[4]]
        a[findall(abs.(a).<1e-8)].=0.0

        b=[py[1] -py[1]+py[2] -py[1]+py[4] py[1]-py[2]+py[3]-py[4]]
        b[findall(abs.(b).<1e-8)].=0.0

        #2. select between the two solutions (to 2nd order
        #non-linear problem below) using polygon interior angles
        angsum=fill(0.0,1)
        PolygonAngle(px,py,angsum)
        sgn=fill(NaN,1)
        isapprox(angsum[1],360.0,rtol=0.01) ? sgn[1]=1.0 : nothing
        isapprox(angsum[1],-360.0,rtol=0.01) ? sgn[1]=-1.0 : nothing
        #ii=findall(isnan.(angsum))
        #length(ii)>0 ? println("warning: edge point was found") : nothing

        #3. solve non-linear problem for `pl,pm` from `px,py` & `a,b`
        # This defines the mapping from physical `x,y` to logical `l,m`

        a=reshape(a,(size(a,1),1,size(a,2))); 
        b=reshape(b,(size(b,1),1,size(b,2))); 

        # quadratic equation coeffs, `aa*mm^2+bb*m+cc=0`
        if ~isempty(ox);
                x=ox; y=oy;
        else;
                x=px; y=py;
                a=repeat(a,1,size(x,2),1);
                b=repeat(b,1,size(x,2),1);
                sgn=repeat(sgn,1,size(x,2));
        end;
        #
        det=fill(0.0,size(x))
        pm=fill(0.0,size(x))
        pl=fill(0.0,size(x))
        for ii=1:size(x,1), jj=1:size(x,2)
                aa = a[ii,jj,4]*b[ii,jj,3]-a[ii,jj,3]*b[ii,jj,4]
                bb = a[ii,jj,4]*b[ii,jj,1]-a[ii,jj,1]*b[ii,jj,4]+a[ii,jj,2]*b[ii,jj,3]-a[ii,jj,3]*b[ii,jj,2]+x[ii,jj]*b[ii,jj,4]-y[ii,jj].*a[ii,jj,4]
                cc = a[ii,jj,2]*b[ii,jj,1]-a[ii,jj,1]*b[ii,jj,2]+x[ii,jj]*b[ii,jj,2]-y[ii,jj].*a[ii,jj,2]
                #compute `pm = (-b+sqrt(b^2-4ac))/(2a)`
                det[ii,jj] = sqrt(bb*bb - 4.0*aa*cc)
                pm[ii,jj]  = (-bb+sgn[ii,jj]*det[ii,jj])/(2.0*aa)
                #compute `pl` by substitution in equation system
                pl[ii,jj]  = (x[ii,jj]-a[ii,jj,1]-a[ii,jj,3]*pm[ii,jj])/(a[ii,jj,2]+a[ii,jj,4]*pm[ii,jj])
        end

        if ~isempty(ox);
                tmp1=(1 .-pl).*(1 .-pm)
                tmp2=pl.*(1 .-pm)
                tmp3=pl.*pm
                tmp4=(1 .-pl).*pm
                ow[:]=cat(tmp1,tmp2,tmp3,tmp4; dims=3)
                #ow[:].=cat(tmp1,tmp2,tmp3,tmp4; dims=3)[:]
                #treat pathological cases if needed
                tmp=ParaCoeffs(px,py,[ox],[oy])
                for kk=1:length(ow)
                        !isfinite(ow[kk]) ? ow[kk]=tmp[kk] : nothing
                end
        end

end

"""
    ParaCoeffs(px,py,ox=[],oy=[])

Compute bilinear interpolation coefficients for `ox,oy` within `px,py`
by remapping these parallelograms to the `unit square`.
- `px,py` are `Mx4` matrices where each line specifies one quadrilateral.
- `ox,oy` are `MxP` position matrices
- `pw` (output) are the `MxPx4` bilinear interpolation weights

```
x=1.0; y=1.0 #Try send the corners to unit square corners?
println(vec(ParaCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',x,y)))
println(vec(QuadCoeffs([0., 2.01, 3., 1.]',[0., 0., 1., 1.]',x,y)))
```
"""
function ParaCoeffs(px,py,ox=[],oy=[])

        tmp1=px[:,1];
        tmp2=-px[:,1]+px[:,2];
        tmp3=-px[:,2]+px[:,3];
        a=[tmp1 tmp2 tmp3];

        tmp1=py[:,1];
        tmp2=-py[:,1]+py[:,2];
        tmp3=-py[:,2]+py[:,3];
        b=[tmp1 tmp2 tmp3];

#        (m,l)=inv([a[1,2] a[1,3];b[1,2] b[1,3]])*[ox[1]-a[1,1]; oy[1]-b[1,1]]
        m=( b[:,3].*(ox-a[:,1])-a[:,3].*(oy-b[:,1]) ) ./(a[:,2].*b[:,3]-a[:,3].*b[:,2])
        l=( -b[:,2].*(ox-a[:,1])+a[:,2].*(oy-b[:,1]) ) ./(a[:,2].*b[:,3]-a[:,3].*b[:,2])

        ow=[];
        tmp1=(1 .-l).*(1 .-m)
        tmp4=l.*(1 .-m)
        tmp3=l.*m
        tmp2=(1 .-l).*m
        ow=cat(tmp1,tmp2,tmp3,tmp4; dims=3)

        return ow
end
