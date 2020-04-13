
using NearestNeighbors
import NearestNeighbors: knn

"""
    knn(xgrid,ygrid::MeshArray,x,y::Array{T,1},k::Int)

Find k nearest neighbors to each point in x,y on xgrid,ygrid

```
lon=collect(0.1:0.5:2.1); lat=collect(0.1:0.5:2.1);
(f,i,j,a)=knn(Γ["XC"],Γ["YC"],lon,lat)
```
"""
function knn(xgrid::MeshArray,ygrid::MeshArray,
        xvec::Array{T,1},yvec::Array{T,1},k=1::Int) where {T}

        #ancillary variables
        Γ=xgrid.grid
        a_f=MeshArray(Γ,Int); [a_f[ii][:,:].=ii for ii=1:Γ.nFaces]
        a_i=MeshArray(Γ,Int); [a_i[ii]=collect(1:Γ.fSize[ii][1])*ones(Int,1,Γ.fSize[ii][2]) for ii=1:Γ.nFaces]
        a_j=MeshArray(Γ,Int); [a_j[ii]=ones(Int,Γ.fSize[ii][1],1)*collect(1:Γ.fSize[ii][2])' for ii=1:Γ.nFaces]

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

"""
    StereographicProjection(XC0,YC0,XC,YC)

Apply stereographic projection that puts `XC0,YC0` at `0.0,0.0`
to target point(s) `XC,YC`

```
lon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)
x,y=StereographicProjection(45.,60.,lon,lat)
```
"""
function StereographicProjection(XC0,YC0,XC,YC)
        #For additional information see e.g.
        #http://www4.ncsu.edu/~franzen/public_html/CH795Z/math/lab_frame/lab_frame.html
        #http://physics.unm.edu/Courses/Finley/p503/handouts/SphereProjectionFINALwFig.pdf

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
    PolygonAngle(px,py,x=[],y=[])

Compute sum of interior angles for polygons or points-to-polygons (when
`px,py,x,y` is provided as input). `px,py` are `MxN` matrices where each line
specifies one polygon. (optional) `x,y` are position vectors.

```
lon=collect(45.:0.1:46.); lat=collect(60.:0.1:61.)
x,y=StereographicProjection(45.,60.,lon,lat)
```
"""
function PolygonAngle(px,py,x=[],y=[])

        M=size(px,1); N=size(px,2); P=1;
        doPointsInPolygon=false
        if length(x)>0;
                doPointsInPolygon=true
                x=reshape(x,1,length(x))
                y=reshape(y,1,length(y))
                P=length(x)
        end;

        angsum=zeros(M,P)
        for ii=0:N-1
                ppx=circshift(px,[0 -ii])
                ppy=circshift(py,[0 -ii])
                if ~doPointsInPolygon
                        #compute sum of corner angles
                        v1x=ppx[:,2]-ppx[:,1]
                        v1y=ppy[:,2]-ppy[:,1]
                        v2x=ppx[:,4]-ppx[:,1]
                        v2y=ppy[:,4]-ppy[:,1]
                else;
                        #compute sum of sector angles
                        v1x=ppx[:,1]*ones(1,P)-ones(M,1)*x
                        v1y=ppy[:,1]*ones(1,P)-ones(M,1)*y
                        v2x=ppx[:,2]*ones(1,P)-ones(M,1)*x
                        v2y=ppy[:,2]*ones(1,P)-ones(M,1)*y
                end
                g_acos=acos.( ( v1x.*v2x+v1y.*v2y )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y ) )
                g_sin= ( v1x.*v2y-v1y.*v2x )./sqrt.( v1x.*v1x+v1y.*v1y )./sqrt.( v2x.*v2x+v2y.*v2y )
                angsum=angsum+rad2deg.(g_acos.*sign.(g_sin));
        end;

        return angsum
end


"""
    QuadCoeffs(px,py,ox=[],oy=[])

Compute bilinear interpolation coefficients for `ox,oy` within `px,py`
by remapping these quadrilaterals to the `unit square`.
- `px,py` are `Mx4` matrices where each line specifies one quadrilateral.
- `ox,oy` are `MxP` position matrices
- `pw` (output) are the `MxPx4` bilinear interpolation weights

```
QuadCoeffs([-1., 8., 13., -4.]',[-1., 3., 11., 8.]',0.,6.)
QuadCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)
```
"""
function QuadCoeffs(px,py,ox=[],oy=[])
        #test case from https://www.particleincell.com/2012/quad-interpolation/
        #  QuadCoeffs([-1., 8., 13., -4.]',[-1., 3., 11., 8.]',0.,6.)
        #However the case of a perfect parallelogram needs special treatment
        #  QuadCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)
        #Deals with this situtation by falling back to ParalCoeffs
        #  ParalCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',0.1,0.1)

        #1. solve linear problem (`a,b` vectors from `px,py`)
        #  A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0]; AI = inv(A);
        #  AI=[1 0 0 0;-1 1 0 0;-1 0 0 1; 1 -1 1 -1];
        #  a = AI*px';
        #  b = AI*py';
        #This defines the mapping from logical `l,m` to physical `x,y` as
        #  x=a(1)+a(2)*l+a(3)*m+a(2)*l*m;
        #  y=b(1)+b(2)*l+b(3)*m+b(2)*l*m;

        tmp1=px[:,1];
        tmp2=-px[:,1]+px[:,2];
        tmp3=-px[:,1]+px[:,4];
        tmp4=px[:,1]-px[:,2]+px[:,3]-px[:,4];
        a=[tmp1 tmp2 tmp3 tmp4];

        tmp1=py[:,1];
        tmp2=-py[:,1]+py[:,2];
        tmp3=-py[:,1]+py[:,4];
        tmp4=py[:,1]-py[:,2]+py[:,3]-py[:,4];
        b=[tmp1 tmp2 tmp3 tmp4];

        #2. select between the two solutions (to 2nd order
        #non-linear problem below) using polygon interior angles
        angsum=PolygonAngle(px,py)
        sgn=NaN*px[:,1];
        ii=findall(abs.(angsum .-360).<1e-3); sgn[ii].=1.
        ii=findall(abs.(angsum .+360).<1e-3); sgn[ii].=-1.
        ii=findall(isnan.(angsum))
        length(ii)>0 ? warning("edge point was found") : nothing

        #3. solve non-linear problem for `pl,pm` from `px,py` & `a,b`
        # This defines the mapping from physical `x,y` to logical `l,m`

        # quadratic equation coeffs, `aa*mm^2+bb*m+cc=0`
        if ~isempty(ox);
                x=[px ox]; y=[py oy];
        else;
                x=px; y=py;
        end;
        a=reshape(a,(size(a,1),1,size(a,2))); a=repeat(a,1,size(x,2),1);
        b=reshape(b,(size(b,1),1,size(b,2))); b=repeat(b,1,size(x,2),1);
        sgn=repeat(sgn,1,size(x,2));
        #
        aa = a[:,:,4].*b[:,:,3] -a[:,:,3].*b[:,:,4]
        bb = a[:,:,4].*b[:,:,1] -a[:,:,1].*b[:,:,4] + a[:,:,2].*b[:,:,3] - a[:,:,3].*b[:,:,2] + x.*b[:,:,4] - y.*a[:,:,4]
        cc = a[:,:,2].*b[:,:,1] -a[:,:,1].*b[:,:,2] + x.*b[:,:,2] - y.*a[:,:,2]

        #compute `pm = (-b+sqrt(b^2-4ac))/(2a)`
        det = sqrt.(bb.*bb - 4.0*aa.*cc)
        pm = (-bb+sgn.*det)./(2.0*aa)
        #compute `pl` by substitution in equation system
        pl = (x-a[:,:,1]-a[:,:,3].*pm)./(a[:,:,2]+a[:,:,4].*pm)

        ow=[];
        if ~isempty(ox);
                tmp1=(1 .-pl[:,5:end]).*(1 .-pm[:,5:end])
                tmp2=pl[:,5:end].*(1 .-pm[:,5:end])
                tmp3=pl[:,5:end].*pm[:,5:end]
                tmp4=(1 .-pl[:,5:end]).*pm[:,5:end]
                ow=cat(tmp1,tmp2,tmp3,tmp4; dims=3)
                #treat pathological cases if needed
                tmp=ParalCoeffs(px,py,ox,oy)
                ow[findall(isnan.(ow))].=tmp[findall(isnan.(ow))]
        end

        return ow
end

"""
    ParalCoeffs(px,py,ox=[],oy=[])

Compute bilinear interpolation coefficients for `ox,oy` within `px,py`
by remapping these parallelograms to the `unit square`.
- `px,py` are `Mx4` matrices where each line specifies one quadrilateral.
- `ox,oy` are `MxP` position matrices
- `pw` (output) are the `MxPx4` bilinear interpolation weights

```
x=1.0; y=1.0 #Try send the corners to unit square corners?
println(vec(ParalCoeffs([0., 2., 3., 1.]',[0., 0., 1., 1.]',x,y)))
println(vec(QuadCoeffs([0., 2.01, 3., 1.]',[0., 0., 1., 1.]',x,y)))
```
"""
function ParalCoeffs(px,py,ox=[],oy=[])

        tmp1=px[:,1];
        tmp2=-px[:,1]+px[:,2];
        tmp3=-px[:,2]+px[:,3];
        a=[tmp1 tmp2 tmp3];

        tmp1=py[:,1];
        tmp2=-py[:,1]+py[:,2];
        tmp3=-py[:,2]+py[:,3];
        b=[tmp1 tmp2 tmp3];

        (m,l)=inv([a[1,2] a[1,3];b[1,2] b[1,3]])*[ox[1]-a[1,1]; oy[1]-b[1,1]]

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
