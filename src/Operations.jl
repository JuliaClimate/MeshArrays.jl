
## basic handling of NaNs that are commonly used in masks

nansum(x) = sum(filter(!isnan,x))
nansum(x,y) = mapslices(nansum,x,dims=y)
nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)
nanmax(x) = maximum(filter(!isnan,x))
nanmax(x,y) = mapslices(nanmax,x,dims=y)
nanmin(x) = minimum(filter(!isnan,x))
nanmin(x,y) = mapslices(nanmin,x,dims=y)

function nanmean(a::Number,b::Number)
	if isnan(a)&&isnan(b)
		NaN
	elseif isnan(a)
		b
	elseif isnan(b)
		a
	else
		(a+b)/2
	end
end

## gradient methods

"""
    gradient(inFLD::AbstractMeshArray,Γ::NamedTuple)

Compute spatial derivatives. Other methods:
```
gradient(inFLD::AbstractMeshArray,Γ::NamedTuple,doDIV::Bool)
gradient(inFLD::AbstractMeshArray,iDXC::AbstractMeshArray,iDYC::AbstractMeshArray)
```
"""
function gradient(inFLD::AbstractMeshArray,Γ::NamedTuple)
(dFLDdx, dFLDdy)=gradient(inFLD,Γ,true)
return dFLDdx, dFLDdy
end

function gradient(inFLD::AbstractMeshArray,Γ::NamedTuple,doDIV::Bool)

exFLD=exchange(inFLD).MA
dFLDdx=similar(inFLD)
dFLDdy=similar(inFLD)

for a=1:inFLD.grid.nFaces
  (s1,s2)=size(exFLD.f[a])
  tmpA=view(exFLD.f[a],2:s1-1,2:s2-1)
  tmpB=tmpA-view(exFLD.f[a],1:s1-2,2:s2-1)
  tmpC=tmpA-view(exFLD.f[a],2:s1-1,1:s2-2)
  if doDIV
    dFLDdx.f[a]=tmpB./Γ.DXC.f[a]
    dFLDdy.f[a]=tmpC./Γ.DYC.f[a]
  else
    dFLDdx.f[a]=tmpB
    dFLDdy.f[a]=tmpC
  end
end

return dFLDdx, dFLDdy
end

function gradient(inFLD::AbstractMeshArray,iDXC::AbstractMeshArray,iDYC::AbstractMeshArray)

exFLD=exchange(inFLD).MA
dFLDdx=similar(inFLD)
dFLDdy=similar(inFLD)

for a=1:inFLD.grid.nFaces
  (s1,s2)=size(exFLD.f[a])
  tmpA=view(exFLD.f[a],2:s1-1,2:s2-1)
  tmpB=tmpA-view(exFLD.f[a],1:s1-2,2:s2-1)
  tmpC=tmpA-view(exFLD.f[a],2:s1-1,1:s2-2)
  dFLDdx.f[a]=tmpB.*iDXC.f[a]
  dFLDdy.f[a]=tmpC.*iDYC.f[a]
end

return dFLDdx, dFLDdy
end

##

to_UV(inFLD::AbstractMeshArray) = 
  length(size(inFLD))==1 ? to_UV_2d(inFLD) : to_UV_3d(inFLD)

function to_UV_3d(inFLD::AbstractMeshArray)
	TatU=similar(inFLD)
	TatV=similar(inFLD)
  to_UV_3d!(inFLD,TatU,TatV)
	TatU,TatV
end

function to_UV_3d!(inFLD::AbstractMeshArray,
  TatU::AbstractMeshArray,TatV::AbstractMeshArray;
  verbose=false)
  tmpU=similar(inFLD[:,1])
  tmpV=similar(inFLD[:,1])
  verbose ? println(size(inFLD)[2]) : nothing
  for k in 1:size(inFLD)[2]
    verbose ? println(k) : nothing
    exFLD=exchange(inFLD[:,k])
    to_UV_2d!(exFLD,tmpU,tmpV)
    for i in 1:size(inFLD)[1]
      TatU[i,k].=tmpU[i]
      TatV[i,k].=tmpV[i]
    end
  end
end

function to_UV_2d(inFLD::AbstractMeshArray)
	TatU=similar(inFLD)
	TatV=similar(inFLD)
  to_UV_2d!(inFLD,TatU,TatV)
	return TatU, TatV
end

function to_UV_2d!(inFLD::AbstractMeshArray,TatU::AbstractMeshArray,TatV::AbstractMeshArray)
  exFLD=exchange(inFLD)
  to_UV_2d!(exFLD,TatU,TatV)
end

function to_UV_2d!(exFLD::MeshArray_wh,TatU::AbstractMeshArray,TatV::AbstractMeshArray)
	for a=1:exFLD.MA.grid.nFaces
		(s1,s2)=size(exFLD.MA.f[a])
		tmpA=view(exFLD.MA.f[a],2:s1-1,2:s2-1)
		TatU.f[a].=0.5*(tmpA+view(exFLD.MA.f[a],1:s1-2,2:s2-1))
		TatV.f[a].=0.5*(tmpA+view(exFLD.MA.f[a],2:s1-1,1:s2-2))
	end
end

##

"""
    curl(u::AbstractMeshArray,v::AbstractMeshArray,Γ::NamedTuple)

Compute curl of a velocity field.
"""
function curl(u::AbstractMeshArray,v::AbstractMeshArray,Γ::NamedTuple)

	uvcurl=similar(Γ.XC)
	fac=exchange(1.0 ./Γ.RAZ)
	(U,V)=exchange_main(u,v,1)
  (DXC,DYC)=exchange_main(Γ.DXC,Γ.DYC,1)
	[DXC.MA[i].=abs.(DXC.MA[i]) for i in eachindex(U.MA)]
	[DYC.MA[i].=abs.(DYC.MA[i]) for i in eachindex(V.MA)]

	for i in eachindex(U.MA)
    ucur=U.MA[i][2:end,:]
    vcur=V.MA[i][:,2:end]        
    tmpcurl=ucur[:,1:end-1]-ucur[:,2:end]
    tmpcurl=tmpcurl-(vcur[1:end-1,:]-vcur[2:end,:])
    tmpcurl=tmpcurl.*fac.MA[i][1:end-1,1:end-1]

		##still needed:
		##- deal with corners
		##- if putCurlOnTpoints

		tmpcurl=1/4*(tmpcurl[1:end-1,2:end]+tmpcurl[1:end-1,1:end-1]+
					tmpcurl[2:end,2:end]+tmpcurl[2:end,1:end-1])

		uvcurl[i]=tmpcurl
	end
	
	return uvcurl
end



export f, β

Ω = 7.2921 * 10^(-5) #rad/s
R = 6.3781 * 10^(6) #m
#Coriolis Parameter
f_Coriolis(φ) = 2Ω * sind(φ)
#f_Coriolis(φ::AbstractMeshArray) = 2Ω * sind.(φ)

"""
    EkmanTrsp(u::AbstractMeshArray,v::AbstractMeshArray,Γ::NamedTuple)

Compute Ekman Transport (in m2/s) from wind stress (in N/m2 , or kg/m/s2).
"""
function EkmanTrsp(u::AbstractMeshArray,v::AbstractMeshArray,Γ::NamedTuple)
	EkX=similar(Γ.DYG)
	EkY=similar(Γ.DXG)
	(U,V)=exchange_main(u,v,1)

	for i in eachindex(U.MA)
    ucur=1/4* (U.MA[i][2:end-1,1:end-2]+U.MA[i][2:end-1,2:end-1]
              +U.MA[i][3:end,1:end-2]+U.MA[i][3:end,2:end-1])
    fcur=f_Coriolis.(Γ.YS[i])
    fcur[abs.(Γ.YS[i]).<10].=NaN
    EkY[i]=-ucur./fcur./1029.0
    ucur=1/4* (V.MA[i][1:end-2,2:end-1]+V.MA[i][1:end-2,2:end-1]
              +V.MA[i][2:end-1,3:end]+V.MA[i][2:end-1,3:end])
    fcur=f_Coriolis.(Γ.YW[i])
    fcur[abs.(Γ.YW[i]).<10].=NaN
    EkX[i]=ucur./fcur./1029.0
	end

	(EkX,EkY)
end

##

function overturning(Utr,Vtr,LC,Γ,mask)
	nz=size(Γ.hFacC,2); nt=12; nl=length(LC)
	ov=Array{Float64,2}(undef,nl,nz)
	#integrate across latitude circles
	for z=1:nz
		UV=Dict("U"=>mask.*Utr[:,z],"V"=>mask.*Vtr[:,z],"dimensions"=>["x","y"])
		[ov[l,z]=ThroughFlow(UV,LC[l],Γ) for l=1:nl]
	end	
	#integrate from bottom
	ov=reverse(cumsum(reverse(ov,dims=2),dims=2),dims=2)
end

function meridional(Utr,Vtr,LC,Γ,mask)
	nl=length(LC)
	nz=size(Utr,2)
	MT=fill(0.0,nl)
	for z=1:nz
		UV=Dict("U"=>mask*Utr[:,z],"V"=>mask*Vtr[:,z],"dimensions"=>["x","y"])
		[MT[l]=MT[l]+ThroughFlow(UV,LC[l],Γ) for l=1:nl]
	end
	MT
end

## mask methods

function mask(fld::AbstractMeshArray)
fldmsk=mask(fld,NaN)
return fldmsk
end

"""
    mask(fld::AbstractMeshArray, val::Number)

Replace non finite values with val. Other methods:
```
mask(fld::AbstractMeshArray)
mask(fld::AbstractMeshArray, val::Number, noval::Number)
```
"""
function mask(fld::AbstractMeshArray, val::Number)
  fldmsk=similar(fld)
  for a=1:fld.grid.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> !isfinite(x) ? val : x, tmp1 )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

function mask(fld::AbstractMeshArray, val::Number, noval::Number)
  fldmsk=similar(fld)
  for a=1:fld.grid.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> x==noval ? val : x, tmp1  )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

"""
    land_mask(m::AbstractMeshArray)

Define land mask from `m` (1 if m>0; NaN if otherwise).
"""
function land_mask(m::AbstractMeshArray)
    μ=m
    μ[findall(μ.>0.0)].=1.0
    μ[findall(μ.==0.0)].=NaN
    μ
end

land_mask(Γ::NamedTuple)=land_mask(Γ.hFacC[:,1])

## convergence methods

"""
    convergence(uFLD::AbstractMeshArray,vFLD::AbstractMeshArray)

Compute convergence of a vector field
"""
function convergence(uFLD::AbstractMeshArray,vFLD::AbstractMeshArray)

#important note:
#  Normally uFLD, vFLD should not contain any NaN;
#  if otherwise then something this may be needed:
#  uFLD=mask(uFLD,0.0); vFLD=mask(vFLD,0.0);

CONV=similar(uFLD)

(tmpU,tmpV)=exch_UV(uFLD,vFLD)
for a=1:tmpU.grid.nFaces
  (s1,s2)=size(uFLD.f[a])
  tmpU1=view(tmpU.f[a],1:s1,1:s2)
  tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
  tmpV1=view(tmpV.f[a],1:s1,1:s2)
  tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
  CONV.f[a]=tmpU1-tmpU2+tmpV1-tmpV2
end

return CONV
end

## smooth function

"""
    smooth(FLD::AbstractMeshArray,DXCsm::AbstractMeshArray,DYCsm::AbstractMeshArray,Γ::NamedTuple)

Smooth out scales below DXCsm / DYCsm via diffusion
"""
function smooth(FLD::AbstractMeshArray,DXCsm::AbstractMeshArray,DYCsm::AbstractMeshArray,Γ::NamedTuple)

#important note:
#input FLD should be land masked (NaN/1) by caller if needed

#get land masks (NaN/1):
mskC=fill(1.0,FLD) + 0.0 * mask(FLD)
(mskW,mskS)=gradient(FLD,Γ,false)
mskW=fill(1.0,FLD) + 0.0 * mask(mskW)
mskS=fill(1.0,FLD) + 0.0 * mask(mskS)

#replace NaN with 0. in FLD and land masks:
FLD=mask(FLD,0.0)
mskC=mask(mskC,0.0)
mskW=mask(mskW,0.0)
mskS=mask(mskS,0.0)

#get inverse grid spacing:
iDXC=similar(FLD)
iDYC=similar(FLD)
for a=1:FLD.grid.nFaces
  iDXC.f[a]=1.0./Γ.DXC.f[a]
  iDYC.f[a]=1.0./Γ.DYC.f[a]
end

#Before scaling the diffusive operator ...
tmp0=DXCsm*iDXC*mskW;
tmp00=maximum(tmp0);
tmp0=DYCsm*iDYC*mskS;
tmp00=max(tmp00,maximum(tmp0));

#... determine a suitable time period:
nbt=ceil(1.1*2*tmp00^2);
dt=1.;
T=nbt*dt;
#println("nbt="*"$nbt")

#diffusion operator times DYG / DXG
KuxFac=mskW*DXCsm*DXCsm/T/2.0*Γ.DYG
KvyFac=mskS*DYCsm*DYCsm/T/2.0*Γ.DXG

#time steping factor:
dtFac=dt*mskC/Γ.RAC

#loop:
for it=1:nbt
  (dTdxAtU,dTdyAtV)=gradient(FLD,iDXC,iDYC);
  tmpU=similar(FLD)
  tmpV=similar(FLD)
  for a=1:FLD.grid.nFaces
      tmpU.f[a]=dTdxAtU.f[a].*KuxFac.f[a];
      tmpV.f[a]=dTdyAtV.f[a].*KvyFac.f[a];
  end
  tmpC=convergence(tmpU,tmpV);
  for a=1:FLD.grid.nFaces
      FLD.f[a]=FLD.f[a]-dtFac.f[a].*tmpC.f[a];
  end
end

#Apply land mask (NaN/1) to end result:
mskC=mask(mskC,NaN,0.0)
FLD=mskC*FLD

return FLD

end

## ThroughFlow function

"""
    ThroughFlow(VectorField,IntegralPath,Γ::NamedTuple)

Compute transport through an integration path
"""
function ThroughFlow(VectorField,IntegralPath,Γ::NamedTuple,msk=[])

    #Note: vertical intergration is not always wanted; left for user to do outside
    MSK=(isempty(msk) ? 1.0 : msk)
    U=MSK .*VectorField["U"] 
    V=MSK .*VectorField["V"] 

    nd=ndims(U)
    #println("nd=$nd and d=$d")

    n=fill(1,4)
    tmp=size(U)
    n[1:nd].=tmp[1:nd]

    haskey(VectorField,"factors") ? f=VectorField["factors"] : f=Array{String,1}(undef,0)
    haskey(VectorField,"dimensions") ? d=VectorField["dimensions"] : d=Array{String,1}(undef,nd)

    #a bit of a hack to distinguish gcmfaces v gcmarray indexing:
    isdefined(U,:fIndex) ? ndoffset=1 : ndoffset=0
    length(d)!=nd+ndoffset ? error("inconsistent specification of dims") : nothing

    trsp=Array{Float64}(undef,1,n[3],n[4])
    do_dz=sum(f.=="dz")
    do_dxory=sum(f.=="dxory")

    for i3=1:n[3]
        #method 1: quite slow
        #mskW=IntegralPath.mskW
        #do_dxory==1 ? mskW=mskW*Γ.DYG : nothing
        #do_dz==1 ? mskW=Γ.DRF[i3]*mskW : nothing
        #mskS=IntegralPath.mskS
        #do_dxory==1 ? mskS=mskS*Γ.DXG : nothing
        #do_dz==1 ? mskS=Γ.DRF[i3]*mskS : nothing
        #
        #method 2: less slow
        tabW=(isa(IntegralPath,NamedTuple) ? IntegralPath.tabW : IntegralPath.W)
        tabS=(isa(IntegralPath,NamedTuple) ? IntegralPath.tabS : IntegralPath.S)
        for i4=1:n[4]
            #method 1: quite slow
            #trsp[1,i3,i4]=sum(mskW*U[:,:,i3,i4])+sum(mskS*V[:,:,i3,i4])
            #
            #method 2: less slow
            trsp[1,i3,i4]=0.0
            for k=1:size(tabW,1)
                (a,i1,i2,w)=tabW[k,:]
                do_dxory==1 ? w=w*Γ.DYG.f[a][i1,i2] : nothing
                do_dz==1 ? w=w*Γ.DRF[i3] : nothing
                isdefined(U,:fIndex) ? u=U.f[a,i3,i4][i1,i2] : u=U.f[a][i1,i2,i3,i4]
                trsp[1,i3,i4]=trsp[1,i3,i4]+w*u
            end
            for k=1:size(tabS,1)
                (a,i1,i2,w)=tabS[k,:]
                do_dxory==1 ? w=w*Γ.DXG.f[a][i1,i2] : nothing
                do_dz==1 ? w=w*Γ.DRF[i3] : nothing
                isdefined(V,:fIndex) ? v=V.f[a,i3,i4][i1,i2] : v=V.f[a][i1,i2,i3,i4]
                trsp[1,i3,i4]=trsp[1,i3,i4]+w*v
            end
        end
    end

    nd+ndoffset<4 ? trsp=dropdims(trsp,dims=3) : nothing
    nd+ndoffset<3 ? trsp=dropdims(trsp,dims=2) : nothing
    nd+ndoffset==2 ? trsp=trsp[1] : nothing

    return trsp
end

## LatitudeCircles function

"""
    LatitudeCircles(LatValues,Γ::NamedTuple; format=:gridpath, range=(0.0,360.0))

Compute integration paths that follow latitude circles, within the specified longitude `range`.
"""
function LatitudeCircles(LatValues,Γ::NamedTuple; 
  format=:gridpath, range=(0.0,360.0))
  T=(format==:NamedTuple ? NamedTuple : gridpath)
  LatitudeCircles=Array{T}(undef,length(LatValues))
    for j=1:length(LatValues)
        LatitudeCircles[j]=LatitudeCircle(LatValues[j],Γ; 
        format=format,range=range)
    end
    (length(LatValues)==1 ? LatitudeCircles[1] : LatitudeCircles)
end

function LatitudeCircle(lat,Γ::NamedTuple; 
  format=:gridpath, range=(0.0,360.0))
      mskCint=1*(Γ.YC .>= lat)
      mskC,mskW,mskS=edge_mask(mskCint)
      restrict_longitudes!(mskC,Γ.XC,range=range)
      restrict_longitudes!(mskS,Γ.XS,range=range)
      restrict_longitudes!(mskW,Γ.XW,range=range)
      LC=if format==:NamedTuple
        (lat=LatValues[j],tabC=MskToTab(mskC),
        tabW=MskToTab(mskW),tabS=MskToTab(mskS))
      else
        gridpath(name="Parallel $lat", grid=Γ,
        C=MskToTab(mskC),W=MskToTab(mskW),S=MskToTab(mskS))
      end

end

is_in_lon_range(x,range)=(range[2].-range[1]>=360)||
  (mod(x-range[1],360).<mod(range[2].-range[1],360))

function restrict_longitudes!(x::AbstractMeshArray,lon::AbstractMeshArray;range=(0.0,360.0))
  for f in 1:x.grid.nFaces
    x[f].=x[f].*is_in_lon_range.(lon[f],Ref(range))
  end
end

function MskToTab(msk::AbstractMeshArray)
  n=Int(sum(msk .!= 0)); k=0
  tab=Array{Int,2}(undef,n,4)
  for i in eachindex(msk)
    a=msk[i]
    b=findall( a .!= 0)
    for ii in eachindex(b)
      k += 1
      tab[k,:]=[i,b[ii][1],b[ii][2],a[b[ii]]]
    end
  end
  return tab
end

"""
    edge_mask(mskCint::AbstractMeshArray)

Compute edge mask (mskC,mskW,mskS) from domain interior mask (mskCint). 
This is used in `LatitudeCircles` and `Transect`.
"""
function edge_mask(mskCint::AbstractMeshArray)
  mskCint=1.0*mskCint

  #treat the case of blank tiles:
  #mskCint[findall(RAC.==0)].=NaN
  
  mskC=similar(mskCint)
  mskW=similar(mskCint)
  mskS=similar(mskCint)

  mskCint=exchange(mskCint).MA

  for i in eachindex(mskCint)
      tmp1=mskCint[i]
      # tracer mask:
      tmp2=tmp1[2:end-1,1:end-2]+tmp1[2:end-1,3:end]+
      tmp1[1:end-2,2:end-1]+tmp1[3:end,2:end-1]
      mskC[i]=1((tmp2.>0).&(tmp1[2:end-1,2:end-1].==0))
      # velocity masks:
      mskW[i]=tmp1[2:end-1,2:end-1] - tmp1[1:end-2,2:end-1]
      mskS[i]=tmp1[2:end-1,2:end-1] - tmp1[2:end-1,1:end-2]
  end

  #treat the case of blank tiles:
  #mskC[findall(isnan.(mskC))].=0.0
  #mskW[findall(isnan.(mskW))].=0.0
  #mskS[findall(isnan.(mskS))].=0.0

  return mskC,mskW,mskS
end

##

"""
    Transect(name,lons,lats,Γ; segment=:short, format=:gridpath)

Compute integration paths that follow a great circle between two geolocations give by `lons`, `lats`.

"""
function Transect(name,lons,lats,Γ; segment=:short, format=:gridpath)
  x0,y0,z0,R=rotate_points(lons,lats)
  x,y,z=rotate_XCYC(Γ,R)
  mskCint=1.0*(z.>0)
  mskCedge,mskWedge,mskSedge=edge_mask(mskCint)

  mskCedge,mskWedge,mskSedge=
  if segment==:short
    shorter_paths!((x,y,z),(x0,y0,z0),(mskCedge,mskWedge,mskSedge))
  elseif segment==:long
    C,W,S=shorter_paths!((x,y,z),(x0,y0,z0),(mskCedge,mskWedge,mskSedge))
    C-mskCedge,W-mskWedge,S-mskSedge
  elseif segment==:circle
    mskCedge,mskWedge,mskSedge
  end
  
  tabC=MskToTab(mskCedge)
  tabW=MskToTab(mskWedge)
  tabS=MskToTab(mskSedge)
  
  if format==:NamedTuple
    (name=name,tabC=tabC,tabW=tabW,tabS=tabS)
  else
    gridpath(name=name,grid=Γ,C=tabC,W=tabW,S=tabS)
  end
end


"""
    edge_path(name,mskCint,Γ)

Compute integration path that follows the outer edge of `mskCint>0`. 

```
Γ=GridLoad(ID=:LLC90)
mask=demo.extended_basin(demo.ocean_basins(),:Pac)
edge=edge_path("Pacific Ocean Edge",mask,Γ)
```
"""
function edge_path(name,mskCint,Γ)
  mskCedge,mskWedge,mskSedge=edge_mask(mskCint)
  tabC=MskToTab(mskCedge)
  tabW=MskToTab(mskWedge)
  tabS=MskToTab(mskSedge)  
  gridpath(name=name,grid=Γ,C=tabC,W=tabW,S=tabS)
end

##

"""
    UVtoUEVN(u,v,G::NamedTuple)

1. Interpolate to grid cell centers (uC,vC)
2. Convert to `Eastward/Northward` components (uE,vN)

Note: land masking `u,v` with `NaN`s preemptively can be adequate.
"""
function UVtoUEVN(u::AbstractMeshArray,v::AbstractMeshArray,G::NamedTuple)
    #u[findall(G.hFacW[:,1].==0)].=NaN
    #v[findall(G.hFacS[:,1].==0)].=NaN

    (u,v)=exch_UV(u,v); uC=similar(u); vC=similar(v)
    for iF=1:u.grid.nFaces
        tmp1=u[iF][1:end-1,:]; tmp2=u[iF][2:end,:]
        uC[iF]=reshape(nanmean([tmp1[:] tmp2[:]],2),size(tmp1))
        tmp1=v[iF][:,1:end-1]; tmp2=v[iF][:,2:end]
        vC[iF]=reshape(nanmean([tmp1[:] tmp2[:]],2),size(tmp1))
    end

    return uC.*G.AngleCS-vC.*G.AngleSN, uC.*G.AngleSN+vC.*G.AngleCS
end

function UVtoSpeed!(uC::AbstractMeshArray,vC::AbstractMeshArray,G::NamedTuple,dD)
  (u,v)=exch_UV(uC,vC)
  for iF=1:u.grid.nFaces
    for i in 1:size(vC[iF],1)
      for j in 1:size(vC[iF],2)
        u0=nanmean(u[iF][i,j],u[iF][i+1,j])
        v0=nanmean(v[iF][i,j],v[iF][i,j+1])
        u1=u0*G.AngleCS[iF][i,j]-v0*G.AngleSN[iF][i,j]
        v1=u0*G.AngleSN[iF][i,j]+v0*G.AngleCS[iF][i,j]
        dD[iF][i,j]=sqrt(u1^2 + v1^2)
      end
    end
  end
end

"""
    UVtoTransport(U,V,G::NamedTuple)

Convert e.g. velocity (m/s) to transport (m^3/s) by multiplying by `DRF` and `DXG`,`DYG`.
"""
function UVtoTransport(U::AbstractMeshArray,V::AbstractMeshArray,G::NamedTuple)
  uTr=deepcopy(U)
  vTr=deepcopy(V)
  UVtoTransport!(uTr,vTr,G)
  return uTr,vTr
end

"""
    UVtoTransport!(U,V,G::NamedTuple)

Convert e.g. velocity (m/s) to transport (m^3/s) by multiplying by `DRF` and `DXG`,`DYG`.

(in place)
"""
function UVtoTransport!(U::AbstractMeshArray,V::AbstractMeshArray,G::NamedTuple)
  for i in eachindex(U)
      for j in eachindex(U[i])
          !isfinite(U[i][j]) ? U[i][j]=0.0 : nothing
          !isfinite(V[i][j]) ? V[i][j]=0.0 : nothing
          U[i][j]=G.DRF[i[2]]*U[i][j].*G.DYG[i[1]][j]
          V[i][j]=G.DRF[i[2]]*V[i][j].*G.DXG[i[1]][j]
      end
  end
end

"""
compute bolus velocty field (bolusU,bolusV,bolusW) 
	from gm streamfunction (GM_PsiX,GM_PsiY)
"""
function calc_bolus(GM_PsiX,GM_PsiY, Γ)
    nr=length(Γ.RC);
    mskW = 0 .*Γ.hFacW; mskS = 0 .*Γ.hFacS; 
    mskC = 0 .*Γ.hFacC;
    for ff in eachindex(mskC)
        mskW.f[ff][Γ.hFacW.f[ff] .> 0] .= 1
        mskS.f[ff][Γ.hFacS.f[ff] .> 0] .= 1;
        mskC.f[ff][Γ.hFacC.f[ff] .> 0] .= 1;
    end

    GM_PsiX[findall((!isfinite).(GM_PsiX))]=0;
    GM_PsiY[findall((!isfinite).(GM_PsiY))]=0;
    
    bolusU=0*Γ.hFacW;
    bolusV=0*Γ.hFacS;
    for k=1:nr-1;
        bolusU.f[:,k].=(GM_PsiX.f[:,k+1].-GM_PsiX.f[:,k])/Γ.DRF[k];
        bolusV.f[:,k].=(GM_PsiY.f[:,k+1].-GM_PsiY.f[:,k])/Γ.DRF[k];
    end;
    bolusU.f[:, nr] .= 0*GM_PsiX.f[:,nr] .-GM_PsiX.f[:,nr]./Γ.DRF[nr];
    bolusV.f[:, nr] .= 0*GM_PsiY.f[:,nr] .-GM_PsiY.f[:,nr]./Γ.DRF[nr];

    bolusU=bolusU.*mskW;
    bolusV=bolusV.*mskS;
    
    #and its vertical part
    #   (seems correct, leading to 0 divergence)
    tmp_x=GM_PsiX.*Γ.DYG
    tmp_y=GM_PsiY.*Γ.DXG
    tmp_w=0*tmp_x

    for k in 1:nr
        (tmpU,tmpV)=exch_UV(tmp_x[:, k],tmp_y[:, k])
        for a=1:tmpU.grid.nFaces
            (s1,s2)=size(tmp_x.f[a])
            tmpU1=view(tmpU.f[a],1:s1,1:s2)
            tmpU2=view(tmpU.f[a],2:s1+1,1:s2)
            tmpV1=view(tmpV.f[a],1:s1,1:s2)
            tmpV2=view(tmpV.f[a],1:s1,2:s2+1)
            tmp_w.f[a, k] = tmpU2 - tmpU1+tmpV2 -tmpV1
            tmp_w.f[a, k] = tmp_w.f[a, k] ./ Γ.RAC.f[a]
        end
    end

    bolusW=tmp_w.*mskC;

    return bolusU, bolusV, bolusW
end
