
## gradient methods

"""
    gradient(inFLD::MeshArray,GridVariables::Dict)

Compute spatial derivatives. Other methods:
```
gradient(inFLD::MeshArray,GridVariables::Dict,doDIV::Bool)
gradient(inFLD::MeshArray,iDXC::MeshArray,iDYC::MeshArray)
```
"""
function gradient(inFLD::MeshArray,GridVariables::Dict)
(dFLDdx, dFLDdy)=gradient(inFLD,GridVariables,true)
return dFLDdx, dFLDdy
end

function gradient(inFLD::MeshArray,GridVariables::Dict,doDIV::Bool)

exFLD=exchange(inFLD,1)
dFLDdx=similar(inFLD)
dFLDdy=similar(inFLD)

for a=1:inFLD.grid.nFaces
  (s1,s2)=size(exFLD.f[a])
  tmpA=view(exFLD.f[a],2:s1-1,2:s2-1)
  tmpB=tmpA-view(exFLD.f[a],1:s1-2,2:s2-1)
  tmpC=tmpA-view(exFLD.f[a],2:s1-1,1:s2-2)
  if doDIV
    dFLDdx.f[a]=tmpB./GridVariables["DXC"].f[a]
    dFLDdy.f[a]=tmpC./GridVariables["DYC"].f[a]
  else
    dFLDdx.f[a]=tmpB
    dFLDdy.f[a]=tmpC
  end
end

return dFLDdx, dFLDdy
end

function gradient(inFLD::MeshArray,iDXC::MeshArray,iDYC::MeshArray)

exFLD=exchange(inFLD,1)
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

## mask methods

function mask(fld::MeshArray)
fldmsk=mask(fld,NaN)
return fldmsk
end

"""
    mask(fld::MeshArray, val::Number)

Replace non finite values with val. Other methods:
```
mask(fld::MeshArray)
mask(fld::MeshArray, val::Number, noval::Number)
```
"""
function mask(fld::MeshArray, val::Number)
  fldmsk=similar(fld)
  for a=1:fld.grid.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> !isfinite(x) ? val : x, tmp1 )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

function mask(fld::MeshArray, val::Number, noval::Number)
  fldmsk=similar(fld)
  for a=1:fld.grid.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> x==noval ? val : x, tmp1  )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

## convergence methods

"""
    convergence(uFLD::MeshArray,vFLD::MeshArray)

Compute convergence of a vector field
"""
function convergence(uFLD::MeshArray,vFLD::MeshArray)

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
    smooth(FLD::MeshArray,DXCsm::MeshArray,DYCsm::MeshArray,GridVariables::Dict)

Smooth out scales below DXCsm / DYCsm via diffusion
"""
function smooth(FLD::MeshArray,DXCsm::MeshArray,DYCsm::MeshArray,GridVariables::Dict)

#important note:
#input FLD should be land masked (NaN/1) by caller if needed

#get land masks (NaN/1):
mskC=fill(1.0,FLD) + 0.0 * mask(FLD)
(mskW,mskS)=gradient(FLD,GridVariables,false)
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
  iDXC.f[a]=1.0./GridVariables["DXC"].f[a]
  iDYC.f[a]=1.0./GridVariables["DYC"].f[a]
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

#diffusion operator times DYG / DXG:
KuxFac=mskW*DXCsm*DXCsm/T/2.0*GridVariables["DYG"];
KvyFac=mskS*DYCsm*DYCsm/T/2.0*GridVariables["DXG"];

#time steping factor:
dtFac=dt*mskC/GridVariables["RAC"];

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
    ThroughFlow(VectorField,IntegralPath,GridVariables::Dict)

Compute transport through an integration path
"""
function ThroughFlow(VectorField,IntegralPath,GridVariables::Dict)

    #Note: vertical intergration is not always wanted; left for user to do outside

    U=VectorField["U"]
    V=VectorField["V"]

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
        #mskW=IntegralPath["mskW"]
        #do_dxory==1 ? mskW=mskW*GridVariables["DYG"] : nothing
        #do_dz==1 ? mskW=GridVariables["DRF"][i3]*mskW : nothing
        #mskS=IntegralPath["mskS"]
        #do_dxory==1 ? mskS=mskS*GridVariables["DXG"] : nothing
        #do_dz==1 ? mskS=GridVariables["DRF"][i3]*mskS : nothing
        #
        #method 2: less slow
        tabW=IntegralPath["tabW"]
        tabS=IntegralPath["tabS"]
        for i4=1:n[4]
            #method 1: quite slow
            #trsp[1,i3,i4]=sum(mskW*U[:,:,i3,i4])+sum(mskS*V[:,:,i3,i4])
            #
            #method 2: less slow
            trsp[1,i3,i4]=0.0
            for k=1:size(tabW,1)
                (a,i1,i2,w)=tabW[k,:]
                do_dxory==1 ? w=w*GridVariables["DYG"].f[a][i1,i2] : nothing
                do_dz==1 ? w=w*GridVariables["DRF"][i3] : nothing
                isdefined(U,:fIndex) ? u=U.f[a,i3,i4][i1,i2] : u=U.f[a][i1,i2,i3,i4]
                trsp[1,i3,i4]=trsp[1,i3,i4]+w*u
            end
            for k=1:size(tabS,1)
                (a,i1,i2,w)=tabS[k,:]
                do_dxory==1 ? w=w*GridVariables["DXG"].f[a][i1,i2] : nothing
                do_dz==1 ? w=w*GridVariables["DRF"][i3] : nothing
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
    LatitudeCircles(LatValues,GridVariables::Dict)

Compute integration paths that follow latitude circles
"""
function LatitudeCircles(LatValues,GridVariables::Dict)

    LatitudeCircles=Array{Dict}(undef,length(LatValues))

    for j=1:length(LatValues)
        mskCint=1*(GridVariables["YC"] .>= LatValues[j])
        mskC=similar(mskCint)
        mskW=similar(mskCint)
        mskS=similar(mskCint)

        mskCint=exchange(mskCint,1)

        for i=1:mskCint.grid.nFaces
            tmp1=mskCint.f[i]
            # tracer mask:
            tmp2=tmp1[2:end-1,1:end-2]+tmp1[2:end-1,3:end]+
            tmp1[1:end-2,2:end-1]+tmp1[3:end,2:end-1]
            mskC.f[i]=1((tmp2.>0).&(tmp1[2:end-1,2:end-1].==0))
            # velocity masks:
            mskW.f[i]=tmp1[2:end-1,2:end-1] - tmp1[1:end-2,2:end-1]
            mskS.f[i]=tmp1[2:end-1,2:end-1] - tmp1[2:end-1,1:end-2]
        end

        function MskToTab(msk::MeshArray)
          n=Int(sum(msk .!= 0)); k=0
          tab=Array{Int,2}(undef,n,4)
          for i=1:msk.grid.nFaces
            a=msk.f[i]
            b=findall( a .!= 0)
            for ii in eachindex(b)
              k += 1
              tab[k,:]=[i,b[ii][1],b[ii][2],a[b[ii]]]
            end
          end
          return tab
        end

        LatitudeCircles[j]=Dict("lat"=>LatValues[j],
        "tabC"=>MskToTab(mskC),"tabW"=>MskToTab(mskW),"tabS"=>MskToTab(mskS))
    end

    return LatitudeCircles

end

##
