
## gradient methods

"""
    gradient(inFLD::gcmfaces)

Compute spatial derivatives.
"""
function gradient(inFLD::gcmfaces)
(dFLDdx, dFLDdy)=gradient(inFLD,true)
return dFLDdx, dFLDdy
end

"""
    gradient(inFLD::gcmfaces,doDIV::Bool)

Compute spatial derivatives with or without dividing by grid scale.
"""
function gradient(inFLD::gcmfaces,doDIV::Bool)

exFLD=exchange(inFLD,1)
dFLDdx=similar(inFLD)
dFLDdy=similar(inFLD)

for a=1:inFLD.nFaces;
  (s1,s2)=fsize(exFLD,a)
  tmpA=view(exFLD.f[a],2:s1-1,2:s2-1)
  tmpB=tmpA-view(exFLD.f[a],1:s1-2,2:s2-1)
  tmpC=tmpA-view(exFLD.f[a],2:s1-1,1:s2-2)
  if doDIV
    dFLDdx.f[a]=tmpB./MeshArrays.DXC.f[a]
    dFLDdy.f[a]=tmpC./MeshArrays.DYC.f[a]
  else
    dFLDdx.f[a]=tmpB
    dFLDdy.f[a]=tmpC
  end
end

return dFLDdx, dFLDdy
end

"""
    gradient(inFLD::gcmfaces,iDXC::gcmfaces,iDYC::gcmfaces)

Compute spatial derivatives and multiply by inverse grid scale.
"""
function gradient(inFLD::gcmfaces,iDXC::gcmfaces,iDYC::gcmfaces)

exFLD=exchange(inFLD,1)
dFLDdx=similar(inFLD)
dFLDdy=similar(inFLD)

for a=1:inFLD.nFaces;
  (s1,s2)=fsize(exFLD,a)
  tmpA=view(exFLD.f[a],2:s1-1,2:s2-1)
  tmpB=tmpA-view(exFLD.f[a],1:s1-2,2:s2-1)
  tmpC=tmpA-view(exFLD.f[a],2:s1-1,1:s2-2)
  dFLDdx.f[a]=tmpB.*iDXC.f[a]
  dFLDdy.f[a]=tmpC.*iDYC.f[a]
end

return dFLDdx, dFLDdy
end

## mask methods

"""
    mask(fld::gcmfaces)

Call mask(fld,NaN)
"""
function mask(fld::gcmfaces)
fldmsk=mask(fld,NaN)
return fldmsk
end

"""
    mask(fld::gcmfaces, val::Number)

Replace non finite values with val
"""
function mask(fld::gcmfaces, val::Number)
  fldmsk=similar(fld)
  for a=1:fld.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> !isfinite(x) ? val : x, tmp1 )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

"""
    mask(fld::gcmfaces, val::Number)

Replace noval instances with val
"""
function mask(fld::gcmfaces, val::Number, noval::Number)
  fldmsk=similar(fld)
  for a=1:fld.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> x==noval ? val : x, tmp1  )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

## convergence methods

"""
    convergence(uFLD::gcmfaces,vFLD::gcmfaces)

Compute convergence of a vector field
"""
function convergence(uFLD::gcmfaces,vFLD::gcmfaces)

#important note:
#  Normally uFLD, vFLD should not contain any NaN;
#  if otherwise then something this may be needed:
#  uFLD=mask(uFLD,0.0); vFLD=mask(vFLD,0.0);

CONV=similar(uFLD)

(tmpU,tmpV)=exch_UV(uFLD,vFLD)
for a=1:tmpU.nFaces
  (s1,s2)=fsize(uFLD,a)
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
    smooth(FLD::gcmfaces,DXCsm::gcmfaces,DYCsm::gcmfaces)

Smooth out scales below threshold
"""
function smooth(FLD::gcmfaces,DXCsm::gcmfaces,DYCsm::gcmfaces)

#important note:
#input FLD should be land masked (NaN/1) by caller if needed

#get land masks (NaN/1):
mskC=fill(1.0,FLD) + 0.0 * mask(FLD)
(mskW,mskS)=gradient(FLD,false)
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
for a=1:FLD.nFaces;
  iDXC.f[a]=1.0./MeshArrays.DXC.f[a]
  iDYC.f[a]=1.0./MeshArrays.DYC.f[a]
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
KuxFac=mskW*DXCsm*DXCsm/T/2.0*MeshArrays.DYG;
KvyFac=mskS*DYCsm*DYCsm/T/2.0*MeshArrays.DXG;

#time steping factor:
dtFac=dt*mskC/MeshArrays.RAC;

#loop:
for it=1:nbt
  (dTdxAtU,dTdyAtV)=gradient(FLD,iDXC,iDYC);
  tmpU=similar(FLD)
  tmpV=similar(FLD)
  for a=1:FLD.nFaces
      tmpU.f[a]=dTdxAtU.f[a].*KuxFac.f[a];
      tmpV.f[a]=dTdyAtV.f[a].*KvyFac.f[a];
  end
  tmpC=convergence(tmpU,tmpV);
  for a=1:FLD.nFaces
      FLD.f[a]=FLD.f[a]-dtFac.f[a].*tmpC.f[a];
  end
end

#Apply land mask (NaN/1) to end result:
mskC=mask(mskC,NaN,0.0)
FLD=mskC*FLD

return FLD

end

## TransportThrough function

"""
    TransportThrough(VectorField,IntegralPath)

Computes transport through in integration path
"""
function TransportThrough(VectorField,IntegralPath)

    #Note: vertical intergration is not always wanted; left for user to do outside

    U=VectorField["U"]
    V=VectorField["V"]

    nd=ndims(U)
    #println("nd=$nd and d=$d")

    n=fill(1,4)
    tmp=size(U)
    n[1:nd].=tmp[1:nd]

    haskey(VectorField,"factors") ? f=VectorField["factors"] : f=Array{String,1}(undef,nd)
    haskey(VectorField,"dimensions") ? d=VectorField["dimensions"] : d=Array{String,1}(undef,nd)

    length(d)!=nd ? error("inconsistent specification of dims") : nothing

    trsp=Array{Float64}(undef,1,n[3],n[4])
    do_dz=sum(f.=="dz")
    do_dxory=sum(f.=="dxory")

    for i3=1:n[3]
        #method 1: quite slow
        #mskW=IntegralPath["mskW"]
        #do_dxory==1 ? mskW=mskW*MeshArrays.DYG : nothing
        #do_dz==1 ? mskW=MeshArrays.DRF[i3]*mskW : nothing
        #mskS=IntegralPath["mskS"]
        #do_dxory==1 ? mskS=mskS*MeshArrays.DXG : nothing
        #do_dz==1 ? mskS=MeshArrays.DRF[i3]*mskS : nothing
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
                do_dxory==1 ? w=w*MeshArrays.DYG.f[a][i1,i2] : nothing
                do_dz==1 ? w=w*MeshArrays.DRF[i3] : nothing
                trsp[1,i3,i4]=trsp[1,i3,i4]+w*U.f[a][i1,i2,i3,i4]
            end
            for k=1:size(tabS,1)
                (a,i1,i2,w)=tabS[k,:]
                do_dxory==1 ? w=w*MeshArrays.DXG.f[a][i1,i2] : nothing
                do_dz==1 ? w=w*MeshArrays.DRF[i3] : nothing
                trsp[1,i3,i4]=trsp[1,i3,i4]+w*V.f[a][i1,i2,i3,i4]
            end
        end
    end

    return trsp
end

## LatCircles function

"""
    LatCircles(LatValues)

Compute integration pathes along latitude circles
"""
function LatCircles(LatValues)

    LatCircles=Array{Dict}(undef,length(LatValues))

    for j=1:length(LatValues)
        mskCint=1*(MeshArrays.YC .>= LatValues[j])
        mskC=similar(mskCint)
        mskW=similar(mskCint)
        mskS=similar(mskCint)

        mskCint=exchange(mskCint,1)

        for i=1:mskCint.nFaces
            tmp1=mskCint.f[i]
            # tracer mask:
            tmp2=tmp1[2:end-1,1:end-2]+tmp1[2:end-1,3:end]+
            tmp1[1:end-2,2:end-1]+tmp1[3:end,2:end-1]
            mskC.f[i]=1((tmp2.>0).&(tmp1[2:end-1,2:end-1].==0))
            # velocity masks:
            mskW.f[i]=tmp1[2:end-1,2:end-1] - tmp1[1:end-2,2:end-1]
            mskS.f[i]=tmp1[2:end-1,2:end-1] - tmp1[2:end-1,1:end-2]
        end

        tmp=vec(collect(mskC[:,1]))
        ind = findall(x -> x!=0, tmp)
        tabC=Array{Int,2}(undef,length(ind),4)
        for j=1:length(ind)
            tabC[j,1:3]=collect(fijind(mskC,ind[j]))
            tabC[j,4]=tmp[ind[j]]
        end

        tmp=vec(collect(mskW[:,1]))
        ind = findall(x -> x!=0, tmp)
        tabW=Array{Int,2}(undef,length(ind),5)
        for j=1:length(ind)
            tabW[j,1:3]=collect(fijind(mskW,ind[j]))
            tabW[j,4]=tmp[ind[j]]
            tabW[j,5]=ind[j]
        end

        tmp=vec(collect(mskS[:,1]))
        ind = findall(x -> x!=0, tmp)
        tabS=Array{Int,2}(undef,length(ind),5)
        for j=1:length(ind)
            tabS[j,1:3]=collect(fijind(mskS,ind[j]))
            tabS[j,4]=tmp[ind[j]]
            tabS[j,5]=ind[j]
        end

        LatCircles[j]=Dict("lat"=>LatValues[j],
        #"mskC"=>mskC,"mskW"=>mskW,"mskS"=>mskS,
        "tabC"=>tabC,"tabW"=>tabW,"tabS"=>tabS)
    end

    return LatCircles

end

##
