
## gradient methods

function gradient(inFLD::gcmfaces)
(dFLDdx, dFLDdy)=gradient(inFLD,true)
return dFLDdx, dFLDdy
end

function gradient(inFLD::gcmfaces,doDIV::Bool)

exFLD=exchange(inFLD,1)
dFLDdx=gcmfaces(inFLD.nFaces,inFLD.grTopo)
dFLDdy=gcmfaces(inFLD.nFaces,inFLD.grTopo)

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

function gradient(inFLD::gcmfaces,iDXC::gcmfaces,iDYC::gcmfaces)

exFLD=exchange(inFLD,1)
dFLDdx=gcmfaces(inFLD.nFaces,inFLD.grTopo)
dFLDdy=gcmfaces(inFLD.nFaces,inFLD.grTopo)

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

function mask(fld::gcmfaces)
fldmsk=mask(fld,NaN)
return fldmsk
end

function mask(fld::gcmfaces, val::Number)
  fldmsk=gcmfaces(fld.nFaces,fld.grTopo)
  for a=1:fld.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> !isfinite(x) ? val : x, tmp1 )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

function mask(fld::gcmfaces, val::Number, noval::Number)
  fldmsk=gcmfaces(fld.nFaces,fld.grTopo)
  for a=1:fld.nFaces
    tmp1=copy(fld.f[a])
    replace!(x -> x==noval ? val : x, tmp1  )
    fldmsk.f[a]=tmp1
  end
  return fldmsk
end

## convergence methods

function convergence(uFLD::gcmfaces,vFLD::gcmfaces);

#important note:
#  Normally uFLD, vFLD should not contain any NaN;
#  if otherwise then something this may be needed:
#  uFLD=mask(uFLD,0.0); vFLD=mask(vFLD,0.0);

CONV=gcmfaces(uFLD.nFaces,uFLD.grTopo)

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
iDXC=gcmfaces(FLD.nFaces,FLD.grTopo)
iDYC=gcmfaces(FLD.nFaces,FLD.grTopo)
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
  tmpU=gcmfaces(FLD.nFaces,FLD.grTopo)
  tmpV=gcmfaces(FLD.nFaces,FLD.grTopo)
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

##
