
function gradient(fld)

  FLD=exchange(fld,1);

  dFLDdx=gcmfaces();
  dFLDdy=gcmfaces();
  for iFace=1:FLD.nFaces;
     tmpA=FLD.f[iFace][2:end-1,2:end-1];
     tmpB=FLD.f[iFace][1:end-2,2:end-1];
#		 tmpC=GCMFaces.DXC.f[iFace];
#		 dFLDdx.f[iFace]=(tmpA-tmpB)./tmpC;
     dFLDdx.f[iFace]=(tmpA-tmpB)./GCMFaces.DXC.f[iFace];
     tmpA=FLD.f[iFace][2:end-1,2:end-1];
     tmpB=FLD.f[iFace][2:end-1,1:end-2];
     dFLDdy.f[iFace]=(tmpA-tmpB)./GCMFaces.DYC.f[iFace];
  end;

  return dFLDdx, dFLDdy

end

##

function mask(fld::gcmfaces)
  fldmsk=mask(fld,(NaN,Inf,0.),NaN)
  return fldmsk
end

function mask(fld::gcmfaces, val::Number)
  fldmsk=mask(fld,(NaN,Inf,0.),val)
  return fldmsk
end

function mask(fld::gcmfaces, cond::Number, val::Number)
  fldmsk=mask(fld,(cond,),val)
  return fldmsk
end

function mask(fld::gcmfaces,cond::Tuple,val::Number)

  fldmsk=gcmfaces(fld.nFaces,fld.grTopo);

  for iFace=1:fld.nFaces;
     tmp1=fld.f[iFace]
     for i=1:length(cond)
       isnan(cond[i]) ? tmp1[findall(isnan.(tmp1))] .= val : nothing
       isinf(cond[i]) ? tmp1[findall(isinf.(tmp1))] .= val : nothing
       isfinite(cond[i]) ? tmp1[findall(tmp1 .== cond[i])] .= val : nothing
     end
     fldmsk.f[iFace]=tmp1
  end

  return fldmsk

end

##

function smooth(fld::gcmfaces,DXCsm::gcmfaces,DYCsm::gcmfaces)

#get land mask:
msk=fill(1.,fld) + 0. * mask(fld,NaN);

#Before scaling the diffusive operator ...
tmp0=DXCsm/GCMFaces.DXC;
tmp0=mask(msk*tmp0,0.);
tmp00=maximum(tmp0);
tmp0=DYCsm/GCMFaces.DYC;
tmp0=mask(msk*tmp0,0.);
tmp00=max(tmp00,maximum(tmp0));

#... determine a suitable time period:
nbt=ceil(1.1*2*tmp00^2);
dt=1.;
T=nbt*dt;
println("nbt in smooth="*"$nbt")

#diffusion operator:
Kux=DXCsm*DXCsm/T/2;
Kvy=DYCsm*DYCsm/T/2;

#loop:
for it=1:nbt
  (dTdxAtU,dTdyAtV)=gradient(fld);
  tmpU=dTdxAtU*Kux*GCMFaces.DYG;
  tmpV=dTdyAtV*Kvy*GCMFaces.DXG;
  tmpC=convergence(tmpU,tmpV);
  fld=fld-dt*msk*tmpC/GCMFaces.RAC;
end

return fld

end

##

function convergence(fldU::gcmfaces,fldV::gcmfaces);

  (tmpU,tmpV)=exch_UV(fldU,fldV);
  tmpU=mask(tmpU,NaN,0);
  tmpV=mask(tmpV,NaN,0);
  tmpC=gcmfaces(tmpU.nFaces,tmpU.grTopo);
  for iFace=1:tmpU.nFaces;
    tmp1=tmpU.f[iFace][1:end-1,:,:,:]-tmpU.f[iFace][2:end,:,:,:];
    tmp2=tmpV.f[iFace][:,1:end-1,:,:]-tmpV.f[iFace][:,2:end,:,:];
    tmpC.f[iFace]=dropdims(tmp1+tmp2,dims=(3,4));
  end;
  return tmpC

end
