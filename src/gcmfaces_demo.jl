
## demo functions:
#
#examples:
# (D,Dexch,Darr,DD)=demo1("LLC90");
# (Rini,Rend,DXCsm,DYCsm)=demo2();
# @time Rend=smooth(Rini,DXCsm,DYCsm);
#
# GCMGridOnes("cs",6,100);
# (Rini,Rend,DXCsm,DYCsm)=demo2();
# @time Rend=smooth(Rini,DXCsm,DYCsm);

function demo1(gridChoice)

GCMGridSpec(gridChoice)

D=read_bin(MeshArrays.grDir*"Depth.data",MeshArrays.ioPrec)

1000+D
D+1000
D+D
D-1000
1000-D
D-D
1000*D
D*1000
D*D
D/1000
1000/D
D/D

Dexch=exchange(D,4)
Darr=convert2array(D)
DD=convert2array(Darr)

GCMGridLoad()

(dFLDdx, dFLDdy)=gradient(MeshArrays.YC)
(dFLDdxEx,dFLDdyEx)=exchange(dFLDdx,dFLDdy,4)

view(MeshArrays.hFacC,:,:,40)
#show(fsize(MeshArrays.hFacC,1))
#show(fsize(view(MeshArrays.hFacC,:,:,40),1))

return (D,Dexch,Darr,DD)

end

##

function demo2()

#Pre-requisite:
#either load predefined grid:
#  demo1("CS32");
#or create a simple one:
#  GCMGridOnes("cs",6,10);

#initialize 2D field of random numbers
tmp1=convert2gcmfaces(MeshArrays.XC);
tmp1=randn(Float32,size(tmp1));
Rini=convert2gcmfaces(tmp1);

#apply land mask
if ndims(MeshArrays.hFacC.f[1])>2
    tmp1=mask(view(MeshArrays.hFacC,:,:,1),NaN,0);
else
    tmp1=mask(MeshArrays.hFacC,NaN,0);
end
msk=fill(1.,tmp1) + 0. *tmp1;
Rini=msk*Rini;

#specify smoothing length scales in x, y directions
DXCsm=3*MeshArrays.DXC; DYCsm=3*MeshArrays.DYC;

#apply smoother
Rend=smooth(Rini,DXCsm,DYCsm);

return (Rini,Rend,DXCsm,DYCsm)

end
