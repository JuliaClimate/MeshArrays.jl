
## This file contains the exchange and exch_UV functions
# along with grid-specific methods (exch_T_N.jl, etc.)

## User Front Ends

"""
    exchange(fld::MeshArray)

Exchange / transfer data between neighboring arrays. Other methods are

    exchange(fld::MeshArray,N::Integer)
    exchange(u::MeshArray,v::MeshArray)
    exchange(u::MeshArray,v::MeshArray,N::Integer)
"""
function exchange(fld::MeshArray)
  MeshArray_wh(exch_T_N(fld,1),1)
end

function exchange(fld::MeshArray,N::Integer)
  MeshArray_wh(exch_T_N(fld,N),N)
end

function exchange(u::MeshArray,v::MeshArray)
  MeshArray_wh.(exch_UV_N(u,v,1),1)
end

function exchange(u::MeshArray,v::MeshArray,N::Integer)
  MeshArray_wh.(exch_UV_N(u,v,N),N)
end

## dispatch over grid types

#note: the "CubeSphere" implementation covers both cs and llc

function exch_T_N(fld,N)

if fld.grid.class=="LatLonCap"
  FLD=exch_T_N_cs(fld,N)
elseif fld.grid.class=="CubeSphere"
  FLD=exch_T_N_cs(fld,N)
elseif fld.grid.class=="PeriodicChannel"
  FLD=exch_T_N_ll(fld,N)
elseif fld.grid.class=="PeriodicDomain"
  FLD=exch_T_N_dpdo(fld,N)
else
  error("unknown grid.class case")
end

return FLD

end

function exch_UV_N(u,v,N)

if u.grid.class=="LatLonCap"
  (uex,vex)=exch_UV_N_cs(u,v,N)
elseif u.grid.class=="CubeSphere"
  (uex,vex)=exch_UV_N_cs(u,v,N)
elseif u.grid.class=="PeriodicChannel"
  (uex,vex)=exch_UV_N_ll(u,v,N)
elseif u.grid.class=="PeriodicDomain"
  (uex,vex)=exch_UV_N_dpdo(u,v,N)
else
  error("unknown grid.class case")
end

return uex,vex

end

function exch_UV(u,v)

if u.grid.class=="LatLonCap"
  (uex,vex)=exch_UV_cs(u,v)
elseif u.grid.class=="CubeSphere"
  (uex,vex)=exch_UV_cs(u,v)
elseif u.grid.class=="PeriodicChannel"
  (uex,vex)=exch_UV_ll(u,v)
elseif u.grid.class=="PeriodicDomain"
  (uex,vex)=exch_UV_dpdo(u,v)
else
  error("unknown grid.class case")
end

return uex,vex

end

## Grid-specific implementations: PeriodicDomain case

function exch_T_N_dpdo(fld::MeshArray,N::Integer)

fillval=0.0

ni,nj=Int.(fld.grid.ioSize[:]./fld.grid.fSize[1][:])
s=fld.fSize
FLD=similar(fld;m=fld.meta);

for i=1:ni
  for j=1:nj
    k=i+ni*(j-1)
    kS=j-1; kS==0 ? kS=nj : nothing; kS=i+ni*(kS-1)
    kN=j+1; kN==nj+1 ? kN=1 : nothing; kN=i+ni*(kN-1)
    kW=i-1; kW==0 ? kW=ni : nothing; kW=kW+ni*(j-1)
    kE=i+1; kE==ni+1 ? kE=1 : nothing; kE=kE+ni*(j-1)

    #step 1

    FLD.f[k]=fill(fillval,s[k].+2N);
    @views FLD.f[k][N+1:N+s[k][1],N+1:N+s[k][2]]=fld.f[k];

    #step 2

    iW=(s[k][1]-N+1:s[k][1],1:s[k][2]);
    iE=(1:N,1:s[k][2]);
    jW=(1:N,N+1:N+s[k][2]);
    jE=(N+1+s[k][1]:2N+s[k][1],N+1:N+s[k][2]);
    FLD.f[k][jW[1],jW[2]]=view(fld.f[kW],iW[1],iW[2])
    FLD.f[k][jE[1],jE[2]]=view(fld.f[kE],iE[1],iE[2])

    #step 3

    iS=(1:s[k][1],s[k][2]-N+1:s[k][2]);
    iN=(1:s[k][1],1:N);
    jS=(N+1:N+s[k][1],1:N);
    jN=(N+1:N+s[k][1],N+1+s[k][2]:2N+s[k][2]);
    FLD.f[k][jS[1],jS[2]]=view(fld.f[kS],iS[1],iS[2])
    FLD.f[k][jN[1],jN[2]]=view(fld.f[kN],iN[1],iN[2])

  end
end

return FLD

end

##

function exch_UV_N_dpdo(fldU,fldV,N);

FLDU=exch_T_N_dpdo(fldU,N);
FLDV=exch_T_N_dpdo(fldV,N);

return FLDU,FLDV

end

##

function exch_UV_dpdo(fldU,fldV);

fillval=0.0

ni,nj=Int.(fldU.grid.ioSize[:]./fldU.grid.fSize[1][:])
s=fldU.fSize
FLDU=similar(fldU;m=fldU.meta)
FLDV=similar(fldV;m=fldV.meta)

for i=1:ni
  for j=1:nj
    k=i+ni*(j-1)
    kS=j-1; kS==0 ? kS=nj : nothing; kS=i+ni*(kS-1)
    kN=j+1; kN==nj+1 ? kN=1 : nothing; kN=i+ni*(kN-1)
    kW=i-1; kW==0 ? kW=ni : nothing; kW=kW+ni*(j-1)
    kE=i+1; kE==ni+1 ? kE=1 : nothing; kE=kE+ni*(j-1)

    #step 1

    FLDU.f[k]=fill(fillval,s[k][1]+1,s[k][2]);
    FLDV.f[k]=fill(fillval,s[k][1],s[k][2]+1);
    @views FLDU.f[k][1:s[k][1],1:s[k][2]]=fldU.f[k];
    @views FLDV.f[k][1:s[k][1],1:s[k][2]]=fldV.f[k];

    #step 2

    FLDU.f[k][s[k][1]+1,1:s[k][2]]=view(fldU.f[kE],1,1:s[k][2])

    #step 3

    FLDV.f[k][1:s[k][1],s[k][2]+1]=view(fldV.f[kN],1:s[k][1],1)

  end
end

return FLDU,FLDV

end

## Grid-specific implementations: PeriodicChannel case

function exch_T_N_ll(fld::MeshArray,N::Integer)

fillval=0.0

#step 1

s=size.(fld.f);
FLD=similar(fld;m=fld.meta);
FLD.f[1]=fill(fillval,s[1].+2N);
@views FLD.f[1][N+1:N+s[1][1],N+1:N+s[1][2]]=fld.f[1];

#step 2

iW=(s[1][1]-N+1:s[1][1],1:s[1][2]);
iE=(1:N,1:s[1][2]);
jW=(1:N,N+1:N+s[1][2]);
jE=(N+1+s[1][1]:2N+s[1][1],N+1:N+s[1][2]);
FLD.f[1][jW[1],jW[2]]=view(fld.f[1],iW[1],iW[2])
FLD.f[1][jE[1],jE[2]]=view(fld.f[1],iE[1],iE[2])

return FLD

end

##

function exch_UV_N_ll(fldU,fldV,N);

FLDU=exch_T_N_ll(fldU,N);
FLDV=exch_T_N_ll(fldV,N);

return FLDU,FLDV

end

##

function exch_UV_ll(fldU,fldV);

fillval=0.0

#step 1

s=size.(fldU.f);
FLDU=similar(fldU;m=fldU.meta);
FLDV=similar(fldV;m=fldV.meta);

FLDU.f[1]=fill(fillval,s[1][1]+1,s[1][2]);
FLDV.f[1]=fill(fillval,s[1][1],s[1][2]+1);
@views FLDU.f[1][1:s[1][1],1:s[1][2]]=fldU.f[1];
@views FLDV.f[1][1:s[1][1],1:s[1][2]]=fldV.f[1];

#step 2

FLDU.f[1][s[1][1]+1,1:s[1][2]]=view(fldU.f[1],1,1:s[1][2])

return FLDU,FLDV

end

## Grid-specific implementations: CubeSphere & LatLonCap case

#note: the "CubeSphere" implementation covers both cs and llc

function exch_T_N_cs(fld::MeshArray,N::Integer)

fillval=0.0

#step 1

s=size.(fld.f)
nf=fld.grid.nFaces
nf==5 ? s=vcat(s,s[3]) : nothing
tp=fld.grid.class
FLD=similar(fld;m=fld.meta)

for i=1:nf; FLD.f[i]=fill(fillval,s[i].+2N); end;
#code below yields strange, seemingly incorrect results:
#for i=1:nf; FLD.f[i]=Array{eltype(fld.f[i])}(undef,s[i].+2N); end;

#all versions below yield same @time and memory (despite diff in allocs)
for i=1:nf;
# FLD.f[i][N+1:end-N,N+1:end-N]=fld.f[i];
 @views FLD.f[i][N+1:N+s[i][1],N+1:N+s[i][2]]=fld.f[i];
end;

#step 2

(ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN)=exch_cs_viewfunctions();

for a=1:nf
(jW, jE, jS, jN)=exch_cs_target(s[a],N)
(aW,aE,aS,aN,iW,iE,iS,iN)=exch_cs_sources(a,s,N)
if !iseven(a)
 aW <= nf ? FLD.f[a][jW[1],jW[2]]=ovfW(fld.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLD.f[a][jE[1],jE[2]]=ovfE(fld.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLD.f[a][jS[1],jS[2]]=ovfS(fld.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLD.f[a][jN[1],jN[2]]=ovfN(fld.f[aN],iN[1],iN[2]) : nothing
else
 aW <= nf ? FLD.f[a][jW[1],jW[2]]=evfW(fld.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLD.f[a][jE[1],jE[2]]=evfE(fld.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLD.f[a][jS[1],jS[2]]=evfS(fld.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLD.f[a][jN[1],jN[2]]=evfN(fld.f[aN],iN[1],iN[2]) : nothing
end
end

return FLD

end

##

function exch_UV_N_cs(fldU::MeshArray,fldV::MeshArray,N::Integer)

fillval=0.0

#step 1

s=size.(fldU.f)
nf=fldU.grid.nFaces
nf==5 ? s=vcat(s,s[3]) : nothing
tp=fldU.grid.class
FLDU=similar(fldU;m=fldU.meta)
FLDV=similar(fldV;m=fldV.meta)

for i=1:nf;
 FLDU.f[i]=fill(fillval,s[i].+2N);
 FLDV.f[i]=fill(fillval,s[i].+2N);
 @views FLDU.f[i][N+1:N+s[i][1],N+1:N+s[i][2]]=fldU.f[i];
 @views FLDV.f[i][N+1:N+s[i][1],N+1:N+s[i][2]]=fldV.f[i];
end;

#step 2

(ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN)=exch_cs_viewfunctions();

for a=1:nf
(jW, jE, jS, jN)=exch_cs_target(s[a],N)
(aW,aE,aS,aN,iW,iE,iS,iN)=exch_cs_sources(a,s,N)
if !iseven(a)
 aW <= nf ? FLDU.f[a][jW[1],jW[2]]=ovfW(fldV.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLDU.f[a][jE[1],jE[2]]=ovfE(fldU.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLDU.f[a][jS[1],jS[2]]=ovfS(fldU.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLDU.f[a][jN[1].+1,jN[2]]=-ovfN(fldV.f[aN],iN[1],iN[2]) : nothing
 aW <= nf ? FLDV.f[a][jW[1],jW[2].+1]=-ovfW(fldU.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLDV.f[a][jE[1],jE[2]]=ovfE(fldV.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLDV.f[a][jS[1],jS[2]]=ovfS(fldV.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLDV.f[a][jN[1],jN[2]]=ovfN(fldU.f[aN],iN[1],iN[2]) : nothing
else
 aW <= nf ? FLDU.f[a][jW[1],jW[2]]=evfW(fldU.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLDU.f[a][jE[1],jE[2]]=evfE(fldV.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLDU.f[a][jS[1].+1,jS[2]]=-evfS(fldV.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLDU.f[a][jN[1],jN[2]]=evfN(fldU.f[aN],iN[1],iN[2]) : nothing
 aW <= nf ? FLDV.f[a][jW[1],jW[2]]=evfW(fldV.f[aW],iW[1],iW[2]) : nothing
 aE <= nf ? FLDV.f[a][jE[1],jE[2].+1]=-evfE(fldU.f[aE],iE[1],iE[2]) : nothing
 aS <= nf ? FLDV.f[a][jS[1],jS[2]]=evfS(fldU.f[aS],iS[1],iS[2]) : nothing
 aN <= nf ? FLDV.f[a][jN[1],jN[2]]=evfN(fldV.f[aN],iN[1],iN[2]) : nothing
end
end

return FLDU,FLDV

end

##

function exch_UV_cs(fldU::MeshArray,fldV::MeshArray)

fillval=0.0

#step 1

s=size.(fldU.f)
nf=fldU.grid.nFaces
nf==5 ? s=vcat(s,s[3]) : nothing
tp=fldU.grid.class
FLDU=similar(fldU;m=fldU.meta)
FLDV=similar(fldV;m=fldV.meta)

for i=1:nf
  FLDU.f[i]=fill(fillval,s[i][1]+1,s[i][2]);
  FLDV.f[i]=fill(fillval,s[i][1],s[i][2]+1);
  @views FLDU.f[i][1:s[i][1],1:s[i][2]]=fldU.f[i];
  @views FLDV.f[i][1:s[i][1],1:s[i][2]]=fldV.f[i];
end

 #step 2

(ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN)=exch_cs_viewfunctions();

for a=1:nf
(jW, jE, jS, jN)=exch_cs_target(s[a],1)
(aW,aE,aS,aN,iW,iE,iS,iN)=exch_cs_sources(a,s,1)
if !iseven(a)
 aE <= nf ? FLDU.f[a][jE[1].-1,jE[2].-1]=ovfE(fldU.f[aE],iE[1],iE[2]) : nothing
 aN <= nf ? FLDV.f[a][jN[1].-1,jN[2].-1]=ovfN(fldU.f[aN],iN[1],iN[2]) : nothing
else
 aE <= nf ? FLDU.f[a][jE[1].-1,jE[2].-1]=evfE(fldV.f[aE],iE[1],iE[2]) : nothing
 aN <= nf ? FLDV.f[a][jN[1].-1,jN[2].-1]=evfN(fldV.f[aN],iN[1],iN[2]) : nothing
end
end

return FLDU,FLDV

end

## Convenience functions used in the cs & llc case

function exch_cs_target(sa::Tuple{Int64,Int64},N::Integer)

    #target array indices
    jW=(1:N,N+1:N+sa[2]);
    jE=(N+1+sa[1]:2N+sa[1],N+1:N+sa[2]);
    jS=(N+1:N+sa[1],1:N);
    jN=(N+1:N+sa[1],N+1+sa[2]:2N+sa[2]);

    return jW, jE, jS, jN

end

function exch_cs_sources(a::Integer,s::Array{Tuple{Int64,Int64},1},N::Integer)

#source array IDs
aW=0; aE=0; aS=0; aN=0;
if a==1;     aW=5; aE=2; aS=6; aN=3;
elseif a==2; aW=1; aE=4; aS=6; aN=3;
elseif a==3; aW=1; aE=4; aS=2; aN=5;
elseif a==4; aW=3; aE=6; aS=2; aN=5;
elseif a==5; aW=3; aE=6; aS=4; aN=1;
elseif a==6; aW=5; aE=2; aS=4; aN=1;
else; error("Array index is out of bounds.");
end;

if !iseven(a)
    #source array indices
    iW=(1:s[aW][1],s[aW][2]-N+1:s[aW][2]);
    iE=(1:N,1:s[aE][2]);
    iS=(1:s[aS][1],s[aS][2]-N+1:s[aS][2]);
    iN=(1:N,1:s[aN][2]);
else
    #source array indices
    iW=(s[aW][1]-N+1:s[aW][1],1:s[aW][2]);
    iE=(1:s[aE][1],1:N);
    iS=(s[aS][1]-N+1:s[aS][1],1:s[aS][2]);
    iN=(1:s[aN][1],1:N);
end

return aW,aE,aS,aN,iW,iE,iS,iN
end

function exch_cs_viewfunctions()
#view functions for odd numbered arrays
ovfW(x,i,j)=PermutedDimsArray(view(x,reverse(i),j),(2,1))
ovfE(x,i,j)=view(x,i,j)
ovfS(x,i,j)=view(x,i,j)
ovfN(x,i,j)=PermutedDimsArray(view(x,i,reverse(j)),(2,1))
#view functions for even numbered arrays
evfW(x,i,j)=view(x,i,j)
evfE(x,i,j)=PermutedDimsArray(view(x,reverse(i),j),(2,1))
evfS(x,i,j)=PermutedDimsArray(view(x,i,reverse(j)),(2,1))
evfN(x,i,j)=view(x,i,j)
#
return ovfW,ovfE,ovfS,ovfN,evfW,evfE,evfS,evfN
end

##
