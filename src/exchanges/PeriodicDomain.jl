
## Grid-specific implementations: PeriodicDomain case

function exch_T_N_dpdo(fld::AbstractMeshArray,N::Integer)

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

function exch_UV_N_dpdo(fldU,fldV,N)
  FLDU=exch_T_N_dpdo(fldU,N)
  FLDV=exch_T_N_dpdo(fldV,N)
  return FLDU,FLDV
end

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
