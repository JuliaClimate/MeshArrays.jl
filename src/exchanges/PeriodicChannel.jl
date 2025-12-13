## Grid-specific implementations: PeriodicChannel case

function exch_T_N_ll(fld::AbstractMeshArray,N::Integer)

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
