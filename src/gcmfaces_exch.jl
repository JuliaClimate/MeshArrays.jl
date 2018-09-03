
## This file contains the following methods:
# exchange.jl
# exch_T_N.jl
# exch_T_N_llc.jl
# exch_T_N_cs.jl
# exch_T_N_ll.jl

#Correspondance to Matlab sym_g are as follows:
#5	rotl90
#6	rot180
#7	rotr90

#an alternative to (or overload of?) rotr90, etc. is needed for mutldim case

## user front end

function exchange(fld::gcmfaces)
  FLD=exch_T_N(fld,1);
end

function exchange(fld::gcmfaces,N::Integer)
  FLD=exch_T_N(fld,N);
end

function exchange(u::gcmfaces,v::gcmfaces)
  (uex,vex)=exch_UV_N(u,v,1);
end

function exchange(u::gcmfaces,v::gcmfaces,N::Integer)
  (uex,vex)=exch_UV_N(u,v,N);
end

## dispatch over grid types

function exch_T_N(fld,N)

if fld.grTopo=="llc";
  FLD=exch_T_N_llc(fld,N);
elseif fld.grTopo=="cs";
  FLD=exch_T_N_cs(fld,N);
elseif fld.grTopo=="ll";
  FLD=exch_T_N_ll(fld,N);
else;
  error("unknown grTopo case");
end;

FLD

end

function exch_UV(u,v)

if u.grTopo=="llc";
  (uex,vex)=exch_UV_llc(u,v);
elseif u.grTopo=="cs";
  (uex,vex)=exch_UV_cs(u,v);
elseif u.grTopo=="ll";
  (uex,vex)=exch_UV_ll(u,v);
else;
  error("unknown grTopo case");
end;

return uex,vex

end

function exch_UV_N(u,v,N)

if u.grTopo=="llc";
  (uex,vex)=exch_UV_N_llc(u,v,N);
elseif u.grTopo=="cs";
  (uex,vex)=exch_UV_N_cs(u,v,N);
elseif u.grTopo=="ll";
  (uex,vex)=exch_UV_N_ll(u,v,N);
else;
  error("unknown grTopo case");
end;

return uex,vex

end

## Grid-specific implementations

function exch_T_N_ll(fld,N)

s=size(fld.f[1]); s=[i for i in s]; s[1:2]=s[1:2].+2*N;
f1=NaN*zeros(Float32,s[1],s[2]);
f1=zeros(Float32,s[1],s[2]);
FLD=gcmfaces(1,"ll",[f1]);

n3=max(size(fld.f[1],3),1); n4=max(size(fld.f[1],4),1);

#initial rotation for "LATLON" faces:
f1=copy(fld.f[1]);

#nan0=NaN*ones(N,size(fld.f[1],2),n3,n4);
nan0=NaN*ones(size(FLD.f[1],1),N,n3,n4);

#face 1:
F1=cat(f1[end-N+1:end,:,:,:],f1,f1[1:N,:,:,:];dims=1);
F1=cat(nan0,F1,nan0;dims=2);
#F1=copy(nan0);

F1=dropdims(F1,dims=(3,4));

#store:
FLD.f[1]=F1;

FLD

end

#

function exch_T_N_llc(fld,N)

s=size(fld.f[1]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f1=NaN*zeros(s[1],s[2]);
s=size(fld.f[2]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f2=NaN*zeros(s[1],s[2]);
s=size(fld.f[3]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f3=NaN*zeros(s[1],s[2]);
s=size(fld.f[4]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f4=NaN*zeros(s[1],s[2]);
s=size(fld.f[5]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f5=NaN*zeros(s[1],s[2]);
FLD=gcmfaces(5,"llc",[f1,f2,f3,f4,f5]);

n3=max(size(fld.f[1],3),1); n4=max(size(fld.f[1],4),1);

#initial rotation for "LATLON" faces:
f1=copy(fld.f[1]);
f2=copy(fld.f[2]);
f4=rotr90(fld.f[4]);
f5=rotr90(fld.f[5]);

nan1=NaN*ones(size(FLD.f[1],1),N,n3,n4);
nan2=NaN*ones(N,N,n3,n4);
#face 1:
F1=cat(f5[end-N+1:end,:,:,:],f1,f2[1:N,:,:,:];dims=1);
f3=rotl90(copy(fld.f[3]));
F1=cat(nan1,F1,cat(nan2,f3[:,1:N,:,:],nan2;dims=1);dims=2);
#face 2:
F2=cat(f1[end-N+1:end,:,:,:],f2,f4[1:N,:,:,:];dims=1);
f3=copy(fld.f[3]);
F2=cat(nan1,F2,cat(nan2,f3[:,1:N,:,:],nan2;dims=1);dims=2);
#face 4:
F4=cat(f2[end-N+1:end,:,:,:],f4,f5[1:N,:,:,:];dims=1);
f3=rotr90(copy(fld.f[3]));
F4=cat(nan1,F4,cat(nan2,f3[:,1:N,:,:],nan2;dims=1);dims=2);
#face 5:
F5=cat(f4[end-N+1:end,:,:,:],f5,f1[1:N,:,:,:];dims=1);
f3=rot180(copy(fld.f[3]));
F5=cat(nan1,F5,cat(nan2,f3[:,1:N,:,:],nan2;dims=1);dims=2);
#face 3:
f3=copy(fld.f[3]); F3=FLD.f[3];
F3[1+N:end-N,1+N:end-N,:,:]=f3;
F3=rotl90(F3); F3[1+N:end-N,1:N,:,:]=f1[:,end-N+1:end,:,:]; F3=rotr90(F3);
F3[1+N:end-N,1:N,:,:]=f2[:,end-N+1:end,:,:];
F3=rotr90(F3); F3[1+N:end-N,1:N,:,:]=f4[:,end-N+1:end,:,:]; F3=rotl90(F3);
F3=rot180(F3); F3[1+N:end-N,1:N,:,:]=f5[:,end-N+1:end,:,:]; F3=rot180(F3);

F1=dropdims(F1,dims=(3,4));
F2=dropdims(F2,dims=(3,4));
#F3=dropdims(F3,dims=(3,4));
F4=dropdims(F4,dims=(3,4));
F5=dropdims(F5,dims=(3,4));

#final rotation for "LATLON" faces:
F4=rotl90(F4);
F5=rotl90(F5);

#store:
FLD.f[1]=F1;
FLD.f[2]=F2;
FLD.f[3]=F3;
FLD.f[4]=F4;
FLD.f[5]=F5;

FLD

end

#

function exch_T_N_cs(fld,N)

s=size(fld.f[1]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f1=NaN*zeros(s[1],s[2]);
s=size(fld.f[2]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f2=NaN*zeros(s[1],s[2]);
s=size(fld.f[3]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f3=NaN*zeros(s[1],s[2]);
s=size(fld.f[4]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f4=NaN*zeros(s[1],s[2]);
s=size(fld.f[5]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f5=NaN*zeros(s[1],s[2]);
s=size(fld.f[6]); s=[i for i in s]; s[1:2]=s[1:2].+2*N; f6=NaN*zeros(s[1],s[2]);
FLD=gcmfaces(6,"cs",[f1,f2,f3,f4,f5,f6]);

n3=max(size(fld.f[1],3),1); n4=max(size(fld.f[1],4),1);

#initial rotation for "LATLON" faces:
f1=copy(fld.f[1]);
f2=copy(fld.f[2]);
f4=rotr90(fld.f[4]);
f5=rotr90(fld.f[5]);

nan1=NaN*ones(size(FLD.f[1],1),N,n3,n4);
nan2=NaN*ones(N,N,n3,n4);
#face 1:
F1=cat(f5[end-N+1:end,:,:,:],f1,f2[1:N,:,:,:];dims=1);
f3=rotl90(copy(fld.f[3]));
f6=copy(fld.f[6]);
F1=cat(cat(nan2,f6[:,end-N+1:end,:,:],nan2;dims=1),F1,cat(nan2,f3[:,1:N,:,:],nan2;dims=1);dims=2);
#face 2:
F2=cat(f1[end-N+1:end,:,:,:],f2,f4[1:N,:,:,:];dims=1);
f3=copy(fld.f[3]);
f6=rotl90(copy(fld.f[6]));
F2=cat(cat(nan2,f6[:,end-N+1:end,:,:],nan2;dims=1),F2,cat(nan2,f3[:,1:N,:,:],nan2;dims=1);dims=2);
#face 4:
F4=cat(f2[end-N+1:end,:,:,:],f4,f5[1:N,:,:,:];dims=1);
f3=rotr90(copy(fld.f[3]));
f6=rot180(copy(fld.f[6]));
F4=cat(cat(nan2,f6[:,end-N+1:end,:,:],nan2;dims=1),F4,cat(nan2,f3[:,1:N,:,:],nan2;dims=1);dims=2);
#face 5:
F5=cat(f4[end-N+1:end,:,:,:],f5,f1[1:N,:,:,:];dims=1);
f3=rot180(copy(fld.f[3]));
f6=rotr90(copy(fld.f[6]));
F5=cat(cat(nan2,f6[:,end-N+1:end,:,:],nan2;dims=1),F5,cat(nan2,f3[:,1:N,:,:],nan2;dims=1);dims=2);
#face 3:
f3=copy(fld.f[3]); F3=FLD.f[3];
F3[1+N:end-N,1+N:end-N,:,:]=f3;
F3=rotl90(F3); F3[1+N:end-N,1:N,:,:]=f1[:,end-N+1:end,:,:]; F3=rotr90(F3);
F3[1+N:end-N,1:N,:,:]=f2[:,end-N+1:end,:,:];
F3=rotr90(F3); F3[1+N:end-N,1:N,:,:]=f4[:,end-N+1:end,:,:]; F3=rotl90(F3);
F3=rot180(F3); F3[1+N:end-N,1:N,:,:]=f5[:,end-N+1:end,:,:]; F3=rot180(F3);
#face 6:
f6=copy(fld.f[6]); F6=FLD.f[6];
F6[1+N:end-N,1+N:end-N,:,:]=f6;
F6[1+N:end-N,end-N+1:end,:,:]=f1[:,1:N,:,:];
F6=rotl90(F6); F6[1+N:end-N,end-N+1:end,:,:]=f2[:,1:N,:,:]; F6=rotr90(F6);
F6=rot180(F6); F6[1+N:end-N,end-N+1:end,:,:]=f4[:,1:N,:,:]; F6=rot180(F6);
F6=rotr90(F6); F6[1+N:end-N,end-N+1:end,:,:]=f5[:,1:N,:,:]; F6=rotl90(F6);

F1=dropdims(F1,dims=(3,4));
F2=dropdims(F2,dims=(3,4));
#F3=dropdims(F3,dims=(3,4));
F4=dropdims(F4,dims=(3,4));
F5=dropdims(F5,dims=(3,4));
#F6=dropdims(F6,dims=(3,4));

#final rotation for "SIDE" faces:
F4=rotl90(F4);
F5=rotl90(F5);

#store:
FLD.f[1]=F1;
FLD.f[2]=F2;
FLD.f[3]=F3;
FLD.f[4]=F4;
FLD.f[5]=F5;
FLD.f[6]=F6;

FLD

end

## function that extends U,V in one direction each

function exch_UV_ll(fldU,fldV);

FLDUtmp=exchange(fldU);
FLDVtmp=exchange(fldV);

FLDU=exchange(fldU);
FLDV=exchange(fldV);

FLDU.f[1]=FLDUtmp.f[1][2:end,2:end-1];
FLDV.f[1]=FLDVtmp.f[1][2:end-1,2:end];

return FLDU,FLDV

end


function exch_UV_llc(fldU,fldV);

FLDUtmp=exchange(fldU);
FLDVtmp=exchange(fldV);

FLDU=exchange(fldU);
FLDV=exchange(fldV);

FLDU.f[1]=FLDUtmp.f[1][2:end,2:end-1];
FLDV.f[1]=FLDVtmp.f[1][2:end-1,2:end];
FLDV.f[1][:,end]=FLDUtmp.f[1][2:end-1,end];

FLDU.f[2]=FLDUtmp.f[2][2:end,2:end-1];
FLDU.f[2][end,:]=FLDVtmp.f[2][end,2:end-1];
FLDV.f[2]=FLDVtmp.f[2][2:end-1,2:end];

FLDU.f[3]=FLDUtmp.f[3][2:end,2:end-1];
FLDV.f[3]=FLDVtmp.f[3][2:end-1,2:end];
FLDV.f[3][:,end]=FLDUtmp.f[3][2:end-1,end];

FLDU.f[4]=FLDUtmp.f[4][2:end,2:end-1];
FLDV.f[4]=FLDVtmp.f[4][2:end-1,2:end];

FLDU.f[5]=FLDUtmp.f[5][2:end,2:end-1];
FLDV.f[5]=FLDVtmp.f[5][2:end-1,2:end];
FLDV.f[5][:,end]=FLDUtmp.f[5][2:end-1,end];

return FLDU,FLDV

end

#

function exch_UV_cs(fldU,fldV);

#!! this differs from regular exhange that adds
#points on all sides -- this just adds one point
#at end of one direction (U or V)
#
#should utlimately be avoided, in favor of
#exch_UV_N_cs, to keep things simple throughout?

FLDUtmp=exchange(fldU);
FLDVtmp=exchange(fldV);

FLDU=exchange(fldU);
FLDV=exchange(fldV);

FLDU.f[1]=FLDUtmp.f[1][2:end,2:end-1];
FLDV.f[1]=FLDVtmp.f[1][2:end-1,2:end];
FLDV.f[1][:,end]=FLDUtmp.f[1][2:end-1,end];

FLDU.f[2]=FLDUtmp.f[2][2:end,2:end-1];
FLDU.f[2][end,:]=FLDVtmp.f[2][end,2:end-1];
FLDV.f[2]=FLDVtmp.f[2][2:end-1,2:end];

FLDU.f[3]=FLDUtmp.f[3][2:end,2:end-1];
FLDV.f[3]=FLDVtmp.f[3][2:end-1,2:end];
FLDV.f[3][:,end]=FLDUtmp.f[3][2:end-1,end];

FLDU.f[4]=FLDUtmp.f[4][2:end,2:end-1];
FLDU.f[4][end,:]=FLDVtmp.f[4][end,2:end-1];
FLDV.f[4]=FLDVtmp.f[4][2:end-1,2:end];
#?? u

FLDU.f[5]=FLDUtmp.f[5][2:end,2:end-1];
FLDV.f[5]=FLDVtmp.f[5][2:end-1,2:end];
FLDV.f[5][:,end]=FLDUtmp.f[5][2:end-1,end];

FLDU.f[6]=FLDUtmp.f[6][2:end,2:end-1];
FLDU.f[6][end,:]=FLDVtmp.f[6][end,2:end-1];
FLDV.f[6]=FLDVtmp.f[6][2:end-1,2:end];
#??u

return FLDU,FLDV

end

## function that extends U, V by N points on all sides

function exch_UV_N_ll(fldU,fldV,N);

FLDU=exchange(fldU,N);
FLDV=exchange(fldV,N);

return FLDU,FLDV

end

#

function exch_UV_N_llc(fldU,fldV,N);

#fill interior of extended arrays
FLDUtmp=exch_T_N(fldU,N); FLDVtmp=exch_T_N(fldV,N);
for iFace=1:fldU.nFaces;
    FLDUtmp.f[iFace][1:N,:].=NaN;
    FLDUtmp.f[iFace][end-N+1:end,:].=NaN;
    FLDUtmp.f[iFace][:,1:N].=NaN;
    FLDUtmp.f[iFace][:,end-N+1:end].=NaN;
    FLDVtmp.f[iFace][1:N,:].=NaN;
    FLDVtmp.f[iFace][end-N+1:end,:].=NaN;
    FLDVtmp.f[iFace][:,1:N].=NaN;
    FLDVtmp.f[iFace][:,end-N+1:end].=NaN;
end;

FLDU=FLDUtmp;
FLDV=FLDVtmp;

#now add the one extra point we need for U and V
(FLDUtmp,FLDVtmp)=exch_UV_llc(fldU,fldV);
#then add the remaining rows and columns...

#U face 1
tmp1=permutedims(FLDVtmp.f[5][:,end-N:end-1,:],[2 1 3]);
FLDU.f[1][1:N,N+1:end-N]=reverse(tmp1,dims=2);
tmp1=FLDUtmp.f[2][1:N,:];
FLDU.f[1][end-N+1:end,N+1:end-N]=tmp1;
#
tmp1=permutedims(-FLDVtmp.f[3][1:N,:,:],[2 1 3]);
FLDU.f[1][N+1:end-N+1,end-N+1:end]=reverse(tmp1,dims=1);

#V face 1
tmp1=permutedims(-FLDUtmp.f[5][:,end-N+1:end,:],[2 1 3]);
FLDV.f[1][1:N,N+1:end-N+1]=reverse(tmp1,dims=2);
tmp1=FLDVtmp.f[2][1:N,:];
FLDV.f[1][end-N+1:end,N+1:end-N+1]=tmp1;
#
tmp1=permutedims(FLDUtmp.f[3][1:N,:,:],[2 1 3]);
FLDV.f[1][N+1:end-N,end-N+1:end]=reverse(tmp1,dims=1);

#U face 2
tmp1=FLDUtmp.f[1][end-N:end-1,:];
FLDU.f[2][1:N,N+1:end-N]=tmp1;
tmp1=permutedims(FLDVtmp.f[4][:,1:N,:],[2 1 3]);
FLDU.f[2][end-N+1:end,N+1:end-N]=reverse(tmp1,dims=2);
#
tmp1=FLDUtmp.f[3][:,1:N];
FLDU.f[2][N+1:end-N+1,end-N+1:end]=tmp1;

#V face 2
tmp1=FLDVtmp.f[1][end-N+1:end,:];
FLDV.f[2][1:N,N+1:end-N+1]=tmp1;
tmp1=permutedims(-FLDUtmp.f[4][:,1:N,:],[2 1 3]);
FLDV.f[2][end-N+1:end,N+1:end-N+1]=reverse(tmp1,dims=2);
#
tmp1=FLDVtmp.f[3][:,1:N];
FLDV.f[2][N+1:end-N,end-N+1:end]=tmp1;

#U face 3
tmp1=permutedims(FLDVtmp.f[1][:,end-N:end-1,:],[2 1 3]);
FLDU.f[3][1:N,N+1:end-N]=reverse(tmp1,dims=2);
tmp1=FLDUtmp.f[4][1:N,:];
FLDU.f[3][end-N+1:end,N+1:end-N]=tmp1;
#
tmp1=FLDUtmp.f[2][:,end-N+1:end];
FLDU.f[3][N+1:end-N+1,1:N]=tmp1;
tmp1=permutedims(-FLDVtmp.f[5][1:N,:,:],[2 1 3]);
FLDU.f[3][N+1:end-N+1,end-N+1:end]=reverse(tmp1,dims=1);

#V face 3
tmp1=permutedims(-FLDUtmp.f[1][:,end-N+1:end,:],[2 1 3]);
FLDV.f[3][1:N,N+1:end-N+1]=reverse(tmp1,dims=2);
tmp1=FLDVtmp.f[4][1:N,:];
FLDV.f[3][end-N+1:end,N+1:end-N+1]=tmp1;
#
tmp1=FLDVtmp.f[2][:,end-N:end-1];
FLDV.f[3][N+1:end-N,1:N]=tmp1;
tmp1=permutedims(FLDUtmp.f[5][1:N,:,:],[2 1 3]);
FLDV.f[3][N+1:end-N,end-N+1:end]=reverse(tmp1,dims=1);

#U face 4
tmp1=FLDUtmp.f[3][end-N:end-1,:];
FLDU.f[4][1:N,N+1:end-N]=tmp1;
#
tmp1=permutedims(-FLDVtmp.f[2][end-N+1:end,:,:],[2 1 3]);
FLDU.f[4][N+1:end-N+1,1:N]=reverse(tmp1,dims=1);
tmp1=FLDUtmp.f[5][:,1:N];
FLDU.f[4][N+1:end-N+1,end-N+1:end]=tmp1;

#V face 4
tmp1=FLDVtmp.f[3][end-N+1:end,:];
FLDV.f[4][1:N,N+1:end-N+1]=tmp1;
#
tmp1=permutedims(FLDUtmp.f[2][end-N:end-1,:,:],[2 1 3]);
FLDV.f[4][N+1:end-N,1:N]=reverse(tmp1,dims=1);
tmp1=FLDVtmp.f[5][:,1:N];
FLDV.f[4][N+1:end-N,end-N+1:end]=tmp1;

#U face 5
tmp1=permutedims(FLDVtmp.f[3][:,end-N:end-1,:],[2 1 3]);
FLDU.f[5][1:N,N+1:end-N]=reverse(tmp1,dims=2);
#
tmp1=FLDUtmp.f[4][:,end-N+1:end];
FLDU.f[5][N+1:end-N+1,1:N]=tmp1;
tmp1=permutedims(-FLDVtmp.f[1][1:N,:,:],[2 1 3]);
FLDU.f[5][N+1:end-N+1,end-N+1:end]=reverse(tmp1,dims=1);

#V face 5
tmp1=permutedims(-FLDUtmp.f[3][:,end-N+1:end,:],[2 1 3]);
FLDV.f[5][1:N,N+1:end-N+1]=reverse(tmp1,dims=2);
#
tmp1=FLDVtmp.f[4][:,end-N:end-1];
FLDV.f[5][N+1:end-N,1:N]=tmp1;
tmp1=permutedims(FLDUtmp.f[1][1:N,:,:],[2 1 3]);
FLDV.f[5][N+1:end-N,end-N+1:end]=reverse(tmp1,dims=1);

return FLDU,FLDV

end

#

function exch_UV_N_cs(fldU,fldV,N);

#fill interior of extended arrays
FLDUtmp=exchange(fldU,N);
FLDVtmp=exchange(fldV,N);
for iFace=1:fldU.nFaces;
    FLDUtmp.f[iFace][1:N,:,:].=NaN;
    FLDUtmp.f[iFace][end-N+1:end,:,:].=NaN;
    FLDUtmp.f[iFace][:,1:N,:].=NaN;
    FLDUtmp.f[iFace][:,end-N+1:end,:].=NaN;
    FLDVtmp.f[iFace][1:N,:,:].=NaN;
    FLDVtmp.f[iFace][end-N+1:end,:,:].=NaN;
    FLDVtmp.f[iFace][:,1:N,:].=NaN;
    FLDVtmp.f[iFace][:,end-N+1:end,:].=NaN;
end;

FLDU=FLDUtmp;
FLDV=FLDVtmp;

#now add the one extra point we need for U and V
(FLDUtmp,FLDVtmp)=exch_UV_cs(fldU,fldV);

#then add the remaining rows and columns
#... got to this point on 2018/08/14

#U face 1
    tmp1=permutedims(FLDVtmp.f[5][:,end-N:end-1,:],[2 1 3]);
    FLDU.f[1][1:N,N+1:end-N,:]=reverse(tmp1,dims=2);
    tmp1=FLDUtmp.f[2][1:N,:,:];
    FLDU.f[1][end-N+1:end,N+1:end-N,:]=tmp1;
#
    tmp1=FLDUtmp.f[6][:,end-N+1:end,:];
    FLDU.f[1][N+1:end-N+1,1:N,:]=tmp1;
    tmp1=permutedims(-FLDVtmp.f[3][1:N,:,:],[2 1 3]);
    FLDU.f[1][N+1:end-N+1,end-N+1:end,:]=reverse(tmp1,dims=1);

#V face 1
    tmp1=permutedims(-FLDUtmp.f[5][:,end-N+1:end,:],[2 1 3]);
    FLDV.f[1][1:N,N+1:end-N+1,:]=reverse(tmp1,dims=2);
    tmp1=FLDVtmp.f[2][1:N,:,:];
    FLDV.f[1][end-N+1:end,N+1:end-N+1,:]=tmp1;
#
    tmp1=FLDVtmp.f[6][:,end-N:end-1,:];
    FLDV.f[1][N+1:end-N,1:N,:]=tmp1;
    tmp1=permutedims(FLDUtmp.f[3][1:N,:,:],[2 1 3]);
    FLDV.f[1][N+1:end-N,end-N+1:end,:]=reverse(tmp1,dims=1);

#U face 2
    tmp1=FLDUtmp.f[1][end-N:end-1,:,:];
    FLDU.f[2][1:N,N+1:end-N,:]=tmp1;
    tmp1=permutedims(FLDVtmp.f[4][:,1:N,:],[2 1 3]);
    FLDU.f[2][end-N+1:end,N+1:end-N,:]=reverse(tmp1,dims=2);
#
    tmp1=permutedims(-FLDVtmp.f[6][end-N+1:end,:,:],[2 1 3]);
    FLDU.f[2][N+1:end-N+1,1:N,:]=reverse(tmp1,dims=1);
    tmp1=FLDUtmp.f[3][:,1:N,:];
    FLDU.f[2][N+1:end-N+1,end-N+1:end,:]=tmp1;

#V face 2
    tmp1=FLDVtmp.f[1][end-N+1:end,:,:];
    FLDV.f[2][1:N,N+1:end-N+1,:]=tmp1;
    tmp1=permutedims(-FLDUtmp.f[4][:,1:N,:],[2 1 3]);
    FLDV.f[2][end-N+1:end,N+1:end-N+1,:]=reverse(tmp1,dims=2);
#
    tmp1=permutedims(FLDUtmp.f[6][end-N:end-1,:,:],[2 1 3]);
    FLDV.f[2][N+1:end-N,1:N,:]=reverse(tmp1,dims=1);
    tmp1=FLDVtmp.f[3][:,1:N,:];
    FLDV.f[2][N+1:end-N,end-N+1:end,:]=tmp1;

#U face 3
    tmp1=permutedims(FLDVtmp.f[1][:,end-N:end-1,:],[2 1 3]);
    FLDU.f[3][1:N,N+1:end-N,:]=reverse(tmp1,dims=2);
    tmp1=FLDUtmp.f[4][1:N,:,:];
    FLDU.f[3][end-N+1:end,N+1:end-N,:]=tmp1;
#
    tmp1=FLDUtmp.f[2][:,end-N+1:end,:];
    FLDU.f[3][N+1:end-N+1,1:N,:]=tmp1;
    tmp1=permutedims(-FLDVtmp.f[5][1:N,:,:],[2 1 3]);
    FLDU.f[3][N+1:end-N+1,end-N+1:end,:]=reverse(tmp1,dims=1);

#V face 3
    tmp1=permutedims(-FLDUtmp.f[1][:,end-N+1:end,:],[2 1 3]);
    FLDV.f[3][1:N,N+1:end-N+1,:]=reverse(tmp1,dims=2);
    tmp1=FLDVtmp.f[4][1:N,:,:];
    FLDV.f[3][end-N+1:end,N+1:end-N+1,:]=tmp1;
#
    tmp1=FLDVtmp.f[2][:,end-N:end-1,:];
    FLDV.f[3][N+1:end-N,1:N,:]=tmp1;
    tmp1=permutedims(FLDUtmp.f[5][1:N,:,:],[2 1 3]);
    FLDV.f[3][N+1:end-N,end-N+1:end,:]=reverse(tmp1,dims=1);

#U face 4
    tmp1=FLDUtmp.f[3][end-N:end-1,:,:];
    FLDU.f[4][1:N,N+1:end-N,:]=tmp1;
    tmp1=permutedims(FLDVtmp.f[6][:,1:N,:],[2 1 3]);
    FLDU.f[4][end-N+1:end,N+1:end-N,:]=reverse(tmp1,dims=2);
#
    tmp1=permutedims(-FLDVtmp.f[2][end-N+1:end,:,:],[2 1 3]);
    FLDU.f[4][N+1:end-N+1,1:N,:]=reverse(tmp1,dims=1);
    tmp1=FLDUtmp.f[5][:,1:N,:];
    FLDU.f[4][N+1:end-N+1,end-N+1:end,:]=tmp1;

#V face 4
    tmp1=FLDVtmp.f[3][end-N+1:end,:,:];
    FLDV.f[4][1:N,N+1:end-N+1,:]=tmp1;
    tmp1=permutedims(-FLDUtmp.f[6][:,1:N,:],[2 1 3]);
    FLDV.f[4][end-N+1:end,N+1:end-N+1,:]=reverse(tmp1,dims=2);
#
    tmp1=permutedims(FLDUtmp.f[2][end-N:end-1,:,:],[2 1 3]);
    FLDV.f[4][N+1:end-N,1:N,:]=reverse(tmp1,dims=1);
    tmp1=FLDVtmp.f[5][:,1:N,:];
    FLDV.f[4][N+1:end-N,end-N+1:end,:]=tmp1;

#U face 5
    tmp1=permutedims(FLDVtmp.f[3][:,end-N:end-1,:],[2 1 3]);
    FLDU.f[5][1:N,N+1:end-N,:]=reverse(tmp1,dims=2);
    tmp1=FLDUtmp.f[6][1:N,:,:];
    FLDU.f[5][end-N+1:end,N+1:end-N,:]=tmp1;
#
    tmp1=FLDUtmp.f[4][:,end-N+1:end,:];
    FLDU.f[5][N+1:end-N+1,1:N,:]=tmp1;
    tmp1=permutedims(-FLDVtmp.f[1][1:N,:,:],[2 1 3]);
    FLDU.f[5][N+1:end-N+1,end-N+1:end,:]=reverse(tmp1,dims=1);

#V face 5
    tmp1=permutedims(-FLDUtmp.f[3][:,end-N+1:end,:],[2 1 3]);
    FLDV.f[5][1:N,N+1:end-N+1,:]=reverse(tmp1,dims=2);
    tmp1=FLDVtmp.f[6][1:N,:,:];
    FLDV.f[5][end-N+1:end,N+1:end-N+1,:]=tmp1;
#
    tmp1=FLDVtmp.f[4][:,end-N:end-1,:];
    FLDV.f[5][N+1:end-N,1:N,:]=tmp1;
    tmp1=permutedims(FLDUtmp.f[1][1:N,:,:],[2 1 3]);
    FLDV.f[5][N+1:end-N,end-N+1:end,:]=reverse(tmp1,dims=1);

#U face 6
    tmp1=FLDUtmp.f[5][end-N:end-1,:,:];
    FLDU.f[6][1:N,N+1:end-N,:]=tmp1;
    tmp1=permutedims(FLDVtmp.f[2][:,1:N,:],[2 1 3]);
    FLDU.f[6][end-N+1:end,N+1:end-N,:]=reverse(tmp1,dims=2);
#
    tmp1=permutedims(-FLDVtmp.f[4][end-N+1:end,:,:],[2 1 3]);
    FLDU.f[6][N+1:end-N+1,1:N,:]=reverse(tmp1,dims=1);
    tmp1=FLDUtmp.f[1][:,1:N,:];
    FLDU.f[6][N+1:end-N+1,end-N+1:end,:]=tmp1;

#V face 6
    tmp1=FLDVtmp.f[5][end-N+1:end,:,:];
    FLDV.f[6][1:N,N+1:end-N+1,:]=tmp1;
    tmp1=permutedims(-FLDUtmp.f[2][:,1:N,:],[2 1 3]);
    FLDV.f[6][end-N+1:end,N+1:end-N+1,:]=reverse(tmp1,dims=2);
#
    tmp1=permutedims(FLDUtmp.f[4][end-N:end-1,:,:],[2 1 3]);
    FLDV.f[6][N+1:end-N,1:N,:]=reverse(tmp1,dims=1);
    tmp1=FLDVtmp.f[1][:,1:N,:];
    FLDV.f[6][N+1:end-N,end-N+1:end,:]=tmp1;

    return FLDU,FLDV

end
