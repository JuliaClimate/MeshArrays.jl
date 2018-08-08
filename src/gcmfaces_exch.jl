# exchange.jl
# exch_T_N.jl
# exch_T_N_llc.jl
# exch_T_N_cs.jl
# exch_T_N_ll.jl
#
#       First Draft Implementation
#
# gaelforget (https://github.com/gaelforget/gcmfaces_jl)
# Julia 0.6.2
# Created: 02.01.18
# Last Edit: 02.01.18

#correspondance to sym_g codes:
#5	rotl90
#6	rot180
#7	rotr90
#
#only works in 2D? do I need to translate sym_g?

function exchange(fld,N)
  FLD=exch_T_N(fld,N);
end

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

function exch_T_N_ll(fld,N)

s=size(fld.f[1]); s=[i for i in s]; s[1:2]=s[1:2]+2*N;
f1=NaN*zeros(Float32,s[1],s[2]);
f1=zeros(Float32,s[1],s[2]);
FLD=gcmfaces(1,"ll",[f1]);

n3=max(size(fld.f[1],3),1); n4=max(size(fld.f[1],4),1);

#initial rotation for "LATLON" faces:
f1=copy(fld.f[1]);

#nan0=NaN*ones(N,size(fld.f[1],2),n3,n4);
nan0=NaN*ones(size(FLD.f[1],1),N,n3,n4);

#face 1:
F1=cat(1,f1[end-N+1:end,:,:,:],f1,f1[1:N,:,:,:]);
F1=cat(2,nan0,F1,nan0);
#F1=copy(nan0);

F1=squeeze(F1,(3,4));

#store:
FLD.f[1]=F1;

FLD

end

function exch_T_N_llc(fld,N)

s=size(fld.f[1]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f1=NaN*zeros(s[1],s[2]);
s=size(fld.f[2]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f2=NaN*zeros(s[1],s[2]);
s=size(fld.f[3]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f3=NaN*zeros(s[1],s[2]);
s=size(fld.f[4]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f4=NaN*zeros(s[1],s[2]);
s=size(fld.f[5]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f5=NaN*zeros(s[1],s[2]);
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
F1=cat(1,f5[end-N+1:end,:,:,:],f1,f2[1:N,:,:,:]);
f3=rotl90(copy(fld.f[3]));
F1=cat(2,nan1,F1,cat(1,nan2,f3[:,1:N,:,:],nan2));
#face 2:
F2=cat(1,f1[end-N+1:end,:,:,:],f2,f4[1:N,:,:,:]);
f3=copy(fld.f[3]);
F2=cat(2,nan1,F2,cat(1,nan2,f3[:,1:N,:,:],nan2));
#face 4:
F4=cat(1,f2[end-N+1:end,:,:,:],f4,f5[1:N,:,:,:]);
f3=rotr90(copy(fld.f[3]));
F4=cat(2,nan1,F4,cat(1,nan2,f3[:,1:N,:,:],nan2));
#face 5:
F5=cat(1,f4[end-N+1:end,:,:,:],f5,f1[1:N,:,:,:]);
f3=rot180(copy(fld.f[3]));
F5=cat(2,nan1,F5,cat(1,nan2,f3[:,1:N,:,:],nan2));
#face 3:
f3=copy(fld.f[3]); F3=FLD.f[3];
F3[1+N:end-N,1+N:end-N,:,:]=f3;
F3=rotl90(F3); F3[1+N:end-N,1:N,:,:]=f1[:,end-N+1:end,:,:]; F3=rotr90(F3);
F3[1+N:end-N,1:N,:,:]=f2[:,end-N+1:end,:,:];
F3=rotr90(F3); F3[1+N:end-N,1:N,:,:]=f4[:,end-N+1:end,:,:]; F3=rotl90(F3);
F3=rot180(F3); F3[1+N:end-N,1:N,:,:]=f5[:,end-N+1:end,:,:]; F3=rot180(F3);

F1=squeeze(F1,(3,4));
F2=squeeze(F2,(3,4));
#F3=squeeze(F3,(3,4));
F4=squeeze(F4,(3,4));
F5=squeeze(F5,(3,4));

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

function exch_T_N_cs(fld,N)

s=size(fld.f[1]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f1=NaN*zeros(s[1],s[2]);
s=size(fld.f[2]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f2=NaN*zeros(s[1],s[2]);
s=size(fld.f[3]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f3=NaN*zeros(s[1],s[2]);
s=size(fld.f[4]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f4=NaN*zeros(s[1],s[2]);
s=size(fld.f[5]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f5=NaN*zeros(s[1],s[2]);
s=size(fld.f[6]); s=[i for i in s]; s[1:2]=s[1:2]+2*N; f6=NaN*zeros(s[1],s[2]);
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
F1=cat(1,f5[end-N+1:end,:,:,:],f1,f2[1:N,:,:,:]);
f3=rotl90(copy(fld.f[3]));
f6=copy(fld.f[6]);
F1=cat(2,cat(1,nan2,f6[:,end-N+1:end,:,:],nan2),F1,cat(1,nan2,f3[:,1:N,:,:],nan2));
#face 2:
F2=cat(1,f1[end-N+1:end,:,:,:],f2,f4[1:N,:,:,:]);
f3=copy(fld.f[3]);
f6=rotl90(copy(fld.f[6]));
F2=cat(2,cat(1,nan2,f6[:,end-N+1:end,:,:],nan2),F2,cat(1,nan2,f3[:,1:N,:,:],nan2));
#face 4:
F4=cat(1,f2[end-N+1:end,:,:,:],f4,f5[1:N,:,:,:]);
f3=rotr90(copy(fld.f[3]));
f6=rot180(copy(fld.f[6]));
F4=cat(2,cat(1,nan2,f6[:,end-N+1:end,:,:],nan2),F4,cat(1,nan2,f3[:,1:N,:,:],nan2));
#face 5:
F5=cat(1,f4[end-N+1:end,:,:,:],f5,f1[1:N,:,:,:]);
f3=rot180(copy(fld.f[3]));
f6=rotr90(copy(fld.f[6]));
F5=cat(2,cat(1,nan2,f6[:,end-N+1:end,:,:],nan2),F5,cat(1,nan2,f3[:,1:N,:,:],nan2));
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

F1=squeeze(F1,(3,4));
F2=squeeze(F2,(3,4));
#F3=squeeze(F3,(3,4));
F4=squeeze(F4,(3,4));
F5=squeeze(F5,(3,4));
#F6=squeeze(F6,(3,4));

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
