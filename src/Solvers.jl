# # Poisson And Laplace Equation Solvers On The Sphere
#
# Compute scalar potential, and divergent component, from vertically integrated transport.

using SparseArrays, Statistics

"""
    MaskWetPoints(TrspCon)

Mask land points with NaN.
```
(TrspCon, mskWet, mskDry)=MaskWetPoints(TrspCon)
```
"""
function MaskWetPoints(TrspCon)
    mskWet=1.0 .+ 0.0 * TrspCon
    mskDry=1.0 * isnan.(mskWet)
    mskDry=mask(mskDry,NaN,0.0)
    #
    tmp1=fill(1.0,mskWet); tmp2=exchange(tmp1);
    for I=1:size(tmp1.f,1)
        tmp3=mskWet[I]; tmp4=tmp2[I];
        tmp4=tmp4[2:end-1,1:end-2]+tmp4[2:end-1,3:end]+tmp4[1:end-2,2:end-1]+tmp4[3:end,2:end-1];
        !isempty(findall(isnan.(tmp4) .& (!isnan).(tmp3))) ? println("warning: modified mask") : nothing
        tmp3[findall(isnan.(tmp4))] .= NaN
        mskWet[I]=tmp3
    end
    #
    TrspCon=mask(TrspCon,0.0)*mskWet;
    return TrspCon, mskWet, mskDry
end

"""
    MapWetPoints(mskWet)

Mapping from global array to global ocean vector.
```
(Kvec,Lvec,Kmap,Lmap)=MapWetPoints(mskWet)
```
"""
function MapWetPoints(mskWet)
    tmp1=write(mskWet)[:]
    kk=findall((!isnan).(tmp1))
    nn=length(kk); s0=size(tmp1); s1=mskWet.grid.ioSize;
    Kvec=fill(0.0,s0...); Kvec[kk]=kk; Kmap=read(reshape(Kvec,s1...),mskWet) #global array indices
    Lvec=fill(0.0,s0...); Lvec[kk]=1:nn; Lmap=read(reshape(Lvec,s1...),mskWet) #global vector indices
    return Kvec,Lvec,Kmap,Lmap
end

"""
    SeedWetPoints(tmp::MeshArray,Kmap::MeshArray,Lmap::MeshArray,I...)

Seed a subset of grid points.
```
(FLDones,FLDkkFROM)=SeedWetPoints(tmp::MeshArray,Kmap::MeshArray,Lmap::MeshArray,I...)
```
"""
function SeedWetPoints(tmp::MeshArray,Kmap::MeshArray,Lmap::MeshArray,I...)
    aa=I[1]
    ii=I[2]
    jj=I[3]

    #1) seed points (FLDones) and neighborhood of influence (FLDkkFROM)
    FLDones=fill(0.0,tmp)
    FLDones[aa][ii:3:end,jj:3:end].=1.0
    FLDones[aa][findall(Kmap[aa].==0.0)].=0.0

    FLDkkFROMtmp=fill(0.0,tmp)
    FLDkkFROMtmp[aa][ii:3:end,jj:3:end]=Kmap[aa][ii:3:end,jj:3:end]
    FLDkkFROMtmp[aa][findall(isnan.(tmp[aa]))].=0.0

    FLDkkFROM=exchange(FLDkkFROMtmp)
    FLDkkFROM=mask(FLDkkFROM,0.0)

    for bb in 1:tmp.grid.nFaces
        tmp1=FLDkkFROM[bb]
        tmp2=zeros(size(tmp1) .- 2)
        for ii2 in 1:3
            for jj2 in 1:3
                tmp2=tmp2+tmp1[ii2:end-3+ii2,jj2:end-3+jj2]
            end
        end
        FLDkkFROM[bb]=tmp2
    end

    return FLDones,FLDkkFROM
end

"""
    MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)

Assemble sparse matrix using mskWet, Kvec, Lvec directly and Kmap, Lmap via SeedWetPoints
```
A=MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)
```
"""
function MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)
    I=Array{Int}(undef,0)
    J=Array{Int}(undef,0)
    V=Array{Float64}(undef,0)

    for aa=1:TrspCon.grid.nFaces;
        for ii=1:3; for jj=1:3;
            #1) compute effect of each point on neighboring target point:
            (FLDones,FLDkkFROM)=SeedWetPoints(TrspCon,Kmap,Lmap,aa,ii,jj);
            (tmpU,tmpV)=gradient(FLDones,Dict(),false)
            dFLDdt=convergence(tmpU,tmpV);
            #2) mask `dFLDdt` since we use a **Neumann** boundary condition.
            #Extrapolation uses a **Dirichlet** boundary condition, so mskFreeze should not be applied then.
            isa(FLDkkFROM,MeshArray) ? FLDkkFROM=write(FLDkkFROM)[:] : nothing;
            #3.1) For wet points -- add contributions in matrix:
            dFLDdtWet=write(dFLDdt.*mskWet)[:];
            tmp1=findall( (dFLDdtWet .!= 0.0) .* (!isnan).(dFLDdtWet));
            tmpV=dFLDdtWet[tmp1]; tmpJ=FLDkkFROM[tmp1]; tmpI=Kvec[tmp1];
            I=[I;Lvec[Int.(tmpI)]]; J=[J;Lvec[Int.(tmpJ)]]; V=[V;tmpV];
            size(tmpV)
            #3.2) For dry points -- This part reflects the `Neumann` boundary condition:
            dFLDdtDry=write(dFLDdt.*mskDry)[:];
            tmp1=findall( (dFLDdtDry .!= 0.0) .* (!isnan).(dFLDdtDry));
            tmpV=dFLDdtDry[tmp1]; tmpIJ=FLDkkFROM[tmp1];
            I=[I;Lvec[Int.(tmpIJ)]]; J=[J;Lvec[Int.(tmpIJ)]]; V=[V;tmpV];
        end; end;
    end;

    nn=sum((!isnan).(mskWet))
    A=sparse(I,J,V,nn,nn)
end

"""
    ScalarPotential(TrspCon)

Scalar potential inversion.
```
TrspPot=ScalarPotential(TrspCon)
```
"""
function ScalarPotential(TrspCon)
    (TrspCon, mskWet, mskDry)=MaskWetPoints(TrspCon)
    (Kvec,Lvec,Kmap,Lmap)=MapWetPoints(mskWet)

    A=MatrixForPoisson(TrspCon,mskWet,mskDry,Kvec,Lvec,Kmap,Lmap)
    yy=write(TrspCon)[:]
    yy=yy[findall(Kvec .!= 0)]

    xx=A\yy
    xx=xx .- median(xx)

    TrspPot=0.0 * write(TrspCon)[:]
    TrspPot[findall(Kvec .!= 0)]=xx
    TrspPot=read(TrspPot,TrspCon)

    return TrspPot
end

"""
    VectorPotential(TrspX,TrspY,GridVariables,method::Int=1)

Vector potential inversion.
```
TrspPot=ScalarPotential(TrspCon)
```
"""
function VectorPotential(TrspX::MeshArray,TrspY::MeshArray,GridVariables::Dict,method::Int=1)

    # 1)  streamfunction face by face:

    (fldU,fldV)=exch_UV(TrspX,TrspY);
    fldU=mask(fldU,0.0); fldV=mask(fldV,0.0);

    psi=similar(fldV)
    for I in eachindex(fldV)
        tmp2=cumsum(fldV[I],dims=1)
        psi[I]=[zeros(1,size(tmp2,2));tmp2]
    end

    # 2a)  reset one land value per face to zero (method 1)

    if method==1
        mskW=mask(1.0 .+ 0.0 * mask(view(GridVariables["hFacW"],:,1),NaN,0.0),0.0)
        mskS=mask(1.0 .+ 0.0 * mask(view(GridVariables["hFacS"],:,1),NaN,0.0),0.0)
        (mskW,mskS)=exch_UV(mskW,mskS); mskW=abs.(mskW); mskS=abs.(mskS)

        for iF=1:TrspX.grid.nFaces
            tmp2=mskS[iF]
            tmp3a=[ones(1,size(tmp2,2));tmp2]
            tmp3b=[tmp2;ones(1,size(tmp2,2))]
            tmpB=tmp3a.*tmp3b
            #
            tmpA=psi[iF];
            I=findall(tmpB.==0)
            ii=I[1][1]; jj=I[1][2]
            tmpA[:,jj]=tmpA[:,jj] .- tmpB[ii,jj]
            #
            for kk=1:jj-1
                tmpE=tmpA[:,kk+1]-tmpA[:,kk]
                tmpE+=fldU[iF][:,kk]
                tmpE=median(tmpE[findall((!isnan).(tmpE))])
                tmpA[:,kk]=tmpA[:,kk] .- tmpE
            end
            #
            for kk=jj+1:size(tmpA,2)
                tmpE=tmpA[:,kk]-tmpA[:,kk-1]
                tmpE+=fldU[iF][:,kk-1]
                tmpE=median(tmpE[findall((!isnan).(tmpE))])
                tmpA[:,kk]=tmpA[:,kk] .- tmpE
            end;
            #
            psi[iF]=tmpA
        end
    end

    # 2b)  subtract divergent flow line by line (method 2)

    if method==2
        for I in eachindex(fldV)
            tmp2=diff(psi[I],dims=2)+fldU[I]
            tmp3=cumsum(mean(tmp2,dims=1),dims=2)
            psi[I]=psi[I]-ones(size(psi[I],1),1)*[0.0 tmp3]
        end
    end

    # 3) match edge values:

    if fldU.grid.nFaces>1
        TMP1=similar(psi)
        for I in eachindex(TMP1); TMP1[I] = fill(I,size(psi[I])); end
        TMP2=exchange(TMP1) #this is a trick

    for I in 1:TrspX.grid.nFaces-1
        tmp2=exchange(psi) #this is a trick
        tmp3=tmp2[I+1]; tmp3[3:end-2,3:end-2].=NaN #mask out interior points
        TMP3=TMP2[I+1]; tmp3[findall(TMP3.>I+1)].=NaN #mask out edges points coming from unadjusted faces
        tmp3[findall(TMP3.==0)].=NaN
        tmp3[:,1]=tmp3[:,1]-tmp3[:,2]; tmp3[:,end]=tmp3[:,end]-tmp3[:,end-1] #compare edge points
        tmp3[1,:]=tmp3[1,:]-tmp3[2,:]; tmp3[end,:]=tmp3[end,:]-tmp3[end-1,:] #compare edge points
        tmp3[2:end-1,2:end-1] .= NaN #mask out remaining interior points
        psi[I+1]=psi[I+1] .+ median(tmp3[findall((!isnan).(tmp3))]) #adjust the face data
    end
    end

    # 4) put streamfunction at cell center

    for I in 1:TrspX.grid.nFaces
        tmp3=psi[I]
        tmp3=(tmp3[:,1:end-1]+tmp3[:,2:end])/2
        tmp3=(tmp3[1:end-1,:]+tmp3[2:end,:])/2
        psi[I]=tmp3
    end

    # 5) reset one land point to zero

    #all values:
    tmp1=write(psi)

    #all land points:
    tmp2=write(view(GridVariables["hFacC"],:,1))
    tmp2=findall( (tmp2 .== 0.0) .& (!isnan).(tmp1))

    #closest to Boston:
    tmp_lon=write(GridVariables["XC"])[tmp2]
    tmp_lat=write(GridVariables["YC"])[tmp2]
    tmp_dis=(tmp_lat .- 42.3601).^2 + (tmp_lon .- 71.0589).^2
    tmp2=tmp2[findall(tmp_dis .== minimum(tmp_dis))]

    #set that point to zero
    psi=psi .- tmp1[tmp2[1]];

    return psi
end
