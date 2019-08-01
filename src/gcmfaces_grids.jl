
## GCMGridSpec function with default gridName argument:

GCMGridSpec() = GCMGridSpec("LLC90")

## GCMGridSpec function with gridName argument:

"""
    GCMGridSpec(gridName)

Set global variables in the module scope for grDir, nFaces, grTopo, ioSize,
facesSize, ioPrec using hard-coded values for LLC90, CS32, LL360 (for now).
"""
function GCMGridSpec(gridName,gridParentDir="./")

global grDir, nFaces, grTopo, ioSize, facesSize, ioPrec;

if gridName=="LLC90";
    grDir=gridParentDir*"GRID_LLC90/";
    nFaces=5;
    grTopo="llc";
    ioSize=[90 1170];
    facesSize=[(90, 270), (90, 270), (90, 90), (270, 90), (270, 90)]
    ioPrec=Float64;
elseif gridName=="CS32";
    grDir=gridParentDir*"GRID_CS32/";
    nFaces=6;
    grTopo="cs";
    ioSize=[32 192];
    facesSize=[(32, 32), (32, 32), (32, 32), (32, 32), (32, 32), (32, 32)]
    ioPrec=Float32;
elseif gridName=="LL360";
    grDir=gridParentDir*"GRID_LL360/";
    nFaces=1;
    grTopo="ll";
    ioSize=[360 160];
    facesSize=[(360, 160)]
    ioPrec=Float32;
else;
    error("unknown gridName case");
end;

mygrid=Dict("grDir" => grDir, "nFaces" => nFaces, "grTopo" => grTopo,
    "ioSize" => ioSize, "facesSize" => facesSize, "ioPrec" => ioPrec)
return mygrid

end

## GCMGridLoad function

"""
    GCMGridLoad(mygrid::Dict)

Loads grid variables from files located in grDir set by GCMGridSpec.

Grid variables are XC, XG, YC, YG, RAC, RAZ, DXC, DXG, DYC, DYG, hFacC,
hFacS, hFacW, Depth based on the MITgcm naming convention.
"""
function GCMGridLoad(mygrid::Dict=Dict())

    #maybe just return as a dictionnary?
    global XC, XG, YC, YG, RAC, RAZ, DXC, DXG, DYC, DYG, hFacC, hFacS, hFacW, Depth
    global AngleCS, AngleSN, RAW, RAS
    global DRC, DRF, RC, RF

    list0=("XC","XG","YC","YG","AngleCS","AngleSN","RAC","RAW","RAS","RAZ",
    "DXC","DXG","DYC","DYG","hFacC","hFacS","hFacW","Depth")
    for ii=1:length(list0);
        tmp1=read_bin(grDir*list0[ii]*".data",ioPrec);
        tmp2=Symbol(list0[ii]);
        @eval (($tmp2) = ($tmp1))
        mygrid[list0[ii]]=tmp1
    end

    list0=("DRC","DRF","RC","RF")
    for ii=1:length(list0);
        fil=grDir*list0[ii]*".data"
        tmp1=stat(fil)
        n3=Int64(tmp1.size/8)

        fid = open(fil);
        tmp1 = Array{Float64,1}(undef,n3);
        read!(fid,tmp1);
        tmp1 = hton.(tmp1);

        tmp2=Symbol(list0[ii]);
        @eval (($tmp2) = ($tmp1))
        mygrid[list0[ii]]=tmp1
    end

    return mygrid

end

"""
    GCMGridOnes(grTp,nF,nP)

Define all-1 grid variables instead of using GCMGridSpec & GCMGridLoad.
"""
function GCMGridOnes(grTp,nF,nP)

    global grDir, nFaces, grTopo, ioSize, facesSize, ioPrec;
    global XC, XG, YC, YG, RAC, RAZ, DXC, DXG, DYC, DYG, hFacC, hFacS, hFacW, Depth;

    grDir=""
    grTopo=grTp
    nFaces=nF
    ioSize=[nP nP*nF]
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(nP,nP)]
    ioPrec=Float32

    list0=("XC","XG","YC","YG","RAC","RAZ","DXC","DXG","DYC","DYG","hFacC","hFacS","hFacW","Depth");
    for ii=1:length(list0);
        tmp1=fill(1.,nP,nP*nF);
        tmp1=convert2gcmfaces(tmp1);
        tmp2=Symbol(list0[ii]);
        @eval (($tmp2) = ($tmp1))
    end

    return "GCMGridOnes: passed"

end


function findtiles(ni,nj,grid="llc90")
    mytiles = Dict()
    if grid=="llc90"
        GCMGridLoad()
    else
        println("Unsupported grid option")
    end
    mytiles["nFaces"]=MeshArrays.nFaces;
    #mytiles.fileFormat=mygrid.fileFormat;
    mytiles["ioSize"]=MeshArrays.ioSize;
    %
    XC=MeshArrays.XC;
    YC=MeshArrays.YC;
    XC11=copy(XC); YC11=copy(XC);
    XCNINJ=copy(XC); YCNINJ=copy(XC);
    iTile=copy(XC); jTile=copy(XC); tileNo=copy(XC);
    tileCount=0;
    for iF=1:XC11.nFaces;
        #global tileCount,XC,YC,XC11,YC11,iTile,jTile,tileNo
        face_XC=XC.f[iF]; face_YC=YC.f[iF];
    #ordering convention that was used in first generation nctile files:
    #    for ii=1:size(face_XC,1)/ni;
    #        for jj=1:size(face_XC,2)/nj;
    #ordering convention that is consistent with MITgcm/pkg/exch2:
        for jj=Int.(1:size(face_XC,2)/nj);
            for ii=Int.(1:size(face_XC,1)/ni);
                tileCount=tileCount+1;
                tmp_i=(1:ni).+ni*(ii-1)
                tmp_j=(1:nj).+nj*(jj-1)
                tmp_XC=face_XC[tmp_i,tmp_j]
                tmp_YC=face_YC[tmp_i,tmp_j]
                XC11.f[iF][tmp_i,tmp_j].=tmp_XC[1,1]
                YC11.f[iF][tmp_i,tmp_j].=tmp_YC[1,1]
                XCNINJ.f[iF][tmp_i,tmp_j].=tmp_XC[end,end]
                YCNINJ.f[iF][tmp_i,tmp_j].=tmp_YC[end,end]
                iTile.f[iF][tmp_i,tmp_j]=collect(1:ni)*ones(Int,1,nj)
                jTile.f[iF][tmp_i,tmp_j]=ones(Int,ni,1)*collect(1:nj)'
                tileNo.f[iF][tmp_i,tmp_j]=tileCount*ones(Int,ni,nj)
            end
        end
    end

    mytiles["XC"] = XC;
    mytiles["YC"] = YC;
    mytiles["XC11"] = XC11;
    mytiles["YC11"] = YC11;
    mytiles["XCNINJ"] = XCNINJ;
    mytiles["YCNINJ"] = YCNINJ;
    mytiles["iTile"] = iTile;
    mytiles["jTile"] = jTile;
    mytiles["tileNo"] = tileNo;

    return mytiles

end
