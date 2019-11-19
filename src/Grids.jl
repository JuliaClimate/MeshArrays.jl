
## GridSpec function with default gridName argument:

GridSpec() = GridSpec("LLC90")

## GridSpec function with gridName argument:

"""
    GridSpec(gridName)

Return a `gmcgrid` specification that provides grid files `path`,
`class`, `nFaces`, `ioSize`, `facesSize`, `ioPrec`, & a `read` function
(not yet) using hard-coded values for `"LLC90"`, `"CS32"`, `"LL360"` (for now).
"""
function GridSpec(gridName,gridParentDir="./")

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
elseif gridName=="FLTXMPL";
    grDir=gridParentDir*"flt_example/";
    nFaces=4;
    grTopo="dpdo";
    ioSize=[80 42];
    facesSize=[(40, 21), (40, 21), (40, 21), (40, 21)]
    ioPrec=Float32;
else;
    error("unknown gridName case");
end;

mygrid=gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)

return mygrid

end

## GridLoad function

"""
    GridLoad(mygrid::gcmgrid)

Return a `Dict` of grid variables read from files located in `mygrid.path` (see `?GridSpec`).

Based on the MITgcm naming convention, grid variables are:

- XC, XG, YC, YG, AngleCS, AngleSN, hFacC, hFacS, hFacW, Depth.
- RAC, RAW, RAS, RAZ, DXC, DXG, DYC, DYG.
- DRC, DRF, RC, RF (one-dimensional)
"""
function GridLoad(mygrid::gcmgrid)

    GridVariables=Dict()

    list0=("XC","XG","YC","YG","AngleCS","AngleSN","RAC","RAW","RAS","RAZ",
    "DXC","DXG","DYC","DYG","Depth")
    for ii=1:length(list0)
        tmp1=mygrid.read(mygrid.path*list0[ii]*".data",MeshArray(mygrid,mygrid.ioPrec))
        tmp2=Symbol(list0[ii])
        @eval (($tmp2) = ($tmp1))
        GridVariables[list0[ii]]=tmp1
    end

    mygrid.ioPrec==Float64 ? reclen=8 : reclen=4

    list0=("DRC","DRF","RC","RF")
    for ii=1:length(list0)
        fil=mygrid.path*list0[ii]*".data"
        tmp1=stat(fil)
        n3=Int64(tmp1.size/reclen)

        fid = open(fil)
        tmp1 = Array{mygrid.ioPrec,1}(undef,n3)
        read!(fid,tmp1)
        tmp1 = hton.(tmp1)

        tmp2=Symbol(list0[ii])
        @eval (($tmp2) = ($tmp1))
        GridVariables[list0[ii]]=tmp1
    end

    list0=("hFacC","hFacS","hFacW")
    n3=length(GridVariables["RC"])
    for ii=1:length(list0)
        tmp1=mygrid.read(mygrid.path*list0[ii]*".data",MeshArray(mygrid,mygrid.ioPrec,n3))
        tmp2=Symbol(list0[ii])
        @eval (($tmp2) = ($tmp1))
        GridVariables[list0[ii]]=tmp1
    end

    return GridVariables

end

"""
    GridOfOnes(grTp,nF,nP)

Define all-ones grid variables instead of using `GridSpec` & `GridLoad`.
"""
function GridOfOnes(grTp,nF,nP)

    grDir=""
    grTopo=grTp
    nFaces=nF
    if grTopo=="llc"
        ioSize=[nP nP*nF]
    elseif grTopo=="cs"
        ioSize=[nP nP*nF]
    elseif grTopo=="ll"
        ioSize=[nP nP]
    elseif grTopo=="dpdo"
        nFsqrt=Int(sqrt(nF))
        ioSize=[nP*nFsqrt nP*nFsqrt]
    end
    facesSize=Array{NTuple{2, Int},1}(undef,nFaces)
    facesSize[:].=[(nP,nP)]
    ioPrec=Float32

    mygrid=gcmgrid(grDir,grTopo,nFaces,facesSize, ioSize, ioPrec, read, write)

    GridVariables=Dict()
    list0=("XC","XG","YC","YG","RAC","RAZ","DXC","DXG","DYC","DYG","hFacC","hFacS","hFacW","Depth");
    for ii=1:length(list0);
        tmp1=fill(1.,nP,nP*nF);
        tmp1=mygrid.read(tmp1,MeshArray(mygrid,Float64));
        tmp2=Symbol(list0[ii]);
        @eval (($tmp2) = ($tmp1))
        GridVariables[list0[ii]]=tmp1
    end

    return GridVariables

end


"""
    findtiles(ni::Int,nj::Int,mygrid::gcmgrid)
    findtiles(ni::Int,nj::Int,grid::String="llc90",gridParentDir="./")

Return a `MeshArray` map of tile indices for tile size `ni,nj`
"""
function findtiles(ni::Int,nj::Int,mygrid::gcmgrid)
    mytiles = Dict()

    GridVariables=GridLoad(mygrid)

    mytiles["nFaces"]=mygrid.nFaces;
    #mytiles.fileFormat=mygrid.fileFormat;
    mytiles["ioSize"]=mygrid.ioSize;

    XC=GridVariables["XC"];
    YC=GridVariables["YC"];
    XC11=copy(XC); YC11=copy(XC);
    XCNINJ=copy(XC); YCNINJ=copy(XC);
    iTile=copy(XC); jTile=copy(XC); tileNo=copy(XC);
    tileCount=0;
    for iF=1:XC11.grid.nFaces
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

findtiles(ni::Int,nj::Int,grid::String="llc90",gridParentDir="./") = findtiles(ni,nj,GridSpec(grid,gridParentDir))
