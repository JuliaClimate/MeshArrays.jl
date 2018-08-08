
#with default gridName argument:

GCMGridSpec() = GCMGridSpec("LLC90")

#with gridName argument:

function GCMGridSpec(gridName)
    global grDir, nFaces, grTopo, ioSize, facesSize, ioPrec;

if gridName=="LLC90";
    grDir="GRID_LLC90/";
    nFaces=5;
    grTopo="llc";
    ioSize=[90 1170];
    facesSize=[[90 270];[90 270];[90 90];[270 90];[270 90]];
    ioPrec=Float64;
elseif gridName=="CS32";
    grDir="GRID_CS32/";
    nFaces=6;
    grTopo="cs";
    ioSize=[32 192];
    facesSize=[[32 32];[32 32];[32 32];[32 32];[32 32];[32 32]];
    ioPrec=Float32;
elseif gridName=="LL360";
    grDir="GRID_LL360/";
    nFaces=1;
    grTopo="ll";
    ioSize=[360 160];
    facesSize=[360 160];
    ioPrec=Float32;
else;
  error("unknown gridName case");
end;

    "Grid has been defined."

end

function GCMGridLoad()

    global XC, XG, YC, YG, RAC, RAZ, DXC, DXG, DYC, DYG, hFacC, hFacS, hFacW, Depth;

    list0=("XC","XG","YC","YG","RAC","RAZ","DXC","DXG","DYC","DYG","hFacC","hFacS","hFacW","Depth");
    for ii=1:length(list0);
        tmp1=read_bin(grDir*list0[ii]*".data",ioPrec);
        tmp2=Symbol(list0[ii]);
        @eval (($tmp2) = ($tmp1))
    end

    "Grid has been loaded."

end

