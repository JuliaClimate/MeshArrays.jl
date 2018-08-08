# gcmfaces_convert.jl
#
#	First Draft Implementation
#
# gaelforget (https://github.com/gaelforget/gcmfaces_jl)
# Julia 0.6.2
# Created: 02.01.18
# Last Edit: 02.01.18

## convert2array ##
function convert2array(fld)

if fld.grTopo=="llc";
    tmp1=cat(1,fld.f[1],fld.f[2],rotr90(fld.f[4]),rotr90(fld.f[5]));
    tmp2=cat(1,rotl90(fld.f[3]),NaN*fld.f[3],NaN*fld.f[3],NaN*fld.f[3]);
    arr=cat(2,tmp1,tmp2);
    arr
elseif fld.grTopo=="cs";
    tmp1=cat(1,fld.f[1],fld.f[2],rotr90(fld.f[4]),rotr90(fld.f[5]));
    tmp2=cat(1,rotl90(fld.f[3]),NaN*fld.f[3],NaN*fld.f[3],NaN*fld.f[3]);
    tmp0=cat(1,NaN*fld.f[3],NaN*fld.f[3],NaN*fld.f[3],rotr90(fld.f[6]));
    arr=cat(2,tmp0,tmp1,tmp2);
    arr
elseif fld.grTopo=="ll";
  arr=fld.f[1];
else;
  error("unknown grTopo case");
end;

end
