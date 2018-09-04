
## read_bin function with full list of argument

function read_bin(fil::String,kt,kk,prec::DataType);

  if ~isempty(kt);
    error("non-empty kt option not implemented yet");
  end;

  if ~isempty(kk);
    error("non-empty kk option not implemented yet");
  end;

  (n1,n2)=GCMFaces.ioSize;
  
  if prec==Float64;
    reclen=8;
  else;
    reclen=4;
  end;
  tmp1=stat(fil);
  n3=Int64(tmp1.size/n1/n2/reclen);

  fid = open(fil);
  fld = Array{prec,1}(undef,(n1*n2*n3));
  read!(fid,fld);
  fld = hton.(fld);

  v0=reshape(fld,(n1,n2,n3));
  convert2gcmfaces(v0);

end

## read_bin function with reduced list of argument

function read_bin(fil::String,prec::DataType);
  read_bin(fil,[],[],prec);
end

## read_bin function with reduced list of argument

function read_bin(fil::String);
  read_bin(fil,[],[],Float32);
end
