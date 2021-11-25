
#bring end points to the equator -> define 3D rotation matrix
function rotate_points(lons,lats)
	#... and of end points
	x0=cosd.(lats).*cosd.(lons)
	y0=cosd.(lats).*sind.(lons)
	z0=sind.(lats)

	#get the rotation matrix:
	#1) rotate around x axis to put first point at z=0
	theta=atan(-z0[1],y0[1])
	R1=[[1;0;0] [0;cos(theta);sin(theta)] [0;-sin(theta);cos(theta)]]
	tmp0=[x0;y0;z0]; tmp1=R1*tmp0; x1=tmp1[1,:]; y1=tmp1[2,:]; z1=tmp1[3,:]
	x0=x1; y0=y1; z0=z1
	#2) rotate around z axis to put first point at y=0
	theta=atan(x0[1],y0[1])
	R2=[[cos(theta);sin(theta);0] [-sin(theta);cos(theta);0] [0;0;1]]
	tmp0=transpose([x0 y0 z0]); tmp1=R2*tmp0; x1=tmp1[1,:]; y1=tmp1[2,:]; z1=tmp1[3,:]
	x0=x1; y0=y1; z0=z1
	#3) rotate around y axis to put second point at z=0
	theta=atan(-z0[2],-x0[2])
	R3=[[cos(theta);0;-sin(theta)] [0;1;0] [sin(theta);0;cos(theta)]]
	tmp0=transpose([x0 y0 z0]); tmp1=R3*tmp0; x1=tmp1[1,:]; y1=tmp1[2,:]; z1=tmp1[3,:]
	x0=x1; y0=y1; z0=z1

	x0,y0,z0,R3*R2*R1
end

function rotate_XCYC(Γ,R)
	#3D carthesian coordinates:
	lon=Γ.XC; lat=Γ.YC
	x=cosd.(lat)*cosd.(lon)
	y=cosd.(lat)*sind.(lon)
	z=sind.(lat)

	#rotation R:
	tmpx=γ.write(x); tmpy=γ.write(y); tmpz=γ.write(z)
	tmp1=findall((!isnan).(tmpx))
	tmpx2=tmpx[tmp1]; tmpy2=tmpy[tmp1]; tmpz2=tmpz[tmp1]
	tmp2=[tmpx2';tmpy2';tmpz2']
	
	tmp3=R*tmp2

	tmpx2=tmp3[1,:]; tmpy2=tmp3[2,:]; tmpz2=tmp3[3,:]
	tmpx[tmp1]=tmpx2; tmpy[tmp1]=tmpy2; tmpz[tmp1]=tmpz2
	x=γ.read(tmpx,Γ.XC); y=γ.read(tmpy,Γ.XC); z=γ.read(tmpz,Γ.XC)
	
	x,y,z
end

##

function shorter_paths!(xyz,xyz0,msk_in)
	(x0,y0,z0)=xyz0[:]
	(x,y,z)=xyz[:]

	#split in two segments:
	theta=zeros(2)
	theta[1]=atan(y0[1],x0[1])
	theta[2]=atan(y0[2],x0[2])

	tmpx=γ.write(x); tmpy=γ.write(y); tmpz=γ.write(z);
	tmptheta=atan.(tmpy,tmpx)
	if theta[2]<0;
		tmp00=findall(tmptheta.<=theta[2])
		tmptheta[tmp00].=tmptheta[tmp00]+2*pi
		theta[2]=theta[2]+2*pi
	end

	msk_out=[]
	for kk in 1:3
		#select field to treat:
		mm=msk_in[kk]
		#select the shorther segment:
		tmpm=γ.write(mm)
       if theta[2]-theta[1]<=pi
            tmpm[findall( (tmptheta.>theta[2]).|(tmptheta.<theta[1]) )].=0.0
        else
            tmpm[findall( (tmptheta.<=theta[2]).&(tmptheta.>=theta[1]) )].=0.0
        end
        mm=γ.read(tmpm,mm);
        #store result:
		push!(msk_out,mm)
    end
	msk_out[:]
end

