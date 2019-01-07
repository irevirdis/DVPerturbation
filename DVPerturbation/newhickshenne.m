#!/usr/bin/octave

pkg load symbolic

syms s
%SC=load('SC');
SC=load('bump_new.txt');
AIR=load('blade_from_le.txt');

N=length(AIR);
sc=length(SC);
t2=6.65;

tmp=[AIR zeros(N,sc+7)];

for i=2:N
	tmp(i,4)=sqrt((tmp(i,2)-tmp(i-1,2))^2+(tmp(i,3)-tmp(i-1,3))^2); # DISTANCE BETEWEEN NEIGHBOURING POINTS ALONG THE SURFACE #
	tmp(i,5)=tmp(i-1,5)+tmp(i,4); # CUMULATIVE CURVILINEAR COORDINATE #
end

tmp(:,6)=tmp(:,5)/tmp(N,5); # NON-DIMENSIONAL CUMULATIVE CURVILINEAR COORDINATE #


for j=1:sc
	jmin=j-1;
	jmax=j+1;
	if (j==1) jmin=sc; end
	if (j==sc) jmax=1; end
	scenter=SC(j,1);
	smin=SC(jmin,1);
	smax=SC(jmax,1);
	amp=SC(j,2)
	diffmin=scenter-smin;
	diffmax=smax-scenter;
	if (diffmin<0) 
		diffmin=diffmin+1;
		smin=smin-1;
        end
	if (diffmax<0) 
		diffmax=diffmax+1;
		smax=smax+1;
	end
	for i=1:N
		s=tmp(i,6);
		if (s-smin>1.0) s=s-1.0; end
		if (smax-s>1.0) s=s+1.0; end
	 	if s>=smin && s<scenter;
			tmp(i,6+j)=amp*(sin(0.5*pi*(s-smin)/diffmin))^t2; # FIRST HALF OF HH BUMP FUNCTION ON CONTROL POINT j STARTING FROM ZERO ON POINT j-1 #
		elseif s>=scenter && s<=smax;
	  		tmp(i,6+j)=amp*(sin(0.5*pi*(smax-s)/diffmax))^t2; # SECOND HALF OF HH BUMP FUNCTION ON CONTROL POINT j REACHING ZERO ON POINT j+1 #
		else
			tmp(i,6+j)=0.0;
		end
	end
end


for j=1:sc
	for i=1:N
		im=i-1;
		ip=i+1;
		if (i==1) im=N; end
		if (i==N) ip=1; end
		xb=tmp(im,2); xm=tmp(i,2); xf=tmp(ip,2);
		yb=tmp(im,3); ym=tmp(i,3); yf=tmp(ip,3);
		DY=yf-yb; DX=xf-xb;
		M=-DX/DY;
		a=atan(M);
		d=tmp(i,6+j);
		if DX>0 && DY<0; # PRIMO QUADRANTE #
			tmp(i,7+sc)=tmp(i,7+sc)+sign(SC(j,2))*(abs(d*cos(atan(M)))); # NEW X VALUE FOR INTERNAL SURFACE NODE #
			tmp(i,8+sc)=tmp(i,8+sc)+sign(SC(j,2))*(abs(d*sin(atan(M)))); # NEW Y VALUE FOR INTERNAL SURFACE NODE #
		elseif DX<0 && DY<0; # SECONDO QUADRANTE #
			tmp(i,7+sc)=tmp(i,7+sc)+sign(SC(j,2))*(abs(d*cos(atan(M)))); # NEW X VALUE FOR INTERNAL SURFACE NODE #
	                tmp(i,8+sc)=tmp(i,8+sc)-sign(SC(j,2))*(abs(d*sin(atan(M)))); # NEW Y VALUE FOR INTERNAL SURFACE NODE #
		elseif DX<0 && DY>0; # TERZO QUADRANTE #
	                tmp(i,7+sc)=tmp(i,7+sc)-sign(SC(j,2))*(abs(d*cos(atan(M)))); # NEW X VALUE FOR INTERNAL SURFACE NODE #
	                tmp(i,8+sc)=tmp(i,8+sc)-sign(SC(j,2))*(abs(d*sin(atan(M)))); # NEW Y VALUE FOR INTERNAL SURFACE NODE #
		elseif DX>0 && DY>0; # QUARTO QUADRANTE #
			tmp(i,7+sc)=tmp(i,7+sc)-sign(SC(j,2))*(abs(d*cos(atan(M)))); # NEW X VALUE FOR INTERNAL SURFACE NODE #
	                tmp(i,8+sc)=tmp(i,8+sc)+sign(SC(j,2))*(abs(d*sin(atan(M)))); # NEW Y VALUE FOR INTERNAL SURFACE NODE #
		end
	end
end

tmp(:,9+sc)=tmp(:,7+sc)+tmp(:,2);
tmp(:,10+sc)=tmp(:,8+sc)+tmp(:,3);

sp=[(tmp(:,1)) tmp(:,9+sc) tmp(:,10+sc)];

MAT=[tmp(:,2:3) tmp(:,1) tmp(:,5) tmp(:,6) sqrt((tmp(:,7+sc)).^2+(tmp(8+sc)).^2)];

id=fopen('surface_positions.dat','w');
for i=1:N
	fprintf(id,'%d %10.6f %10.6f \n', sp(i,:));
end
save ('-ascii','TMP.dat','tmp')
%save ('-ascii','-double','surface_positions.dat','sp')
