function [Sw,w]=specden(Xt,dt,Navg) 
%[Sw,w]=specden(Xt,dt) 
% gives the two sided spectral density of time series Xt of length Nt sampled
% at time interval dt. Only the positive freq. portion is contained in Sw.
% The length of Sw is Nt/2+1 when Nt is even and (Nt+1)/2 when Nt is odd, 
% as the symmetric part of Sw in the original fft is omitted. w is a vector 
% containing the corresponding frequencies.
%
% Mathematically, if Xt is one-dimensional
%  Sw = 1/T*abs(fft(Xt)*dt).^2 
% 
% If Xt is multi-dimensional, the cross spectrum is returned
% in which case Sw is a Nx*Nx cell of column vectors
%  Sw{i,j}=1/T*fft(Xt(:,i))*conj(fft(Xt(:,j)))
%       
% Note that theoretically and by construction,
%  Sw(i,j)=conj(Sw(j,i)), dw = (2*pi)/T
%
% where T=(Nt-1)*dt is the duration of the time series.
%
% By construction, 2*sum(Sw)*dw = std(Xt)^2 (Parseval Equality)
%
% ... = specden(Xt,dt,Navg) computes the averaged spectral density
% by dividing the data into Navg time segments. In this case, individual
% segments are detrended before FFT is calculated.
%
% Example: 
%   [Sw,w]=specden(Xt,dt); [r,tau]=spec2cov(Sw,dw)
% Then the following holds
%       2*sum(Sw)*dw = std(Xt)^2 = r(1)
%
% See also SPEC2COV.

% 150330, remove factor 2*pi in (2*pi*T)
%         revised pltspec, pltband accordingly
% Revision 1.4, 7/11/08, CityU, allow Navg
% Revision 1.3, 25/10/08, CityU, change Sw as a 3-D array
% Revision 1.2, 8/10/97, Caltech, documentation
% Revision 1.1, 27/3/97, HKUST
% Written by Sui-Kui Au, 3/27/97, HKUST
% 091010 - modified to save space

if nargin>2
    Sw = 0;
    nt = floor(size(Xt,1)/Navg);
    for k = 1:Navg
      I = (k-1)*nt + [1:nt];
      if k==1
        [Sw,w] = specden(detrend(Xt(I,:)),dt);
      else
        Sw = Sw + specden(detrend(Xt(I,:)),dt);
      end
    end
    Sw = Sw/Navg;
    return;
end

[Nt,Nx]=size(Xt);
Nw=floor((Nt+2)/2);
% T=(Nt-1)*dt;
T=Nt*dt;
dw=2*pi/Nt/dt; % note that the lowest nonzero freq. is not 2*pi/(Nt-1)/dt

w = [0:dw:pi/dt];
ncol = size(Xt,2);
if ncol==1
    w=w.';
end    

%Sw=1/(2*pi*dt*Nt)*abs(fft(Xt)*dt).^2;
fftX=fft(Xt)*dt;	% process by column
clear Xt; % to save space

Sw = zeros(Nw,Nx,Nx);
for ix=1:Nx
  for jx=ix:Nx
    if Nx==1
     Sw=1/(dt*Nt)*fftX(1:Nw,ix).*conj(fftX(1:Nw,jx));
    else
     Sw(:,ix,jx) = 1/(dt*Nt)*fftX(1:Nw,ix).*conj(fftX(1:Nw,jx));
     Sw(:,jx,ix) = conj(Sw(:,ix,jx));
    end              
    if ix==jx
      Sw(:,ix,ix)=real(Sw(:,ix,ix)); % remove complex due to round-off
    end
  end
end

