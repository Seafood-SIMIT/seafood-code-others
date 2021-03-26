% FFT Improvement Factor
% ----------------------

% Chebyshev weights are default.
% To change weights add % to Chebyshev weights 
% and remove % from desired weights

% fft_impfac.m

clear;clc;
j=sqrt(-1);

n=8;             % FFT Points
for K=1:n;       % Fiter Number - 1,2....n
%PRF= 216;       % ARSR-4 Radar PRF - Hz 
PRF=1200;        % ASR-9 Radar PRF - Hz
T=1/PRF;         % Radar Period - s
%sig_f=5.68;     % ARSR-4 Clutter Spectrum- sd -Hz
sig_f=15.3;		  % ASR-9 Clutter Spectrum- sd -Hz


% Find FFT Weights

for k=1:n;
   for i=1:n;
        w(k,i)=exp(-j*2*pi*(i-1)*(k-1)/n);
     end;
end;

ynn=chebwin(n,25);
yn=hamming(n);
wx=ynn'.*w(K,:);	 % Chebyshev Weights
%wx=yn'.*w(K,:);   % Hamming Weights
%wx=w(K,:);        % Uniform Weights

% Clutter Covariance Matrix

x=sig_f/PRF;
omega=2*pi*x;

rho1=exp( -(omega^2)/2);

for ii=1:n;
   for jj=1:n;
      Rn(ii,jj)=rho1^((ii-jj)^2);
   end;
end;

Pn= wx*Rn*wx';

% Signal Covariance Matrix

N=512;
for m=1:N;
   fd=PRF*(m-1)/(N-1)+1e-12;
   xx=2*pi*fd/PRF;
   for nx=1:n;
      for mx=1:n
         Ms(nx,mx)=exp(-j*xx*(mx-nx));
      end;
   end; 
 Ps(K,m)=wx*Ms*wx';     % Transfer Function 
 Imp(K,m)=Ps(K,m)/Pn;         % FFT Improvement Factor
 I_db(K,m)=10*log10(Imp(K,m));
end;
end;

Ps_db=10*log10(Ps);

% Frequency Vector

for l=1:N;
  fd(l)=PRF*(l-1)/(N-1)+1e-12;
 end;
 
figure(1);
%plot(fd,Imp,'k');grid;
%axis([ 0 1200  0 2500]);
plot(fd,I_db);grid;
axis([ 0 1200  0 40]);
title('FFT IMPROVEMENT FACTOR');  
xlabel('Frequency - Hz');
%ylabel('Magnitude');
ylabel('Magnitude - dB');

figure(2);
plot(fd,Ps);grid;
%axis([ 0 1200  0 70]);				% Uniform Weights
axis([ 0 1200  0 40]);           % Chebyshev Weights
%plot(fd,Ps_db,'k');grid;
%axis([ 0 216  -30 20]);
title('FFT TRANSFER FUNCTION');  
xlabel('Frequency - Hz');
ylabel('Magnitude');
%ylabel('Magnitude - dB');
hold on

% Uniform Grid Lines x=0-1200, space=100;y=0-70;space=5
%xl=0;xh=1200;yl=0;yh=70;nx=13;ny=15;
%hx=gridline(xl,xh,nx,yl,yh,ny);

% Chebyshev Grid Lines x=0-1200, space=100;y=0-40;space=5
xl=0;xh=1200;yl=0;yh=40;nx=13;ny=17;
hx=gridline(xl,xh,nx,yl,yh,ny);
hold off


 