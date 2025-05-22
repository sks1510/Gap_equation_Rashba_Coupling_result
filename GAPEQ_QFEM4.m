% This code solves the gap equation matrix and search temperature (T_c) where he eigenvalue is 1.

format("long");

tic


% initializing squared frequency of insulator, omega_0^2
om0sqin=-0.0;
om0sq=0.0;        

%pl=parpool(18);

for ilpomg=1:50      % loop for squared frequency of insulator omega_0^2

om0sq=om0sqin-ilpomg*0.2;


fileID=fopen(['eigendata_loop1_m' num2str(abs(om0sq)) '.dat'],'w');    % opening datafile to write
fprintf(fileID,'%1s %12s    %8s \n', '#',  'ne*10^18','Tc(K)');



%setting the momentum integration and frequency sum
nsm=200;
lambda=1*pi;  % upper momentum cutoff
M=100;  % Matsubara frequencies


% initializing carrier density 
nex=0;
delnex=1;


for ine=1:1000       % carrier density loop

nex=nex+delnex;
ned=nex*10^(18);    %carrier density in cm^-3

disp(nex);

%parameters list1
alt=3.9;    % lattice constant in Angstrom
cbs=7.5*10^5*10^8;    %phonon velocity in  Angstrom/s
hcut=6.582119*10^(-13);   %Planck constant in meV-s
kboltz=8.617333*10^(-5);    %Boltzmann constant in eV/K
melec=0.51099*10^6/(3*10^10)^2;      %mass of electron in eV/(cm/s)^2
kfa=(3*pi^2*ned)^(1/3)*alt*10^(-8);  % density

% defining linear bare coupling constant, gto
gfac=1.44;     %Angstrom
t4=27;         %meV
xi=28.0/3.0;     %meV
tau=115;          %meV/Angstrom
gto=tau*sqrt(2)*sin(kfa/sqrt(2));
gto=gto/sqrt(1+(4*t4*(sin(kfa/sqrt(2)))^2/xi)^2);
gto=gto*gfac;        %meV


% DOS of STO
ndos=0.524805*(kfa^1.33358)*10^(-3);  %meV^-1
melecef=(ndos*2*pi^2*hcut^2)*10^(-3)/((3*pi^2*ned)^(1/3)*(alt*10^(-8))^3); %eV/(cm/s)^2
mratio=melecef/melec;  %effective mass

% defining phonon frequency at finite carrier density, 
omegatogsq=(om0sq+20.837*(ned/10^(20)));    %meV^2
omegatogsqnl=omegatogsq;

%  phonon frequency in the ordered state
if omegatogsq < 0
omegatogsq=2*abs(omegatogsq);    % phonon frequency in the FE state in meV
end
omegatog=sqrt(omegatogsq);   % phonon frequency in the PE state in meV


% dimensionless "r" in the bosonic propagator
ar1=omegatog^2*(alt/(hcut*cbs))^2;
rc=ar1;



%% linear effective Coupling Constant
ggm=gto^2*ndos;
ggm=ggm/1.2392;    % rescaling the Coupling Constant


%% Truncated linear bare Coupling constant
%alpha=xxx;  %parameters
%beta=xxx;    %parameters
%xne=xxx;     %parameters
%x0=xxx;      %parameters
%y0=xxx;      %parameters
%gtomd=gto*(1/(1+exp((xne-alpha*x0)/beta)));
%ggm=gtomd^2*ndos;
%ggm=ggm/1.2392;




%%  nonlinear coupling
ggm2=0.0;  
%omegaL=xxx;    %meV  (parameter sqrt(b))
%ar2=0.0d0;
%if om0sq < 0
%aneqcp=-om0sq/20.837d0;
%if (ned/10^(20)) < aneqcp
%ar2=(-omegatogsqnl/omegaL^2);
%end
%end

%gto2=xxxx;     %meV; bare nonlinear coupling 
%ggm2=gto2^2*ndos*ar2;
%ggm2=ggm2/1.2392;   %rescaling of coupling constant





%parameters list2
kfermi=(3*pi^2*ned)^(1/3);  %Fermi momentum in cm^-1
efermi=hcut^2*10^(-6)*kfermi^2/(2.0*melecef); %Fermi energy in eV
efermi_red=efermi/kboltz;   %Fermi energy in K
vf=sqrt(2.0*efermi/melecef);   %Fermi velocity in cm/s
kavf=(vf/alt)*10^(8);          % v_F/a in 1/s
kbbyhct=(kboltz/hcut)*10^3;    % K_B/hcut in 1/(K-s)


%effective coupling constant in gap equation
fac=ggm+ggm2;
facpi=ggm*8*kfa^2;



% initializing temperature
temp=0.0;      % K
  nt=10^10;
  delnt=0.000001;
ncnt=0;
ncnt1=0;
ncnt2=0;
ncnt3=0;
ncnt4=0;

  for ll=1:nt   % temperature Loop

%  reducing number of temperature iterations and simultaneously increasing the precision of the T_c value	  
	  
if ncnt==0
    delnt=delnt*10;
end
      
      if ncnt==1
          if ncnt1==0
      delnt = delnt/10.0;
    temp=0;
    ncnt1=ncnt1+1;
          end
     end
if ncnt==2
    if ncnt2==0
    temp=temp-delnt;
    delnt=delnt/10; 
    ncnt2=ncnt2+1;
    end
end
if ncnt==3
    if ncnt3==0
    temp=temp-delnt;
    delnt=delnt/10;
    ncnt3=ncnt3+1;
    end
end

if ncnt==4
    if ncnt4==0
    temp=temp-delnt;
    delnt=delnt/10;
    ncnt4=ncnt4+1;
    end
end

if ncnt==5
xfrmt=' %10.5f  %18.13f \n ';
fprintf(fileID,xfrmt,ned*10^(-18),temp-delnt);    % writing T_c vs. carrier density in datafile
break
end

tempold=temp;
  
if ncnt==0
    temp=delnt;
else
temp=temp+delnt;
end


% frequency rescaling

   fener=efermi_red;
mxkl=fix(0.5*(efermi_red/(pi*temp)-1));
ds=0.1;
if mxkl>M
   imll=fix(0.1*M);
        Mq=M;
       dimq=M+1;
aal=(log(mxkl-Mq)-log(ds))/((Mq-imll)*1.0);
abl=(imll*log(mxkl-Mq)-Mq*log(ds))/((Mq-imll)*1.0);
else
  imll=mxkl+1;
        Mq=mxkl;
       dimq=mxkl+1;
        aal=0;    
        abl=0;    
end



   % Initialize Matrix
   mat=zeros(dimq,dimq);


% Matrix formation
for ii=1:dimq      % frequency loop

 nopii=0;
 if ii >= (1+imll)
    nopii=1;
end
fk=(2*(ii-1)+1)*pi*temp+nopii*2*exp(aal*(ii-1)-abl)*pi*temp;    

sumq0=0.0;
for iq=1:dimq      % frequency loop
 nopiq=0;

if iq >= (1+imll)
    nopiq=1;
end

fq=(2*(iq-1)+1)*pi*temp+nopiq*2*exp(aal*(iq-1)-abl)*pi*temp;


if ne(ii,iq)==1
dfnqp=0.0;
dfnqm=0.0;
dh=1.0/nsm;
qlm=0.000001;

for irr=1:nsm+1  % momentum integration of bosonic function

bpropm=qlm/(rc*lambda^(-2)+qlm^2+lambda^(-2)*(fq-fk)^2*(kbbyhct*alt/cbs)^2+lambda^(-3)*facpi*(abs(fq-fk)/qlm)*(kbbyhct*(1/2)/kavf)*atan(lambda*(kavf/kbbyhct)*qlm/abs(fq-fk)));

bpropp=qlm/(rc*lambda^(-2)+qlm^2+lambda^(-2)*(fq+fk)^2*(kbbyhct*alt/cbs)^2+lambda^(-3)*facpi*(abs(fq+fk)/qlm)*(kbbyhct*(1/2)/kavf)*atan(lambda*(kavf/kbbyhct)*qlm/abs(fq+fk)));


if irr==1
    dfnqp=dfnqp+bpropp;
    dfnqm=dfnqm+bpropm;
elseif irr==nsm+1
    dfnqp=dfnqp+bpropp;
    dfnqm=dfnqm+bpropm;
else
if mod(irr,2)==0
    dfnqp=dfnqp+bpropp*4.0;
    dfnqm=dfnqm+bpropm*4.0;
else
    dfnqp=dfnqp+bpropp*2.0;
    dfnqm=dfnqm+bpropm*2.0;
end
end
qlm=qlm+dh;
end    % momentum  integration ends

if  iq < (imll+1)
sumq0=sumq0+(-dfnqp+dfnqm)*dh/3.0;
else
sumq0=sumq0+((-dfnqp+dfnqm)*dh/3.0)*(1+exp(aal*(iq-1)-abl)*(exp(aal)-1));
end

end

end
sumq0=1+fac*temp*sumq0/fk;        %denominator in the gap equation    


for j=1:dimq       % frequency loop


nopj=0;
if j >= (1+imll)
    nopj=1;
end
fp=(2*(j-1)+1)*pi*temp+nopj*2*exp(aal*(j-1)-abl)*pi*temp;

dfnp=0.0;
dfnm=0.0;
dh=1.0/nsm;
plm=0.000001;

for ir=1:nsm+1  % momentum integration bosonic function

bpropm=plm/(rc*lambda^(-2)+plm^2+lambda^(-2)*(fp-fk)^2*(kbbyhct*alt/cbs)^2+lambda^(-3)*facpi*(abs(fp-fk)/plm)*(kbbyhct*(1/2)/kavf)*atan(lambda*(kavf/kbbyhct)*plm/abs(fp-fk)));

bpropp=plm/(rc*lambda^(-2)+plm^2+lambda^(-2)*(fp+fk)^2*(kbbyhct*alt/cbs)^2+lambda^(-3)*facpi*(abs(fp+fk)/plm)*(kbbyhct*(1/2)/kavf)*atan(lambda*(kavf/kbbyhct)*plm/abs(fp+fk)));


if ir==1
    dfnp=dfnp+bpropp;
    dfnm=dfnm+bpropm;
elseif ir==nsm+1
    dfnp=dfnp+bpropp;
    dfnm=dfnm+bpropm;
else
if mod(ir,2)==0
    dfnp=dfnp+bpropp*4.0;
    dfnm=dfnm+bpropm*4.0;
else
    dfnp=dfnp+bpropp*2.0;
    dfnm=dfnm+bpropm*2.0;
end
end
plm=plm+dh;
end    % momentum integration ends

dfnp=dfnp*dh/3.0;
dfnm=dfnm*dh/3.0;


% matrix elements
matelem=0.0;
if ne(ii,j)==1
matelem=matelem+((dfnm+dfnp)/fp)/sumq0;
end

      if  j < (imll+1)
               mat(ii,j)=fac*matelem*temp;
              else
              mat(ii,j)=fac*matelem*temp*(1+exp(aal*(j-1)-abl)*(exp(aal)-1));    %matrix(k0 x p0) elements
       end


        end
        end
%mat;
[V,D]=eigs(mat,1);
Y=[temp,D];
if D<1
    ncnt=ncnt+1;
end
%disp(Y);
    

  
  
  end     %  temperature loop end
 
  
end   % density loop end  
  
  
  fclose(fileID);

end

%delete(pl)

toc
