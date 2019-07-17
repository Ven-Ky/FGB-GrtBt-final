%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						
% Model for diffusional exchange of Fe-Mg between two minerals Grt and
% Bt (in which diffusion is assumed to be sufficiently rapid to allow 
% the Btal to continuously homogenize).  Minerals are approximated as 
% spheres, with isotropic diffusion, and equilibrium is maintained with 
% the Btal at the surfaces of the grains.
%
% m --> g.b.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;	 %XCHG2P_en_ab_Si_isotopes

tic

P = 1.2e9;                      % pressure (Pa)
Pkbar = P/1e8;                  % presure (Kbar)
To = 800+273;					% initial temperature (K)
nr = 300;						% number of radius steps
dTdt = 10/3.1536e13;            % Cooling rate (K/s); (3.1536e13 seconds = 1 Myr)
Tfin = 550+273;                 % Final temperature
tmax = (To-Tfin)/dTdt;          % total time (s)
nt = 100;						% # timesteps
deltat = tmax/nt;				% time step (s)
t = [0:deltat:(nt-1)*deltat];	% time vector 
size(t)

T = zeros(1,nt);
T(1:nt) = [To:-dTdt*deltat:(Tfin+dTdt*deltat)];

%%%%%%% Beta factors for Fe56/Mg24 in Gt and Bt %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                             

%Beta_Gt=2104./T.^2;         %orthopyroxene; simplified from Wang et al. (2017)
%Beta_Bt=1458./T.^2;          %clinopyroxene; simplified from Wang et al. (2017)
%Beta_gb=1800./T.^2;          %completely made up (and doesn't matter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %alphaGtgb = exp((Beta_Gt-Beta_gb));
 %alphaBtgb = exp((Beta_Bt-Beta_gb));

 %a = 2.7635e6;
 %b = -1.1561e4;
 %c = 3.8431;
 
 %%%%%%% Initial composition of Grt using Bt ,Kd,Mass_balance_eq%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XmgBt = 0.6470;                    %Initial Bt mg molefractions
XfeBt = 1-XmgBt;                   %Initial Bt Fe molefractions
fgb = 1e-6;                      % Mass prop of total Bt in rock BtF
fgrt = 0.5; %Mass prop of total Grt in rock
fBt = 1-fgrt-fgb;
lnKd1 = (19.506).*(T)-52107.536
lnkd2 = 3*8.314.*(T); % Kd at 800C intital
Kd = exp(lnKd1./lnkd2);
R = Kd(1).*(XmgBt/XfeBt); %Ratio of Fe/Mg in Grt 
XmgGrt =R/(R+1);        % Through mass balance and Xfebulk+Xmgbulk = 1
XfeGrt = 1-XmgGrt; 
Xmgbulk = XmgBt*fBt+XmgGrt*fgrt; % Bulk concentration of the system 
Xfebulk = 1-Xmgbulk;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% calculation of Kmg and KFe partition coffecients %%%%%%%%%%%%%%
%solving quadratic equation from mass balance%
asolve = Kd*(fgrt)^2;
bsolve = (fgrt)*(fBt*(1+Kd)- (Xmgbulk + (Kd*(Xfebulk))));
csolve = fBt*(fBt-1);
Kfe = (-bsolve +sqrt( (bsolve).^2- 4*asolve*csolve))./(2.*asolve);  % equilibrium Grt/Bt partition coefficient
Kmg = Kd.*Kfe; % equilibrium Grt/Bt partition coefficient
% Nur1 = ((Xfebulk/XfeBt)-fBt)/fgrt;
% Kmg =  (Nur1).*Kd;
% Nur2 = ((Xmgbulk/XmgBt)-fBt)/fgrt;
% Kfe = Nur2./Kd;    
kMg24Bt = ones(size(t));            % equilibrium Bt/gb partition coefficient
kFe56Bt = ones(size(t));
kMg24Gt = Kmg;       % equilibrium Gt/gb partition coefficient
kFe56Gt = Kfe;					

mMg24 = 24;                         % atomic mass of Mg24
mFe56 = 56;                         % atomic mass of Fe56
lambdaMg24 = 0;                     % decay constant of Mg24
lambdaFe56 = 0;                  

BtF = fBt;          % volume fraction Bt (= Bt/(g.b.+Gt+Bt)
mF = fgb;                      % gb fraction (= g.b./(g.b.+Gt+Bt))
GtF = fgrt;               % volume fraction gt

RBt = 3e-4;					% Bt grain radius (m)
RGt = 3e-4;					% Gt grain radius (m)
V = (4/3)*pi*RBt^3/BtF;		% system volume (m^3) (system defined as containing
								% one Bt sphere)
VBt = V*BtF;                  % Bt volume (m^3)
VGt = V*GtF;                  % Gt volume
Vm = V*mF;
NGt = (GtF*RBt^3)/(BtF*RGt^3);	% # of Gt spheres

CoMg24 = Xmgbulk;
CoMg24m = CoMg24/(mF+BtF+GtF*kMg24Gt(1));
CoMg24Bt =CoMg24m;        %error
CoMg24Gt =CoMg24Bt*kMg24Gt(1);         %error  

CoFe56 = Xfebulk;
CoFe56m = CoFe56/(mF+BtF+GtF*kFe56Gt(1));
CoFe56Bt = CoFe56m;     %error
CoFe56Gt = CoFe56m*kFe56Gt(1);      %error

DoFegt = 1.64e-10;
DoMggt = 2.72e-10;
QFegt = 226900;
QMggt = 228300;
VFegt = 5.6e-6;
VMggt = 5.3e-6;
DFegt(1) = DoFegt*exp(-(QFegt+P*VFegt)/(8.3145*T(1)));
DMggt(1) = DoMggt*exp(-(QMggt+P*VMggt)/(8.3145*T(1)));
xFe_gt = 0.5;

DoFe56Bt = DoFegt*10^3;		% pre-exponential term for Fe-Mg diffusion in Bt Mueller et al (2013)
VFe56Bt = 5.6e-6;            % activation volume (m^3/mol)
QFe56Bt = 226900;			% activation energy (J/mol) (fitting low T ferromagnetic data from Lubbenhusen and Mehrer 1990)
DFe56Bt(1) = DoFe56Bt*exp(-(QFe56Bt+P*VFe56Bt)/(8.3145*To));
%fK_Bt = 0.5;               % correlation coefficient * coupling constant
DMg24Bt(1) = DFe56Bt(1);

%DoFe56Gt = 6.93e-6;		% pre-exponential term for Ca diffusion in gt (Borinski et al., 2012)
%VFe56Gt = 5e-6;            % activation volume (m^3/mol)
%QFe56Gt = 384000;		% activation energy (J/mol) 
DFe56Gt(1) = DFegt(1)*DMggt(1)/(xFe_gt*DFegt(1)+(1-xFe_gt)*DMggt(1));
%fK_Gt = 0.3;               % correlation coefficient * coupling constant
DMg24Gt(1) = DFe56Gt(1);

massBtFe56 = zeros(1,nt);
massGtFe56 = zeros(1,nt);
massBtMg24 = zeros(1,nt);
massGtMg24 = zeros(1,nt);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 					    INITIAL CONDITIONS						
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
massBtFe56(1) = CoFe56Bt*4/3*pi*(RBt^3);
massGtFe56(1) = NGt*CoFe56Gt*4/3*pi*(RGt^3);
massBtMg24(1) = CoMg24Bt*4/3*pi*(RBt^3);
massGtMg24(1) = NGt*CoMg24Gt*4/3*pi*(RGt^3);
massmFe56(1) = CoFe56m*Vm;
massmMg24(1) = CoMg24m*Vm;

CsMg24(1:nr,1) = CoMg24Bt.*ones(nr,1);
CsFe56(1:nr,1) = CoFe56Bt.*ones(nr,1);	
CsFe56(nr,1) = CoFe56m*kFe56Bt(1);
CsMg24(nr+1:2*nr,1) = CoMg24Gt.*ones(nr,1);
CsFe56(nr+1:2*nr,1) = CoFe56Gt.*ones(nr,1);
CsFe56(2*nr,1) = CoFe56m*kFe56Gt(1);

rBt(1:nr,1) = [0:RBt/(nr-1):RBt]';
rGt(1:nr,1) = [0:RGt/(nr-1):RGt]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%					   BEGIN TIME STEPPING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 2:nt

	DFe56Bt(j) = DoFe56Bt*exp(-(QFe56Bt+P*VFe56Bt)/(8.3145*T(j)));
    DMg24Bt(j) = DFe56Bt(j);
    DFegt(j) = DoFegt*exp(-(QFegt+P*VFegt)/(8.3145*T(j)));
    DMggt(j) = DoMggt*exp(-(QMggt+P*VMggt)/(8.3145*T(j)));
    DFe56Gt(j) = DFegt(j)*DMggt(j)/(xFe_gt*DFegt(j)+(1-xFe_gt)*DMggt(j));
    DMg24Gt(j) = DFe56Gt(j);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ve = zeros(size(t));					% melt extraction rate vector

    centerAMg24Bt = (1+(DMg24Bt(j)/(RBt/nr)^2+0.5*lambdaMg24)*deltat).*ones(nr,1);
	centerAMg24Bt(1) = 1+(3*DMg24Bt(j)/(RBt/nr)^2+0.5*lambdaMg24)*deltat;
	AdMg24Bt = diag(centerAMg24Bt);
	centerAMg24Gt = (1+(DMg24Gt(j)/(RGt/nr)^2+0.5*lambdaMg24)*deltat).*ones(nr,1);
	centerAMg24Gt(1) = 1+(3*DMg24Gt(j)/(RGt/nr)^2+0.5*lambdaMg24)*deltat;
	AdMg24Gt = diag(centerAMg24Gt);
	centerBMg24Bt = (1-(DMg24Bt(j-1)/(RBt/nr)^2+0.5*lambdaMg24)*deltat).*ones(nr,1);
	centerBMg24Bt(1) = 1-(3*DMg24Bt(j-1)/(RBt/nr)^2+0.5*lambdaMg24)*deltat;
	BdMg24Bt = diag(centerBMg24Bt);
	centerBMg24Gt = (1-(DMg24Gt(j-1)/(RGt/nr)^2+0.5*lambdaMg24)*deltat).*ones(nr,1);
	centerBMg24Gt(1) = 1-(3*DMg24Gt(j-1)/(RGt/nr)^2+0.5*lambdaMg24)*deltat;
	BdMg24Gt = diag(centerBMg24Gt);
	upperAMg24Bt(1) = -3*DMg24Bt(j)*deltat/(RBt/nr)^2;
	lowerAMg24Bt(1) = 0;
	upperBMg24Bt(1) = 3*DMg24Bt(j-1)*deltat/(RBt/nr)^2;
	lowerBMg24Bt(1) = 0;
	upperAMg24Gt(1) = -3*DMg24Gt(j)*deltat/(RGt/nr)^2;
	lowerAMg24Gt(1) = 0;
	upperBMg24Gt(1) = 3*DMg24Gt(j-1)*deltat/(RGt/nr)^2;
	lowerBMg24Gt(1) = 0;
	for i = 2:nr-1;
		upperAMg24Bt(i) = -(i)*DMg24Bt(j)*deltat/(2*(i-1)*(RBt/nr)^2);
		lowerAMg24Bt(i) = -(i-1)*DMg24Bt(j)*deltat/(2*i*(RBt/nr)^2);
		upperBMg24Bt(i) = (i)*DMg24Bt(j-1)*deltat/(2*(i-1)*(RBt/nr)^2);
		lowerBMg24Bt(i) = (i-1)*DMg24Bt(j-1)*deltat/(2*i*(RBt/nr)^2);
		upperAMg24Gt(i) = -(i)*DMg24Gt(j)*deltat/(2*(i-1)*(RGt/nr)^2);
		lowerAMg24Gt(i) = -(i-1)*DMg24Gt(j)*deltat/(2*i*(RGt/nr)^2);
		upperBMg24Gt(i) = (i)*DMg24Gt(j-1)*deltat/(2*(i-1)*(RGt/nr)^2);
		lowerBMg24Gt(i) = (i-1)*DMg24Gt(j-1)*deltat/(2*i*(RGt/nr)^2);
	end
	AudMg24Bt = diag(upperAMg24Bt,1);
	AldMg24Bt = diag(lowerAMg24Bt,-1);
	BudMg24Bt = diag(upperBMg24Bt,1);
	BldMg24Bt = diag(lowerBMg24Bt,-1); 
	AudMg24Gt = diag(upperAMg24Gt,1);
	AldMg24Gt = diag(lowerAMg24Gt,-1);
	BudMg24Gt = diag(upperBMg24Gt,1);
	BldMg24Gt = diag(lowerBMg24Gt,-1);
	AMg24Bt = AdMg24Bt + AudMg24Bt + AldMg24Bt;
	BMg24Bt = BdMg24Bt + BudMg24Bt + BldMg24Bt;
	AMg24Gt = AdMg24Gt + AudMg24Gt + AldMg24Gt;
	BMg24Gt = BdMg24Gt + BudMg24Gt + BldMg24Gt;
	bottomAMg24Bt(1) = (1+0.5*lambdaMg24*deltat)*0.206*4/3*pi*(rBt(2)^3-rBt(1)^3);
	bottomBMg24Bt(1) = (1-0.5*lambdaMg24*deltat)*0.206*4/3*pi*(rBt(2)^3-rBt(1)^3);
	bottomAMg24Bt(nr) = (1+0.5*lambdaMg24*deltat)*0.794*4/3*pi*(rBt(nr)^3-rBt(nr-1)^3);
	bottomBMg24Bt(nr) = (1-0.5*lambdaMg24*deltat)*0.794*4/3*pi*(rBt(nr)^3-rBt(nr-1)^3);
	bottomAMg24Gt(1) = (1+0.5*lambdaMg24*deltat)*NGt*0.206*4/3*pi*(rGt(2)^3-rGt(1)^3);
	bottomBMg24Gt(1) = (1-0.5*lambdaMg24*deltat)*NGt*0.206*4/3*pi*(rGt(2)^3-rGt(1)^3);
	bottomAMg24Gt(nr) = Vm/kMg24Gt(j)+1/(2*kMg24Gt(j))*(ve(j)+lambdaMg24*Vm)*deltat+(1+0.5*lambdaMg24*deltat)*NGt*0.794*4/3*pi*(rGt(nr)^3-rGt(nr-1)^3);
	bottomBMg24Gt(nr) = Vm/kMg24Gt(j-1)-1/(2*kMg24Gt(j-1))*(ve(j-1)+lambdaMg24*Vm)*deltat+(1-0.5*lambdaMg24*deltat)*NGt*0.794*4/3*pi*(rGt(nr)^3-rGt(nr-1)^3);
	middleAMg24Bt(nr) = Vm/kMg24Bt(j)+1/(2*kMg24Bt(j))*(ve(j)+lambdaMg24*Vm)*deltat+(1+0.5*lambdaMg24*deltat)*0.794*4/3*pi*(rBt(nr)^3-rBt(nr-1)^3);
	middleBMg24Bt(nr) = Vm/kMg24Bt(j-1)-1/(2*kMg24Bt(j-1))*(ve(j-1)+lambdaMg24*Vm)*deltat+(1-0.5*lambdaMg24*deltat)*0.794*4/3*pi*(rBt(nr)^3-rBt(nr-1)^3);
	middleAMg24Gt(nr) = (1+0.5*lambdaMg24*deltat)*NGt*0.794*4/3*pi*(rGt(nr)^3-rGt(nr-1)^3);
	middleBMg24Gt(nr) = (1-0.5*lambdaMg24*deltat)*NGt*0.794*4/3*pi*(rGt(nr)^3-rGt(nr-1)^3);
	for i = 2:nr-1
		bottomAMg24Bt(i) = (1+0.5*lambdaMg24*deltat)*4/3*pi*(0.206*rBt(i+1)^3+0.588*rBt(i)^3-0.794*rBt(i-1)^3);
		bottomBMg24Bt(i) = (1-0.5*lambdaMg24*deltat)*4/3*pi*(0.206*rBt(i+1)^3+0.588*rBt(i)^3-0.794*rBt(i-1)^3);
		bottomAMg24Gt(i) = (1+0.5*lambdaMg24*deltat)*NGt*4/3*pi*(0.206*rGt(i+1)^3+0.588*rGt(i)^3-0.794*rGt(i-1)^3);
		bottomBMg24Gt(i) = (1-0.5*lambdaMg24*deltat)*NGt*4/3*pi*(0.206*rGt(i+1)^3+0.588*rGt(i)^3-0.794*rGt(i-1)^3);
	end
	AMg24 = zeros(2*nr,2*nr);
	AMg24(1:nr,1:nr) = AMg24Bt;
	AMg24(nr+1:2*nr,nr+1:2*nr) = AMg24Gt;
	AMg24(nr,1:nr-1) = bottomAMg24Bt(1:nr-1);
	AMg24(nr,nr) = middleAMg24Bt(nr);
	AMg24(nr,nr+1:2*nr-1) = bottomAMg24Gt(1:nr-1);
	AMg24(nr,2*nr) = middleAMg24Gt(nr);
	AMg24(2*nr,1:nr) = bottomAMg24Bt;
	AMg24(2*nr,nr+1:2*nr) = bottomAMg24Gt;
	BMg24 = zeros(2*nr,2*nr);
	BMg24(1:nr,1:nr) = BMg24Bt;
	BMg24(nr+1:2*nr,nr+1:2*nr) = BMg24Gt;
	BMg24(nr,1:nr-1) = bottomBMg24Bt(1:nr-1);
	BMg24(nr,nr) = middleBMg24Bt(nr);
	BMg24(nr,nr+1:2*nr-1) = bottomBMg24Gt(1:nr-1);
	BMg24(nr,2*nr) = middleBMg24Gt(nr);
	BMg24(2*nr,1:nr) = bottomBMg24Bt;
	BMg24(2*nr,nr+1:2*nr) = bottomBMg24Gt;
	CsMg24(:,j) = inv(AMg24)*BMg24*CsMg24(:,j-1);
	massmMg24(j) = CsMg24(nr,j)*Vm/kMg24Bt(j);
	
	% Correction term for flux of Th230 into melt, subtracting the amount of Th230 decayed from the solid and adding that produced in the solid
	% B matrices: (1-0.5*deltat*(lambdaTh-lambdaU*mTh/mU/CsTh(n,j-1)*(CsU(n,j)+CsU(n,j-1))))*
	% A matrices: (1+0.5*lambdaTh*deltat)*
	
	centerAFe56Bt = (1+(DFe56Bt(j)/(RBt/nr)^2+0.5*lambdaFe56)*deltat).*ones(nr,1);
	centerAFe56Bt(1) = 1+(3*DFe56Bt(j)/(RBt/nr)^2+0.5*lambdaFe56)*deltat;
	AdFe56Bt = diag(centerAFe56Bt);
	centerAFe56Gt = (1+(DFe56Gt(j)/(RGt/nr)^2+0.5*lambdaFe56)*deltat).*ones(nr,1);
	centerAFe56Gt(1) = 1+(3*DFe56Gt(j)/(RGt/nr)^2+0.5*lambdaFe56)*deltat;
	AdFe56Gt = diag(centerAFe56Gt);
	centerBFe56Bt = (1-(DFe56Bt(j-1)/(RBt/nr)^2+0.5*lambdaFe56)*deltat).*ones(nr,1)...
			+(0.5*lambdaMg24*deltat).*(CsMg24(1:nr,j)+CsMg24(1:nr,j-1))./CsFe56(1:nr,j-1);
	centerBFe56Bt(1) = 1-(3*DFe56Bt(j-1)/(RBt/nr)^2+0.5*lambdaFe56...
			-0.5*lambdaMg24*(CsMg24(1,j)+CsMg24(1,j-1))/CsFe56(1,j-1))*deltat;
	BdFe56Bt = diag(centerBFe56Bt);
	centerBFe56Gt = (1-(DFe56Gt(j-1)/(RGt/nr)^2+0.5*lambdaFe56)*deltat).*ones(nr,1)...
			+(0.5*lambdaMg24*deltat).*(CsMg24(nr+1:2*nr,j)+CsMg24(nr+1:2*nr,j-1))./CsFe56(nr+1:2*nr,j-1);
	centerBFe56Gt(1) = 1-(3*DFe56Gt(j-1)/(RGt/nr)^2+0.5*lambdaFe56...
			-0.5*lambdaMg24*(CsMg24(nr+1,j)+CsMg24(nr+1,j-1))/CsFe56(nr+1,j-1))*deltat;
	BdFe56Gt = diag(centerBFe56Gt);
	upperAFe56Bt(1) = -3*DFe56Bt(j)*deltat/(RBt/nr)^2;
	lowerAFe56Bt(1) = 0;
	upperBFe56Bt(1) = 3*DFe56Bt(j-1)*deltat/(RBt/nr)^2;
	lowerBFe56Bt(1) = 0;
	upperAFe56Gt(1) = -3*DFe56Gt(j)*deltat/(RGt/nr)^2;
	lowerAFe56Gt(1) = 0;
	upperBFe56Gt(1) = 3*DFe56Gt(j-1)*deltat/(RGt/nr)^2;
	lowerBFe56Gt(1) = 0;
	for i = 2:nr-1;
		upperAFe56Bt(i) = -(i)*DFe56Bt(j)*deltat/(2*(i-1)*(RBt/nr)^2);
		lowerAFe56Bt(i) = -(i-1)*DFe56Bt(j)*deltat/(2*i*(RBt/nr)^2);
		upperBFe56Bt(i) = (i)*DFe56Bt(j-1)*deltat/(2*(i-1)*(RBt/nr)^2);
		lowerBFe56Bt(i) = (i-1)*DFe56Bt(j-1)*deltat/(2*i*(RBt/nr)^2);
		upperAFe56Gt(i) = -(i)*DFe56Gt(j)*deltat/(2*(i-1)*(RGt/nr)^2);
		lowerAFe56Gt(i) = -(i-1)*DFe56Gt(j)*deltat/(2*i*(RGt/nr)^2);
		upperBFe56Gt(i) = (i)*DFe56Gt(j-1)*deltat/(2*(i-1)*(RGt/nr)^2);
		lowerBFe56Gt(i) = (i-1)*DFe56Gt(j-1)*deltat/(2*i*(RGt/nr)^2);
	end
	AudFe56Bt = diag(upperAFe56Bt,1);
	AldFe56Bt = diag(lowerAFe56Bt,-1);
	BudFe56Bt = diag(upperBFe56Bt,1);
	BldFe56Bt = diag(lowerBFe56Bt,-1); 
	AudFe56Gt = diag(upperAFe56Gt,1);
	AldFe56Gt = diag(lowerAFe56Gt,-1);
	BudFe56Gt = diag(upperBFe56Gt,1);
	BldFe56Gt = diag(lowerBFe56Gt,-1);
	AFe56Bt = AdFe56Bt + AudFe56Bt + AldFe56Bt;
	BFe56Bt = BdFe56Bt + BudFe56Bt + BldFe56Bt;
	AFe56Gt = AdFe56Gt + AudFe56Gt + AldFe56Gt;
	BFe56Gt = BdFe56Gt + BudFe56Gt + BldFe56Gt;
	bottomAFe56Bt(1) = (1+0.5*lambdaFe56*deltat)*0.206*4/3*pi*(rBt(2)^3-rBt(1)^3);
	bottomBFe56Bt(1) = (1-0.5*deltat*(lambdaFe56-lambdaMg24/CsFe56(1,j-1)*(CsMg24(1,j)+CsMg24(1,j-1))))*0.206*4/3*pi*(rBt(2)^3-rBt(1)^3);
	bottomAFe56Bt(nr) = (1+0.5*lambdaFe56*deltat)*0.794*4/3*pi*(rBt(nr)^3-rBt(nr-1)^3);
	bottomBFe56Bt(nr) = (1-0.5*deltat*(lambdaFe56-lambdaMg24/CsFe56(nr,j-1)*(CsMg24(nr,j)+CsMg24(nr,j-1))))*0.794*4/3*pi*(rBt(nr)^3-rBt(nr-1)^3);
	bottomAFe56Gt(1) = (1+0.5*lambdaFe56*deltat)*NGt*0.206*4/3*pi*(rGt(2)^3-rGt(1)^3);
	bottomBFe56Gt(1) = (1-0.5*deltat*(lambdaFe56-lambdaMg24/CsFe56(nr+1,j-1)*(CsMg24(nr+1,j)+CsMg24(nr+1,j-1))))*NGt*0.206*4/3*pi*(rGt(2)^3-rGt(1)^3);
	bottomAFe56Gt(nr) = Vm/kFe56Gt(j)+1/(2*kFe56Gt(j))*(lambdaFe56*Vm)*deltat+(1+0.5*lambdaFe56*deltat)*NGt*0.794*4/3*pi*(rGt(nr)^3-rGt(nr-1)^3);
	bottomBFe56Gt(nr) = Vm/kFe56Gt(j-1)-1/(2*kFe56Gt(j-1))*(lambdaFe56*Vm-lambdaMg24*kFe56Gt(j-1)/kMg24Gt(j-1)/CsFe56(2*nr,j-1)*(Vm*CsMg24(2*nr,j)+Vm*CsMg24(2*nr,j-1)))*deltat+...
		(1-0.5*deltat*(lambdaFe56-lambdaMg24/CsFe56(2*nr,j-1)*(CsMg24(2*nr,j)+CsMg24(2*nr,j-1))))*NGt*0.794*4/3*pi*(rGt(nr)^3-rGt(nr-1)^3);
	middleAFe56Bt(nr) = Vm/kFe56Bt(j)+1/(2*kFe56Bt(j))*(lambdaFe56*Vm)*deltat+(1+0.5*lambdaFe56*deltat)*0.794*4/3*pi*(rBt(nr)^3-rBt(nr-1)^3);
	middleBFe56Bt(nr) = Vm/kFe56Bt(j-1)-1/(2*kFe56Bt(j-1))*(lambdaFe56*Vm-lambdaMg24*kFe56Bt(j-1)/kMg24Bt(j-1)/CsFe56(nr,j-1)*(Vm*CsMg24(nr,j)+Vm*CsMg24(nr,j-1)))*deltat+...
		(1-0.5*deltat*(lambdaFe56-lambdaMg24/CsFe56(nr,j-1)*(CsMg24(nr,j)+CsMg24(nr,j-1))))*0.794*4/3*pi*(rBt(nr)^3-rBt(nr-1)^3);
	middleAFe56Gt(nr) = (1+0.5*lambdaFe56*deltat)*NGt*0.794*4/3*pi*(rGt(nr)^3-rGt(nr-1)^3);
	middleBFe56Gt(nr) = (1-0.5*deltat*(lambdaFe56-lambdaMg24/CsFe56(2*nr,j-1)*(CsMg24(2*nr,j)+CsMg24(2*nr,j-1))))*NGt*0.794*4/3*pi*(rGt(nr)^3-rGt(nr-1)^3);
	for i = 2:nr-1
		bottomAFe56Bt(i) = (1+0.5*lambdaFe56*deltat)*4/3*pi*(0.206*rBt(i+1)^3+0.588*rBt(i)^3-0.794*rBt(i-1)^3);
		bottomBFe56Bt(i) = (1-0.5*deltat*(lambdaFe56-lambdaMg24/CsFe56(i,j-1)*(CsMg24(i,j)+CsMg24(i,j-1))))*4/3*pi*(0.206*rBt(i+1)^3+0.588*rBt(i)^3-0.794*rBt(i-1)^3);
		bottomAFe56Gt(i) = (1+0.5*lambdaFe56*deltat)*NGt*4/3*pi*(0.206*rGt(i+1)^3+0.588*rGt(i)^3-0.794*rGt(i-1)^3);
		bottomBFe56Gt(i) = (1-0.5*deltat*(lambdaFe56-lambdaMg24/CsFe56(nr+i,j-1)*(CsMg24(nr+i,j)+CsMg24(nr+i,j-1))))*NGt*4/3*pi*(0.206*rGt(i+1)^3+0.588*rGt(i)^3-0.794*rGt(i-1)^3);
	end
	AFe56 = zeros(2*nr,2*nr);
	AFe56(1:nr,1:nr) = AFe56Bt;
	AFe56(nr+1:2*nr,nr+1:2*nr) = AFe56Gt;
	AFe56(nr,1:nr-1) = bottomAFe56Bt(1:nr-1);
	AFe56(nr,nr) = middleAFe56Bt(nr);
	AFe56(nr,nr+1:2*nr-1) = bottomAFe56Gt(1:nr-1);
	AFe56(nr,2*nr) = middleAFe56Gt(nr);
	AFe56(2*nr,1:nr) = bottomAFe56Bt;
	AFe56(2*nr,nr+1:2*nr) = bottomAFe56Gt;
	BFe56 = zeros(2*nr,2*nr);
	BFe56(1:nr,1:nr) = BFe56Bt;
	BFe56(nr+1:2*nr,nr+1:2*nr) = BFe56Gt;
	BFe56(nr,1:nr-1) = bottomBFe56Bt(1:nr-1);
	BFe56(nr,nr) = middleBFe56Bt(nr);
	BFe56(nr,nr+1:2*nr-1) = bottomBFe56Gt(1:nr-1);
	BFe56(nr,2*nr) = middleBFe56Gt(nr);
	BFe56(2*nr,1:nr) = bottomBFe56Bt;
	BFe56(2*nr,nr+1:2*nr) = bottomBFe56Gt;
	CsFe56(:,j) = inv(AFe56)*BFe56*CsFe56(:,j-1);
	massmFe56(j) = CsFe56(nr,j)*Vm/kFe56Bt(j);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % To include entire crystal, use 1:(nr-1); to exclude rim, use 1:(nr-2)
    %
	for i = 1:(nr-2)
		massBtMg24(j) = massBtMg24(j)+(0.794*CsMg24(i+1,j)+(1-0.794)*CsMg24(i,j))*(4/3*pi)*((rBt(i+1)^3)-(rBt(i)^3));
		massBtFe56(j) = massBtFe56(j)+(0.794*CsFe56(i+1,j)+(1-0.794)*CsFe56(i,j))*(4/3*pi)*((rBt(i+1)^3)-(rBt(i)^3));
		massGtMg24(j) = massGtMg24(j)+NGt*(0.794*CsMg24(nr+i+1,j)+(1-0.794)*CsMg24(nr+i,j))*(4/3*pi)*((rGt(i+1)^3)-(rGt(i)^3));
		massGtFe56(j) = massGtFe56(j)+NGt*(0.794*CsFe56(nr+i+1,j)+(1-0.794)*CsFe56(nr+i,j))*(4/3*pi)*((rGt(i+1)^3)-(rGt(i)^3));
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%						Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=10;					% # timesteps to plot
avCsMg24 = zeros(1,n);
avCsFe56 = zeros(1,n);

for j = 1:n
	avCBtMg24(j) = massBtMg24(j*nt/n)/(4/3*pi*((RBt*(1-1/nr))^3));
	avCBtFe56(j) = massBtFe56(j*nt/n)/(4/3*pi*((RBt*(1-1/nr))^3));
	avCGtMg24(j) = massGtMg24(j*nt/n)/NGt/(4/3*pi*((RGt*(1-1/nr))^3));
	avCGtFe56(j) = massGtFe56(j*nt/n)/NGt/(4/3*pi*((RGt*(1-1/nr))^3));
	avCsMg24(j) = (massBtMg24(j*nt/n)+massGtMg24(j*nt/n))/V;
	avCsFe56(j) = (massBtFe56(j*nt/n)+massGtFe56(j*nt/n))/V;
end

Bt = [avCBtMg24;avCBtFe56]';
Gt = [avCGtMg24;avCGtFe56]';
Sol = [avCsMg24;avCsFe56]';

toc

massMg24 = massBtMg24 + massGtMg24 + massmMg24;
massFe56 = massBtFe56 + massGtFe56 + massmFe56;

alpha_app = (massGtFe56./massGtMg24)./(massBtFe56./massBtMg24); % apparent fractionation factor (vector)
T_isotope = ((2105-1458)./log(alpha_app)).^0.5 - 273*ones(size(alpha_app)); % apparent temperature based on bulk distribution of 44Ca/40Ca between Gt and Bt
alpha_app_final = (massGtFe56(nt)/massGtMg24(nt))/(massBtFe56(nt)/massBtMg24(nt));
%alpha_app_final - 1
T_isotope_final = ((2105-1458)/log(alpha_app_final))^0.5 - 273

K_1 = (massGtFe56./VGt)./(massBtFe56./VBt);
K_2 = (massGtMg24./VGt)./(massBtMg24./VBt);
K_app = K_2./K_1;
for j = 1:nt
   T_app(j) =((-52107.536)./(3*8.314*log(K_app(j))-19.506))-273;
end
% K_app_final = (massGtMg24(nt)/VGt)/(massBtMg24(nt)/VBt);
% T_app_final = 1/((-b-sqrt((b^2)-(4*a*(c-log(K_app_final)))))/(2*a)) - 273

%Delta44Gt_Bt = 1000*log((massGtFe56./massGtMg24)./(massBtFe56./massBtMg24));
%Delta44Gt_Bt_final = Delta44Gt_Bt(nt)
%DCa_app_final = (massGtMg24(nt)/VGt)/(massBtMg24(nt)/VBt)


%figure(1)
%plot(t./3.1536e13,massBtMg24,'g-',t./3.1536e13,massmMg24,'b-')
%xlabel('Time (Myr)')
%ylabel('28Si in Gt (green), Btal (blue)')

figure(2)
hold off
plot(rBt-RBt*ones(size(rBt)),((CsFe56(1:nr,1)./CsMg24(1:nr,1))/(CoFe56/CoMg24)-1)*1000,'k-')
hold on
plot(RGt*ones(size(rGt))-rGt,((CsFe56(nr+1:2*nr,1)./CsMg24(nr+1:2*nr,1))/(CoFe56/CoMg24)-1)*1000,'k-')
for j = 1:n
    plot(rBt-RBt*ones(size(rBt)),((CsFe56(1:nr,j*nt/n)./CsMg24(1:nr,j*nt/n))/(CoFe56/CoMg24)-1)*1000,'k-')
    plot(RGt*ones(size(rGt))-rGt,((CsFe56(nr+1:2*nr,j*nt/n)./CsMg24(nr+1:2*nr,j*nt/n))/(CoFe56/CoMg24)-1)*1000,'k-')
end
xlabel('Distance from Rim')
ylabel('d^{44}Ca')
title('Bt                                   Gt')

figure(3)
hold off
semilogy(rBt-RBt*ones(size(rBt)),CsMg24(1:nr,1),'k-',RGt*ones(size(rGt))-rGt,CsMg24(nr+1:2*nr,1),'k-')
hold on
for j = 1:n
	semilogy(rBt-RBt*ones(size(rBt)),CsMg24(1:nr,j*nt/n),'k-',RGt*ones(size(rGt))-rGt,CsMg24(nr+1:2*nr,j*nt/n),'k-')
end

xlabel('Radial distance')
ylabel('Mg concentration')
%axis([-2e-5 0 0 5])
%hold off
title('Bt                                   Gt')


figure(4)
hold off
semilogy(rBt-RBt*ones(size(rBt)),CsFe56(1:nr,1),'k-',RGt*ones(size(rGt))-rGt,CsFe56(nr+1:2*nr,1),'k-')
hold on
for j = 1:n
	semilogy(rBt-RBt*ones(size(rBt)),CsFe56(1:nr,j*nt/n),'k-',RGt*ones(size(rGt))-rGt,CsFe56(nr+1:2*nr,j*nt/n),'k-')
end

xlabel('Radial distance')
ylabel('Fe concentration')
%axis([-2e-5 0 0 5])
%hold off
title('Bt                                   Gt')


% figure(4)
% hold off
% plot(t./3.1536e13,1000*log((massGtFe56./massGtMg24)./(massBtFe56./massBtMg24)),'k-',t./3.1536e13,1000*log(alphaGtgb./alphaBtgb),'k--')
% xlabel('Time (Myr)')
% ylabel('Delta ^{44}Ca/^{40}Ca Gt/Bt')
% 
figure(5)
hold off
plot(T-273*ones(size(T)),T_app,'b-',[T(nt)-273 T(1)-273],[T(nt)-273 T(1)-273],'k--')
xlabel('Real Temperature (^{o}C)')
ylabel('Apparent Temperature (^{o}C)')
% 
% figure(6)
% hold off
% plot(T_app,Delta44Gt_Bt,'k-')
% xlabel('Ca Gt/Bt temperature')
% ylabel('Delta ^{44}Ca/^{40}Ca Gt/Bt')
% 
% figure(7)
% hold on
% plot(VBt/VGt*massGtMg24./massBtMg24,Delta44Gt_Bt,'k-')
% %hold on
% %plot(Perid_DCa,Perid_Delta44,'ko')
% xlabel('Ca Gt/Bt apparent partition coefficient')
% ylabel('Delta ^{44}Ca/^{40}Ca Gt/Bt')
% 
% error=std(((CsFe56(1:nr,nt)./CsMg24(1:nr,nt))/(CoFe56/CoMg24)-1)*1000)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Comparison statements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Mg Bt----');
disp(XmgBt);
disp(CoMg24Bt);
disp('Fe Bt----');
disp(XfeBt);
disp(CoFe56Bt);
disp('Mg Grt----');
disp(XmgGrt);
disp(CoMg24Gt);
disp('fe Grt----');
disp(XfeGrt);
disp(CoFe56Gt);
