function [fCO2,pH,CO2,HCO3,CO3,omega_cal,omega_arag] = carbonate_eq(S,T,press,DIC,Alk)


% carbonate_eq   Concentrations of carbonate species
%=========================================================================
% carbonate_eq  Version 1.4 (June 14, 2010)
%
% USAGE:  [fCO2,pH,CO2,HCO3,CO3,omega_cal,omega_arg] = carbonate_eq(S,T,press,DIC,Alk)
%
% DESCRIPTION:
%    Concentrations of carbonate species including fCO2 and pH from DIC and
%    Alkalinity as a function of temperature, salinity and depth.
%    Also calculates saturation state for calcite and aragonite.
%    Considers only carbon, boron, and water ions in alkalinity.
%
% INPUT:  (if inputs are not singular they must have same dimensions)
%   S = salinity    [PSS]
%   T = temperature [degree C]
%   press = pressure   [dbar]
%   DIC = Dissolved Inorganic Carbon   [micromol/kg]
%   Alk = Alkalinity   [microeq/kg]
%
% OUTPUT:
%   fCO2 = fugacity of CO2  [uatm]
%   pH = -log10 of hydrogen ion concentration
%   CO2 = concentration of dissolved CO2 and H2CO3  [micromol/kg]
%   HCO3 = concentration of dissolved HCO3  [micromol/kg]
%   CO3 = concentration of dissolved CO3  [micromol/kg]
%   omega_cal = degree of saturation for calcite (omega < 1 is
%   undersaturated)
%   omega_arag = degree of saturation for aragonite (omega < 1 is
%   undersaturated)
% 
% AUTHOR:  Roberta Hamme (University of Victoria) rhamme@uvic.ca
%
% REFERENCE:
%    This program uses equilibrium constants and equations from Dickson, A.G., 
%    Sabine, C.L. and Christian, J.R. (Eds.) 2007. Guide to best practices for 
%    ocean CO2 measurements. PICES Special Publication 3, 191 pp.
%    Depth dependence from Millero (Geochim. Cosmochim. Acta 43, 1979).
%    Calcite and aragonite solubility from Mucci "The solubility of calcite 
%    and aragonite in seawater at various salinities, temperatures, and one 
%    atmosphere total pressure" (Am. J. Sci., 283, 1983) with depth
%    dependence of calcite from Ingle (Mar. Chem. 3, 1975) and 
%    depth dependence of aragonite from Millero (Geochim. Cosmochim. Acta 43, 1979).
%
% DISCLAIMER:
%    This software is provided "as is" without warranty of any kind.  
%=========================================================================

% CALLER: general purpose
% CALLEE: none

%----------------------
% Check input parameters
%----------------------
if nargin ~=5
   error('carbonate_eq.m: Must pass 5 parameters')
end %if

% Determine variable dimensions
[ms,ns] = size(S);
[mt,nt] = size(T);
[mz,nz] = size(press);
[md,nd] = size(DIC);
[ma,na] = size(Alk);
  
% Check that inputs have the same shape or are singular
if (ms+ns)>2
    mall=ms; nall=ns;
end

if (mt+nt)>2
    if exist('mall','var')
        if (mt ~= mall || nt ~= nall)
            error('carbonate_eq.m: inputs must have same dimensions or be singular')
        end
    else mall=mt; nall=nt;
    end
end

if (mz+nz)>2
    if exist('mall','var')
        if (mz ~= mall || nz ~= nall)
            error('carbonate_eq.m: inputs must have same dimensions or be singular')
        end
    else mall=mz; nall=nz;
    end
end

if (md+nd)>2
    if exist('mall','var')
        if (md ~= mall || nd ~= nall)
            error('carbonate_eq.m: inputs must have same dimensions or be singular')
        end
    else mall=md; nall=nd;
    end
end

if (ma+na)>2 && exist('mall','var') && (ma ~= mall || na ~= nall)
    error('carbonate_eq.m: inputs must have same dimensions or be singular')
end

if ~exist('mall','var')
    mall = ma; nall=na;
end

%----------------------
% Do some basic unit conversion
%----------------------

TK = T + 273.15;        % Convert temperature to Kelvin
Pr = press ./ 10;       % Conversion of dbar to bar
Alk = Alk * .000001;    % Convert DIC to mol/kg
DIC = DIC * .000001;    % Convert Alkalinity to eq/kg
R = 83.14472;             % Gas constant

%----------------------
% Calculate equilibrium constants
%----------------------

% Calculate total borate (BT) from salinity (Uppström (1974))

BT = .000416 .* S ./ 35.0;

% Calculate total calcium from salinity (Riley and Tongudai 1967)

Ca = 0.01028 .* S ./ 35.0;  %(mol/kg)

% Calculate Henry's Law coeff (KH) from temp & sal (Weiss, 1974, Mar. Chem vol2, 203-215)

U1 = -60.2409 + 93.4517 * (100 ./ TK) + 23.3585 * log(TK ./ 100);
U2 = S .* (0.023517 - 0.023656 .* (TK ./ 100) + 0.0047036 * (TK ./ 100) .^ 2);
KH = exp(U1 + U2);

% Calculate borate equil constant (KB) from temp & sal (Dickson, 1990, Deep-Sea Res vol37, 755-766)

KB = exp((-8966.9 - 2890.53 * S.^0.5 - 77.942 * S + 1.728 * S.^1.5 - 0.0996 * S.^2)./TK...
   + 148.0248 + 137.1942 * S.^0.5 + 1.62142 * S - (24.4344 + 25.085 * S.^0.5 +...
   0.2474 * S) .* log(TK) + 0.053105 * S.^0.5 .* TK);

% Calculate carbonate equil constants (K1 & K2) from temp & sal (Lueker et al., 2000, Mar Chem vol70, 105-119)

K1 = 10.^(-(3633.86./TK - 61.2172 + 9.67770 * log(TK)- 0.011555*S + 0.0001152 * S.^2));
K2 = 10.^(-(471.78./TK + 25.9290 - 3.16967 * log(TK) - 0.01781*S + 0.0001122 * S.^2));

% Calculate effect of pressure on K1, K2, and KB (Millero, 1979, GCA vol43, 1651-1661)

dvB = -29.48 + 0.295 * (S-34.8) + 0.1622 * T - .002608 * T.^2;
dv1 = -25.50 - 0.151 * (S-34.8) + 0.1271 * T;
dv2 = -15.82 + 0.321 * (S-34.8) - 0.0219 * T;

dkB = -0.00284 + 3.54e-4 * (S-34.8);
dk1 = -0.00308 - 5.78e-4 * (S-34.8) + 0.0000877 * T;
dk2 = 0.00113 - 3.14e-4 * (S-34.8) - 0.0001475  * T;

KB = (exp(-(dvB./(R*TK)).*Pr + (0.5 * dkB./(R*TK)).*Pr.^2)) .* KB;
K1 = (exp(-(dv1./(R*TK)).*Pr + (0.5 * dk1./(R*TK)).*Pr.^2)) .* K1;
K2 = (exp(-(dv2./(R*TK)).*Pr + (0.5 * dk2./(R*TK)).*Pr.^2)) .* K2;

% Calculate water equil constant (KW) from temp & sal (Millero, 1995, GCA vol59, 661-677)

KW1 = 148.9652-13847.26./TK-23.6521*log(TK);
KW2 = (118.67./TK-5.977+1.0495*log(TK)).*S.^0.5-0.01615*S;
KW = exp(KW1+KW2);

% Calculate effect of pressure on KW from CO2SYS.m function (This is from
% Millero, 1983 and his programs CO2ROY(T).BAS.)

dvW  = -20.02 + 0.1119.*T - 0.001409.*T.^2;
dkW = (-5.13 + 0.0794.*T)./1000; % Millero, 1983
KW = (exp(-(dvW./(R*TK)).*Pr + (0.5 * dkW./(R*TK)).*Pr.^2)) .* KW;

% Calculate Ksp for calcite and argonite at total pressure of 1 atm 

Ksp_calcite =  10 .^(-171.9065 - 0.077993.*TK + 2839.319./TK + 71.595.*log10(TK) + (-0.77712 + 0.0028426.*TK + 178.34./TK) .* S.^0.5 - 0.07711.*S + 0.0041249.*S.^1.5);
Ksp_aragonite = 10 .^(-171.945 - 0.077993.*TK + 2903.293./TK + 71.595.*log10(TK) + (-0.068393 + 0.0017276.*TK + 88.135./TK) .* S.^0.5 - 0.10018.*S + 0.0059415.*S.^1.5);

% Correct Ksp for calcite and argonite for pressure effect

deltaVK_calcite = -48.76 + 0.5304 .* T;
KappaK_calcite  = (-11.76 + 0.3692 .* T)./1000;
lnKCafac  = (-deltaVK_calcite + 0.5.*KappaK_calcite.*Pr).*Pr./ (R .* TK);
Ksp_calcite = Ksp_calcite .* exp(lnKCafac);
       
deltaVK_aragonite = deltaVK_calcite + 2.8;
KappaK_aragonite  = KappaK_calcite;
lnKArfac  = (-deltaVK_aragonite + 0.5.*KappaK_aragonite.*Pr).*Pr./ (R .* TK);
Ksp_aragonite = Ksp_aragonite .* exp(lnKArfac);


%----------------------
% Solve for H ion concentration using a minimization routine
% first getting inputs into correct dimensions
%----------------------

h = zeros(mall,nall);
for i = 1:mall
    for j = 1:nall
        if ma+na > 2
            Alk_loop = Alk(i,j);
        else Alk_loop = Alk;
        end
        if md+nd > 2
            DIC_loop = DIC(i,j);
        else DIC_loop = DIC;
        end
        if (ms+ns > 2)
            BT_loop = BT(i,j);
        else BT_loop = BT;
        end
        if (ms+ns > 2) || (mt+nt > 2) || (mz+nz > 2)
            KW_loop = KW(i,j);
        else KW_loop = KW;
        end
        if (ms+ns > 2) || (mt+nt > 2) || (mz+nz > 2)
            K1_loop = K1(i,j);
            K2_loop = K2(i,j);
            KB_loop = KB(i,j);
        else K1_loop = K1;
            K2_loop = K2;
            KB_loop = KB;
        end
        if sum(isnan([Alk_loop DIC_loop K1_loop])) == 0
            h(i,j) = fminsearch(@(H) H_fitsearch(H,Alk_loop,DIC_loop,K1_loop,K2_loop,KB_loop,KW_loop,BT_loop),1e-10);
        else h(i,j) = NaN;
        end        

    end
end

%----------------------
% calculate HCO3, CO3 and CO2, fCO2 and pH
%----------------------

HCO3 = DIC./(1 + h./K1 + K2./h) * 1e6;
CO3 = DIC./(1 + h./K2 + h.*h./(K1.*K2)) * 1e6;
CO2 = DIC./(1 + K1./h + K1.*K2./(h.*h)) * 1e6;
fCO2 = CO2 ./ KH;
pH=-log10(h);
omega_cal = Ca .* (CO3/1e6) ./ Ksp_calcite;
omega_arag = Ca .* (CO3/1e6) ./ Ksp_aragonite;

%% Subfunction called by primary function to perform fit

function [minimize] = H_fitsearch(H,Alk,DIC,K1,K2,KB,KW,BT)

%----------------------
% Subfunction of carbonate_eq used to calculate least squares residuals for
% current guess of H
%----------------------

minimize = Alk - ((DIC .* K1 .* H ./ (H.^2 + K1 .* H + K1 .* K2)) + ...
    (2 * DIC .* K1 .* K2 ./ (H.^2 + K1 .* H + K1 .* K2)) + (KB .* BT ./ (H + KB)) + (KW ./ H) - H);

minimize = abs(minimize) * 1e6;