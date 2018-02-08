clear; clc; close all;
format longG;
set(0,'DefaultFigureWindowStyle','docked')
%% TRUE MODEL PARAMETERS
rho=[500 1e2 200];   % LAYER RESISTIVITIES
h=[10 5];            % LAYER THICKNESSES
TxA=1200;            % TRANSMITTER AREA (SQUARE OR CIRCULAR)
Ic=12;               % TRANSMITTER CURRENT



%% FREQUENCIES
FREQ=logspace(0,5,6); % FREQUENCIES
 
%% J=sensitivity matrix or sensitivity matrix
J=ANALYTIC_SENSITIVITY_FREQ_DOMAIN_EM(rho,h,TxA,Ic,FREQ);
%% MATLAB FIGURES 2;5;7 AND 10 ARE THE ANALYTICAL RESULTS IN FIGURE 3a, 3c, 3b and 3d RESPECTIVELY IN THE PAPER

Jreal=real(J); % REAL      COMPONENT OF THE SENSITIVITY MATRIX
Jimag=imag(J); % IMAGINARY COMPONENT OF THE SENSITIVITY MATRIX


NL=length(rho);
%NP=4*NL-2; % Number of parameters
for i=1:2*NL-1
    figure;
    semilogx(FREQ,Jreal(:,i),'ko-','linewidth',1.5)
    hold on;    
    xlabel('Frequencies (Hz)')
    ylabel('Sensitivity')
    if i<=NL
        title(['Real component sensitivity for layer ',num2str(i),' conductivity'])
    else
        title(['Real component sensitivity for layer ',num2str(mod(i,NL)),' thickness'])
    end
%     end
    grid on;
    set(gca,'fontsize',14,'fontweight','bold')
    


end


for i=1:2*NL-1
    figure;
    semilogx(FREQ,Jimag(:,i),'ko-','linewidth',1.5)
    hold on;    
    xlabel('Frequencies (Hz)')
    ylabel('Sensitivity')
    if i<=NL
        title(['Imaginary component sensitivity for layer ',num2str(i),' conductivity'])
    else
        title(['Imaginary component sensitivity for layer ',num2str(mod(i,NL)),' thickness'])
    end
%     end
    grid on;
    set(gca,'fontsize',14,'fontweight','bold')
    


end



