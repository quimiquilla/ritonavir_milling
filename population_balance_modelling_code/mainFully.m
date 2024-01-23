% run this to start the simulation
clear
close all

coolors = [[205, 11, 22]; [249, 65, 74]; [252, 168, 173]; [0 0 0]; [100 100 100]]/255;
lineW = [1, 1.5 ,2, 2.5,3];
lineS = {'-', ':'};


%% parameters
% (1) 
% distribution, values relative to x0
ID = 4; 
mewVec = [200]  ; % at 20 then xinit = 200 microns
sigmaVec = [1];
hurVec = [5e-2]; %%  Nature: [1.25e-2 2.5e-2 5e-2 1e-1 2e-1]% hold up ratio, I think is the amount of seeds % 0.15

% grid specs
dy = 0.1; % Nature simulations run at dy = 0.4
Lmax = 205; % scaled

% mill constant, arbitrary
pmillVec = [2.5e-4 5e-4 1e-3 2e-3 4e-3]; % Nature: [2.5e-4 5e-4 1e-3 2e-3 4e-3];

% conditions to run mill growth or diss
millCond = 1;
growthCond =1;
dissCond = 1;

% nothing happened: -15
Si = 1; % initial supersaturation

% growth rate mechanism
beta = [0]; % growth rate equation parameter, 0 for reaction-limited, -1 for diffusion-limited
aVec = [1e-2];%Nature: [1e-2];   % alpha asterisk in Martin's paper

% time to run for or until to reach equilibrium
Duration = 1e4; % s / s, amount of scaled time, simulations = 1e6
condEquilibrium = 0;    % 1 to stop when tolSuperaturation is met, simulations = 1
tolSuperaturation = 1e-9;   % the dS/dt required to reach equilibrium

% fac = kv * rho / cb; fac and cb will change if cond_keep_m3_constant = 1
kv = 0.5; rho = 2000; cb = 10;
facVec = [kv*rho./cb]; %99.2; % when cond_keep_m3_constant = 1 this is calculated later

cond_keep_m3_constant_vec = [1]; % 1 to change fac (solubility), 0 to change m3
nParticles = 1;

% breakage variables
HzVec = 1; % 1/s, mill frequency

% plotting
condHold = 0;
cMech = 1; % 1 for reaction-limited growth, just for plotting


% Growth constant for checks only, otherwise = 1
kg = 1;

% 1 to remove negatives
condNegs = 0;

% plotting
condPlot = 1;
approxSteps = 1000000;  % also max steps
plotLive = 0; % plot the live distriubiton
plotVol = 0; % 1 to plot vol-weighted
pSteps = 10; % interval of steps until next plot, including saving



maxCFL = 0.9; % multiplies the dt from the CFL condition
maxDS = 1e-5; % max supersaturation change in each step
maxdt = 1e0; 

Conditions = combvec(mewVec,sigmaVec,hurVec,HzVec ,pmillVec ,aVec ,cond_keep_m3_constant_vec ,facVec);
Folder_save = 'Results';
save(fullfile(Folder_save,'Conditions.mat'),'Conditions')
for cIter = 1:size(Conditions,2)
    disp(cIter)
    tic
    
    
    if ~condHold
        close all
    end
    clearvars -except condHold sigmaVec mewVec rhoSVec cIter aVec plotLive...
        approxSteps Duration ID Lmax Nbins coolors A1 dMinus1 dPlus1 daughter selection...
        millCond growthCond pmillVec condNegs pSteps hurVec plotVol lineW fac beta...
        maxCFL dissCond Si maxDS y ym dy Geq Nbins a nskip beta betaVec cMech...
        lineS x0 t0 HzVec dyVec kg aVec a maxdt facVec cond_keep_m3_constant_vec...
        Conditions cLine condPlot kv rho cb nParticles...
        condEquilibrium tolSuperaturation Folder_save
    
    
    % choose parameter
    mew = Conditions(1,cIter);
    sigma = Conditions(2,cIter);
    hur = Conditions(3,cIter);
    Hz = Conditions(4,cIter);
    pmill = Conditions(5,cIter);
    a = Conditions(6,cIter);
    cond_keep_m3 = Conditions(7,cIter);
    
    if ~cond_keep_m3
        fac = Conditions(8,cIter);
    end
    %% grid (rescaled)
    y = 0:dy:Lmax;
    ym = y(1:end-1)+dy/2;
    Nbins = length(ym);
    
    %% non-uniform smooth grid
    %         h = 5e-3;
    %         i = 1:600;
    %         gamma_min = 0;
    %         dy = exp(gamma_min + (i-1)*h) * (exp(h) - 1);
    %         y = [0 cumsum(dy)];
    %         ym = y(2:end) - dy/2;
    %         Nbins = length(dy);
    %% growth/dissolution parameters
    %     %physical
    %     rho = 1760; % density, kg/m3 % 1760
    %     kv = pi/6; % shape factor
    %     % solubility
    %     cb = 1;  % at constant temperature, kg/kg % 0.1252
    
    % rescaling
    x0 = 1e-5; % m
    
    %% INITIAL CONDITIONS
    % initialize distribution and growth
    Q = zeros(1,Nbins);
    G = zeros(1,Nbins);
    
    % define initial distribution, scaled
    Q(1,:) = 1./(sigma * sqrt(2*pi)) * exp(-0.5 .* ((ym - mew)/sigma).^2);
    
    % volume of seeds, from hold up ratio equation
    %     m3i = hur * cb / kv / rho;
    
    % code to change the hur by changing m3
    if ~cond_keep_m3
        m3i = hur / fac;
        m3 = calcMoms(3,Q,ym,dy) * x0^3;  % third moment of non scaled distribution
        scF = m3i/m3;
        Q = Q *  scF;    % convert distribution to have m3i
    else
        m0i = calcMoms(0,Q,ym,dy);
        scF = m0i/nParticles;
        Q = Q *  scF;    % increase the number of particles
        m3i = calcMoms(3,Q,ym,dy) * x0^3;
        fac = hur/m3i;
        cb = rho*kv/fac;
        Conditions(8,cIter) = fac;
    end
    
    %% DEFINE BREAKAGE ELEMENTS
    
    millingElements
    
    % to check the balance
    m3i = calcMoms(3,Q,ym,dy);
    %%%% test moments
    %     Q = [0 2 3 5 1 6 0];
    
    
    % plotting
    nR = 3;  % number of rows of subplots
    nC = 4;  % number of cols of subplots
    
    figure(1)
    
    subplot(nR,nC,1)
    if plotVol
        plot(ym, Q .* ym.^3 ,'Color',coolors(cIter,:), 'LineWidth', lineW(cIter), 'LineStyle', lineS{cMech}); title('Initial Distribution'); xlabel('y'); ylabel('n'); hold on
    else
        plot(ym, Q ,'Color',coolors(cIter,:), 'LineWidth', lineW(cIter), 'LineStyle', lineS{cMech}); title('Initial Distribution'); xlabel('y'); ylabel('n'); hold on
    end
    
    
    m0i = calcMoms(0,Q,ym,dy); % calculate initial number of particles, #/kg
    
    
    %% initial values
    
    
    tprof = zeros(approxSteps,Nbins);
    Gprof = zeros(approxSteps,Nbins);
    Qprof = zeros(approxSteps,Nbins);
    Sprof = zeros(approxSteps,1);
    stepprof = zeros(approxSteps,1);
    dQprofMill = zeros(approxSteps,1);
    dQprofGrowth = zeros(approxSteps,1);
    dQprofDiss = zeros(approxSteps,1);
    
    % initial values
    tnew = 0;  % time
    cSteps = 0; % counter of steps
    dtbad = 0; % the number of bad iteration for a step, resets to 0 if the step is fine
    % Initial growth rate
    G =  kg * ym.^(beta) .* (Si - 1 - a./ym) ; %Geq(Si);
    
    m3new = calcMoms(3,Q,ym,dy);
    
    Qprof(1,:) = Q;
    tprof(1) = tnew;
    Sprof(1) = Si;
    Gprof(1,:) = G;
    dt = 1;
    Snew = Si;
    Qnew = Q;
    
    
    
    nbadsteps = 0;
    cSaved = 0; % number of points saved
    %% loop
    while Duration > tnew
        cSteps = cSteps+1;
        if condNegs
            Qnew(Qnew<0) = 0;
        end
        
        m3old = m3new;
        % supersaturation
        Sold = Snew;
        % distribution, Q is Qold
        Q = Qnew;
        % time
        told = tnew;
        
        %% CFL condition, G*dt = dy
        dt = min(min(dt*1.5,maxdt), maxCFL*min(dy./abs(G))) / exp(dtbad); % exp(dtbad) = for each bad iteration the dt is decreased exponentially
        if told + dt > Duration
            dt = Duration-told;
        end
        % first step must be 1e-3 for plotting
        if cSteps == 1 || tnew < 10
            dt = min(dt, 1e-4);
        elseif tnew < 100
            dt = min(dt, 1e-2);
        elseif tnew < 1000
            dt = min(dt, 1e-1);
        elseif tnew < 10000
            dt = min(dt, 1e-0);
        end
        tnew = told + dt;
        
        %% GROWTH/DISSOLUTION
        dQdt = 0;
        if growthCond
            % calculate fluxes for growing bins
            rplus = (Q(3:end) - Q(2:end-1) + 1e-9)./ (Q(2:end-1) - Q(1:end-2) + 1e-9);
            phiplus = max(0, min(2*rplus, min(1/3+2/3*rplus, 2)));
            F = zeros(1,Nbins+1);
            F(1) = 0; % no nucleation
            F(2) = 0.5 * G(1) * (Q(1)+Q(2));
            F(end) =  G(end) * (Q(end) + 0.5*(Q(end) - Q(end-1))); % leaves the grid
            F(3:end-1) = G(2:end-1) .* (Q(2:end-1) + 0.5*phiplus.*(Q(2:end-1) - Q(1:end-2)));
            % keep the fluxes from positive rates
            F([false G<=0]) = 0; % put 1 first for the nucleation flux (it's zero anw)
            dQdt = (F(2:end) - F(1:end-1))./dy; % change in number density
        end
        
        if dissCond
            % calculate fluxes for dissolving bins
            rminus = (Q(1:end-2) - Q(2:end-1) + 1e-9) ./ (Q(2:end-1) - Q(3:end) + 1e-9);
            phiminus = max(0, min(2*rminus, min(1/3+2/3*rminus, 2)));
            F = zeros(1,Nbins+1);
            F(1) = G(1) * (Q(1) + 0.5*(Q(1) - Q(2)));
            F(end-1) = 0.5 * G(end) * (Q(end) + Q(end-1));
            F(end) = 0; % nothing comes from the end
            F(2:end-2) = G(2:end-1) .* (Q(2:end-1) + 0.5*phiminus .* (Q(2:end-1) - Q(3:end)));
            % keep the fluxes from negative rates
            F([G>0 false]) = 0;
            dQdt = dQdt + (F(2:end) - F(1:end-1))./dy; % change in number density
        end
        %% MASS BALANCE
        if growthCond | dissCond
            % update PSD and Conc
            Qnew = Q - dQdt*dt;
            % change in volume - rescaled
            
            m3new = sum(ym.*ym.*ym .* Qnew .* dy);
            dm3 = m3new - m3old;
            
            dS = - fac * dm3 * x0^3; % rescale
            Snew = Sold + dS;
        end
        % save dQdt
        Qsave1 = Qnew;
        %% BREAKAGE
        if millCond
            % convert number to mass, per length
            g_L = Qnew .* ym.^3;
            % convert per length to per volume
            g_x = g_L.*dy./dx;
            % mass density per length, transpose
            g_x = g_x';
            % breakage happening
            g_x = g_x + A1 * g_x * dt;
            %                 figure(4); plot(ym,A1 * g * dt)
            
            g_x = g_x' ;
            %disp(calcMoms(-1,g_x,xm,dx));
            % convert to mass per length
            g_L = g_x.*dx./dy;
            %disp(calcMoms(-3,g_L,ym,dy));
            % convert back to number
            Qnew = g_L ./ ym.^3;
            %disp(calcMoms(0,Qnew,ym,dy));
        end
        Qsave2 = Qnew;
        %% TIME STEP CHECK AND LOG DATA
        % check whether the step was huge or not, by checking the change is supersaturation
        % check if supersaturation went under or above one a lot
        if abs(Snew-Sold) < maxDS  % not a lot = allow iteration
            dtbad = 0;
            %%%%
            % save profiles
            G = kg * ym.^(beta) .* (Snew - 1 - a./ym) ;%Geq(Snew);
            % save first 10 steps then every thousand and last step
            if mod(cSteps,pSteps) == 0 | cSteps < 1000 | Duration == tnew
                cSaved = cSaved + 1;
                Gprof(cSaved+1,:) = G;
                Sprof(cSaved+1) = Snew;
                Qprof(cSaved+1,:) = Qnew;
                tprof(cSaved+1) = tnew;
                stepprof(cSaved) = dt;
                dQsave = dQdt .* dy ;
                dQprofMill(cSaved) = -sum((Qsave1 - Qsave2).*dy / dt);
                dQprofGrowth(cSaved) = sum(dQsave(dQsave>0)) ;
                dQprofDiss(cSaved) = sum(dQsave(dQsave<0)) ;
            end
            % break if everything dissolved
            if m3new/m3i < 1e-5
                disp('99.999% mass dissolved')
                break
            end
            % break if equilibrium is reached
            if condEquilibrium
                if (Snew-Sold)/dt < tolSuperaturation
                    disp(['Equilibrium reached at ' num2str(round(tnew,1))])
                    break
                end
            end
            
            
            % plot evolution
            if mod(cSteps,pSteps)==0  || Duration == tnew
                if plotLive
                    % plot milling effect
                    figure(4); plot(xm,Qnew  .*dy./dx)
                    xlim([0 1e3]); title('Number density per volume, #/L^3'); %set(gca, 'XScale', 'log');
                    figure(2); title('Number density per volume, #/L')
                    if plotVol
                        plot(ym,Qnew .*ym.^3)
                    else
                        plot(ym, Qnew )
                    end
                    %                 set(gca, 'YScale', 'log');
                    hold on;
                    [V,I] = min(abs(G));
                    if isempty(I)
                        I = length(ym);
                    end
                    xline(ym(I),'r-')
                    hold off
                end
            end
            % condition that dt was succesful
            dtbad = max(0,dtbad - 1);
        else  % a lot of change = redo iteration
            % condition to recalculate dt
            dtbad = dtbad + 1;
            % total bad steps
            nbadsteps = nbadsteps + 1;
            % assign values as in previous step
            Qnew = Q;       % Q is Qold
            Snew = Sold;
            tnew = told;
            m3new = m3old;
            cSteps = cSteps - 1;
        end
        %%%%
        % display values
        if mod(cSteps,pSteps*50) == 0 
            disp(['time: ' num2str(tnew) ', dt: ' num2str(dt) ', Supersaturation: ' num2str(Sold) ', dtbad: ' num2str(dtbad)])
        end
    end
    %% remove zeros from pre-defined variables
    
    tprof(cSaved+2:end) = [];
    Sprof(cSaved+2:end) = [];
    Gprof(cSaved+2:end,:) = [];
    Qprof(cSaved+2:end,:) = [];
    stepprof(cSaved+1:end,:) = [];
    dQprofMill(cSaved+1:end) = [];
    dQprofGrowth(cSaved+1:end) = [];
    dQprofDiss(cSaved+1:end) = [];
    
    %% calculate distribution properties, central moments and moments
    % Log data
    m0 = calcMoms(0,Qprof,ym,dy);
    m1 = calcMoms(1,Qprof,ym,dy);
    m2 = calcMoms(2,Qprof,ym,dy);
    m3 = calcMoms(3,Qprof,ym,dy);
    m4 = calcMoms(4,Qprof,ym,dy);
    cm2 = calcCentralMoms(2,Qprof,ym,dy); % ym and dy are vector inputs
    cm3 = calcCentralMoms(3,Qprof,ym,dy);
    cm4 = calcCentralMoms(4,Qprof,ym,dy);
    
    Mean = m1./m0;
    STD = sqrt(cm2./m0); % std
    cov = STD.*m0./m1; % coefficient of variation
    Skew = cm3./(m0.*STD.^3); % skewness
    Kurt = cm4./ (m0 .* STD.^4) - 3;
    %
    %
    %% plot profiles
    if condPlot
        disp(['Bad steps: ' num2str(nbadsteps / (nbadsteps+cSteps)*100) '%'])
        plot_profiles
        %% number and mass check
        % number check
        m0f = calcMoms(0,Q,ym,dy);
        disp(['Number change: ' num2str(round((m0i-m0f)/m0i*100,3)) '%'])
        % mass check
        Merr = ((Sprof(end) + fac * m3(end) * x0^3) - (Sprof(1) + fac * m3i * x0^3))/(Sprof(1) + fac * m3i * x0^3);
        disp(['Mass error: ' num2str( Merr ) '%'])
        
        %% save figs
        fh = figure(1);
        fh.WindowState = 'maximized';
        savefig(fullfile(Folder_save,[num2str(ID) num2str(char(96+cIter)) '_mew = ' num2str(mew) ', sigma = ' num2str(sigma) ', hur = ' num2str(hur) ', Hz = ' num2str(Hz) ', pmill = ' num2str(pmill) ', a = ' num2str(a) ', cond_keep_m3 = ' num2str(cond_keep_m3) ', fac = ' num2str(fac) '.fig']))
        close all
    end
    %% save results
    save(fullfile(Folder_save,'Conditions.mat'),'Conditions')
    save(fullfile(Folder_save,[num2str(ID) num2str(char(96+cIter)) '_mew = ' num2str(mew) ', sigma = ' num2str(sigma) ', hur = ' num2str(hur) ', Hz = ' num2str(Hz) ', pmill = ' num2str(pmill) ', a = ' num2str(a) ', cond_keep_m3 = ' num2str(cond_keep_m3) ', fac = ' num2str(fac) '.mat']))
    cIter = cIter + 1;
    toc
end


function mom = calcMoms(i,Q,ym,dy)
% calculates moments
ymesh = repmat(ym,size(Q,1),1); % make 2D grid
dymesh = repmat(dy,size(Q,1),1);
mom = sum(ymesh.^i .* Q .* dymesh,2);
end

function cmom = calcCentralMoms(i,Q,ym,dy)
% calculates central moments
ymesh = repmat(ym,size(Q,1),1); % make 2D grid
dymesh = repmat(dy,size(Q,1),1);
meanL = calcMoms(1,Q,ym,dy)./calcMoms(0,Q,ym,dy); % find mean
cmom = sum((ymesh - meanL).^i .* Q .* dymesh,2); % find central moment
end



