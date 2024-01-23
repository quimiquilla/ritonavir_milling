% defining elements for milling

%% Daughter and Selection Functions



% parameters
xref = 1e-6; % no units, reference for numerical stability
kb = 1.^3; % no units, minimum breakage size

t0 = rho * x0^(1-beta) / kv / cb; % 45, I made k = kv, it doesn't matter, but k is growth rate constant
cb;

% breakage rate, x is volume
zeta = @(xk) xk .* (1 + exp(- (xk - kb)./xref)) .^ (-1) ;
daughter = @(L,xk) 2./xk + 0*L; % simple
selection = @(xk) pmill .* Hz^2 .* zeta(xk) *x0 .* t0; % not sure if *x0 but it's constant so units are fine

% convert grid to volume
xm = ym.^3;
dx(1) = xm(1)*2;
for cx = 2:length(xm)
    dx(cx) = (xm(cx)-(sum(dx(1:cx-1))))*2;
end

% calculate matrix
if millCond
    if isfile(fullfile(Folder_save,['A_Hz' num2str(Hz) '_pmill' num2str(pmill) '_t0' num2str(t0) '.mat']))
        A = load(fullfile(Folder_save,['A_Hz' num2str(Hz) '_pmill' num2str(pmill) '_t0' num2str(t0) '.mat']));
        A1 = A.A1;
        maxdtmill = A.maxdtmill;
        
    else
        [A1, maxdtmill] = breakage_matrix(Nbins, xm, dx, selection, daughter);
        save(fullfile(Folder_save,['A_Hz' num2str(Hz) '_pmill' num2str(pmill) '_t0' num2str(t0) '.mat']), 'A1', 'maxdtmill')
    end
    maxdt = min(maxdt,maxdtmill);
else
    maxdt = min(maxdt,1);
end
function [A1, maxdtmill] = breakage_matrix(Nbins, ym, dy, selection, daughter)

%% BREAKAGE MATRIX
% define matrices
A1 = zeros(Nbins,Nbins);
dMinus1 = zeros(Nbins,Nbins);
dPlus1 = zeros(Nbins,Nbins);

% breakage fluxes
warning('on','all')

parfor i = 1:Nbins % this can be parfor
    % warnings come from division by 0, but the result is zero which fine
    warning('on','all')
    yMinus = ym(i)-dy(i)/2; %dy(i)
    yPlus = ym(i)+dy(i)/2;
    for k = 1:Nbins
        ykMinus = ym(k) - dy(k)/2; %dy(k)
        ykPlus = ym(k) + dy(k)/2;
        
        % fluxes for jsut length, maybe wrong
        %             selection_integral = integral(@(y) selection(y)./y.^3, ykMinus, ykPlus,'AbsTol', 1e-16, 'RelTol', 1e-12) % ,'AbsTol', 1e-16, 'RelTol', 1e-12
        %             dMinus1(i,k) = - selection_integral  *   integral(@(L) L.^3 .* daughter(L,ym(k)), 0, yMinus,'AbsTol', 1e-16, 'RelTol', 1e-12); % max because it's nan on the yMinus=0, otherwise is positive
        %             dPlus1(i,k) = - selection_integral  * integral(@(L) L.^3 .* daughter(L,ym(k)), 0, yPlus,'AbsTol', 1e-16, 'RelTol', 1e-12);
        
        selection_integral = integral(@(y) selection(y)./y, ykMinus, ykPlus,'AbsTol', 1e-16, 'RelTol', 1e-12) % ,'AbsTol', 1e-16, 'RelTol', 1e-12
        dMinus1(i,k) = - selection_integral  *   integral(@(L) L .* daughter(L,ym(k)), 0, yMinus,'AbsTol', 1e-16, 'RelTol', 1e-12); % max because it's nan on the yMinus=0, otherwise is positive
        dPlus1(i,k) = - selection_integral  * integral(@(L) L .* daughter(L,ym(k)), 0, yPlus,'AbsTol', 1e-16, 'RelTol', 1e-12);
        
        % the befores, might be wrong or only for S(e) = e and b = 1/x
        %                             dMinus1(i,k) = - integral(selection, ykMinus, ykPlus)  *   integral(@(L) daughter(L,ym(k)), 0, yMinus); % max because it's nan on the yMinus=0, otherwise is positive
        %                             dPlus1(i,k) = - integral(selection, ykMinus, ykPlus)  * integral(@(L) daughter(L,ym(k)), 0, yPlus);
        
    end
end
%         dMinus1(isnan(dMinus1)) = 0;
%         dPlus1(isnan(dPlus1)) = 0;
% breakage matrix
for i = 1:Nbins
    A1(i,i) =  1./dy(i) * dMinus1(i,i); %dy(i)
    A1(i,i+1:Nbins) = 1./ dy(i).*( dMinus1(i,i+1:Nbins) - dPlus1(i,i+1:Nbins) );  %dy(i)
    
end

% maximum timestep
for k = 1:length(ym)
    dtmill(k) = 1/ ( 1/dy(k) * sum(selection(ym(k)) / ym(k) * dy(k)) * sum(ym(1:k) .* daughter(ym(1:k), ym(k)) .* dy(k))); % 1/ ( 1/dy(i) * sum(selection(ym(i)) / ym(i) * dy(i)* sum(ym(1:i-1) .* daughter(ym(1:i-1), ym(i)) .* dy(1:i-1))));
end
maxdtmill = min(dtmill);


end