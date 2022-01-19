function [T2_combinations] =  TRIR_Inspect_separate_TIs(data, mask, nr_echos, TIs, fitted_T1, fitted_M0, fitted_IE)
% This function performs TRIR fitting for each TI separately
% Copyright 2021, Koen Baas

test_TIs = nr_echos;

% Pre-allocate some variables
T2_combinations = NaN(test_TIs,1);
Errors_rmse = NaN(test_TIs);
Errors_sse = NaN(test_TIs);
Errors_rsquare  =NaN(test_TIs);
Errors_dfe = NaN(test_TIs);
Errors_adjrsquare = NaN(test_TIs);
Errors_CI_T2 = NaN(test_TIs);

Averages = zeros(nr_echos,4);

%% Calculate averages
for i = 1:4
    for echo = 1:nr_echos
            
        slice = data(:,:,(i-1)*nr_echos + echo);
        Averages(echo,i) =  mean(slice(mask));
        
    end
end

cnt = 0;

for use_nr_echos = 1:test_TIs
    
    z1 = Averages(use_nr_echos,1);
    z1(2) = Averages(use_nr_echos,2);
    z1(3) = Averages(use_nr_echos,3);
    z1(4) = Averages(use_nr_echos,4);
    cnt = cnt +1;
    
    y1 = [0 40 80 160]';
    TI = TIs(use_nr_echos);
    
    %Only fit T2
    H = @(T2,x)(fitted_M0*abs(1-(1+exp(-(x./T2))*fitted_IE).*exp(-(TI/fitted_T1))));
       
    try

        [thefit,fit_info] = fit(y1,z1',H,'StartPoint',[60],'Upper',[300],'Lower',[30]);
        
        confidenceint = confint(thefit);
        T2_combinations(use_nr_echos) = thefit.T2;
        
        Errors_rmse(use_nr_echos) = fit_info.rmse;
        Errors_sse(use_nr_echos) = fit_info.sse;
        Errors_rsquare(use_nr_echos) = fit_info.rsquare;
        Errors_dfe(use_nr_echos) = fit_info.dfe;
        Errors_adjrsquare(use_nr_echos) = fit_info.adjrsquare;
        %Errors_CI_T2(use_nr_echos) = (confidenceint(4)-confidenceint(3))/thefit.T2;
        Errors_CI_T2(use_nr_echos) = (confidenceint(2)-confidenceint(1))/thefit.T2;
   
        
        
    catch
        T2_combinations(use_nr_echos) = NaN;
        Errors_rmse(use_nr_echos) = Inf;
        Errors_CI_T2(use_nr_echos) = Inf;
    end
    
 % Plot the first 25 results   
 if use_nr_echos <26
    subplot(5,5,use_nr_echos)
    plot(thefit,y1,z1)
    axis([0 160 0 max(z1)*1.1])
    title(['TI = ' num2str(TI) 'ms' ', T2 = ' sprintf('%1.1f',thefit.T2) 'ms ± ' sprintf('%1.1f',(confidenceint(2)-confidenceint(1)))])
    xlabel('eTE [ms]')
    ylabel('Fitted T2 [ms]')
 end
 
end


%% Smooth curve and get metrics           
figure,
plot(1:use_nr_echos,T2_combinations)
title('T2 per TI')
xlabel('TI used')
ylabel('T2 [ms]')

figure
subplot(4,2,1)
plot(1:use_nr_echos,Errors_CI_T2)
title('Relative CI T2 [%]')
xlabel('TI used')
ylabel('Relative error')
axis([0 test_TIs 0 0.5])

subplot(4,2,2)
plot(1:use_nr_echos,Errors_rmse)
title('rmse')
xlabel('TI used')
ylabel('rmse')
axis([0 test_TIs 0 max(Errors_rmse(:))*1.1])

subplot(4,2,3)
plot(1:use_nr_echos,Errors_sse)
title('sse')
xlabel('TI used')
ylabel('sse')

subplot(4,2,4)
plot(1:use_nr_echos,Errors_rsquare)
title('rsquare')
xlabel('TI used')
ylabel('rsquare')
axis([0 test_TIs 0 1])


subplot(4,2,5)
plot(1:use_nr_echos,Errors_dfe)
title('dfe')
xlabel('TI used')
ylabel('dfe')

subplot(4,2,6)
plot(1:use_nr_echos,Errors_adjrsquare)
title('adjrsquare')
xlabel('TI used')
ylabel('adjrsquare')
axis([0 test_TIs 0.1 max(max(Errors_adjrsquare))])

subplot(4,2,7)
plot(1:use_nr_echos,T2_combinations)
title('T2 combinations')
xlabel('TI used')
ylabel('T2 [ms]')






end

