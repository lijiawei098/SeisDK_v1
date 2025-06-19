function [Sc,testStatistics] = Cal_CriticalValue(loop, Mag, alpha, setting, Statistics)
    
    % Parameters
    mbin = 0.1;  mc = setting.mc;  mt = setting.mt;  b = setting.b; 
    setting.loop = loop;

    % Number of simulations
    numSimulations = 10000; % 10000;
    
    % Initialize array to store test statistics
    testStatistics = zeros(1, numSimulations);

    % Number of samples per simulation
    Num_sam = length(Mag);
    
    % Perform Monte Carlo simulations
    % figure; hold on;
    for i = 1:numSimulations
        % Generate random samples from the Gutenberg-Richter PDF
        m = [];  
        m = GR_truncated_simulator_ver2(Num_sam,8.5,mc-mbin/2,b);
        % m = randomGutenbergRichter(beta, mc-mbin/2, Num_sam);  m =m';
        % m = round(m, 1);
        m = BinMags(m, mc, mbin);
        m = sort(m);
        
        % Plot simulating power law
        % Plot_SimGR( m, path, setting, 0 );
        
        % Calculate test statistic (e.g., max-robust-sum test statistic, MRS)
        testStatistics(i) = TestStatistic(m,m(1:setting.Num),loop,Statistics); 
        
    end
    % hold off;
    % title(['No.',num2str(loop),',$\beta e^{-\beta (m - m_{\mathrm{c}})}, \beta = ', num2str(round(setting.beta, 2)), '$'], 'Interpreter', 'latex');
    
    % Determine the critical value as the (1 - alpha) percentile of the test statistics distribution
    Sc = prctile(testStatistics, 100 * (1 - alpha));
end




