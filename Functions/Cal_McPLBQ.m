function [Mc] = Cal_McPLBQ(mag) 
    Mc = 4.0;
    Mmin = 2.0;% floor(min(mag)); 
    Mmax = Mc;
    % Mmax = floor(max(mag));
    mc = [Mmin:0.1:Mmax]';
    
    a=[];b=[];dm=0.1;nsim=1000;pval=0.01;%0.05
    
    % 求第一个满足h = 0的mc
    % ind = -1 ;
    for n = 1:length(mc) % 下标
        [stat,h,dobs] = ClausetTest(mag,a,b,mc(n,1),dm,nsim,pval);
        mc(n,2) = h;
%         disp(['mc = ',num2str(mc(n)),', stat metrix = ',num2str(stat)]);
	    if h == 0
		    Mc = mc(n,1);
            break
        end
    end
    
%     if ind == -1  % 所有mc的h都为1的情况
% %         fprintf('No Mc is satisfied!!!');
%         return
%     end
%     
%     % 在更小的区域找Mc
%     if ind == 1
%         Mc = Mmin; % 即Mc = mc(1, 1);
%     else
%         MC = [mc(ind-1,1):dm:mc(ind,1)]';
%         Mc = mc(ind,1);
%         for N = 1:length(MC)
%             [stat,h,~] = ClausetTest(mag,a,b,MC(N,1),dm,nsim,pval);
%             MC(N,2) = h;
% %             disp(['MC = ',num2str(MC(N)),', stat metrix = ',num2str(stat)]);
% 		    if h == 0
% 			    Mc = MC(N,1);
%                 break
%             end
%         end
%     end        
end
