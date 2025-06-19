function [ mc, mt, b ] = Cal_mc_mt_New( MAG )

            % Calculating mc
            disp('Calculating mc ... ')
            mc = Cal_McPLBQ(MAG);
            
            % Complete catalog
            Mag = MAG(MAG >= mc);
            
            % Calculating beta
            disp('Calculating beta ... ')
            [b,~] = Tinti_aval(Mag,mc,0.1,0); 
            bval_temp = b;
            
            % Calculating mt
            disp('Calculating mt ... ')
            [res] = Cal_fmd(Mag,0.1,mc); res.mi = round(res.mi,1);
            ind_mc = find(res.mi == round(mc,1));
            mt_loop = res.mi(ind_mc:end);
            for t = 1:size(mt_loop,2)
                m = Mag( Mag > mc &  Mag<=mt_loop(t));
                N = size(m,1);
                S = 0;
                for i = 1:N
                    EM_CDF = 2*i / (2*N+1);
                    F_CDF  = -(10.^(-bval_temp .* m(i)) - 10.^(-bval_temp .* mc)) ./ (10.^(-bval_temp .* mc) - 10.^(-bval_temp .* max(Mag)));
                    S_temp = 0;
                    S_temp = EM_CDF * log(F_CDF) + (1 - EM_CDF) * log(1 - F_CDF);
                    if isempty(S_temp)
                        break
                    end
                    S = S+S_temp;
                end
                D(t) = -N-2*S;
            end
            D(end) = D(end-1);  D(1) = D(2); % Excluding 0 and inf
            [~,ind_min] = min(D);
            mt = res.mi(ind_min);

            % Make sure
            if size(Mag(Mag>=mc & Mag < mt),1) >= 100
                mt = mt;
            else
                mt = res.mi(end);
            end

            % Calculating mc
            disp('Calculating mc ... ')
            mag_temp = []; mag_temp = MAG(MAG < mt);
            mc = Cal_McPLBQ(mag_temp);
                
            % Re-calculating beta
            disp('Calculating b ... ')
            mag = []; mag = MAG(MAG >= mc & MAG < mt);
            [b,~] = Tinti_aval(mag,mc,0.1,0); 

            mc
            mt
        
end