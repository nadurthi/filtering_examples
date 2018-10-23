function [xt,F1,F2,F3,PF,GF] = runFilters(iRun,model,prior,time,measAvail,pf,inc,myfilter,domain)
%
% Run all the 4 filters for one particular MC run. Filters:
%   1. no update when measurements are available
%   2. Bayesian classical update
%   3. Discrete update based on CKE
%   4. Particle Filter
%   5. Grid-based Filter
%
% for each Fx/PF there will be saved:
%   Fx.x{k}                     - estimate first moment (mean) - over time
%   Fx.P{k}                     - estimate second moment (covariance) - over time
%   Fx.w{k}                     - weights for gaussian components - over time
%   Fx.ll{k}                    - log-likelihood of the particles - over time
%   Fx.ISD{k}                   - integral square difference between filter Fx and the Grid based filter - over time
%
% Gabriel Terejanu (terejanu@buffalo.edu)

%% ------------------------------------------------------------------------
% init
%--------------------------------------------------------------------------
path(path, genpath(inc));
opt = odeset('reltol',1e-8,'abstol',1e-8);

fprintf('- Run # %d\n',iRun);

iRun = 9;
      
%% ------------------------------------------------------------------------
% create truth & measurement
%--------------------------------------------------------------------------
fprintf('  - create truth & measurements\n');
xt = cell(1, time.nSteps);
y_meas = cell(1, time.nSteps);

% select a gaussian component
rand('state',iRun*100);
gs_sel = 1;
u = rand;
u_total = 0;
for j = 1 : prior.n
    if ((u >= u_total) && (u < u_total + prior.weig(j)))
        gs_sel = j;
        break;
    end
    u_total = u_total + prior.weig(j);
end

% draw a sample from the chosen gaussian component
randn('state',iRun*100);
xt{1} = prior.mu{gs_sel} + chol(prior.sig{gs_sel})' * randn(model.fn,1); 
y_meas{1} = feval(model.hx,xt{1}) + chol(model.R)' * randn(model.hn,1);

% get the trajectory of the sample over time
for k = 2 : time.nSteps
    y = feval(model.fx, time.dt, k, xt{k-1});
    xt{k} = y + model.sQ * randn(model.fn,1);
    y_meas{k} = feval(model.hx, xt{k}) + chol(model.R)' * randn(model.hn,1);
end

%% ------------------------------------------------------------------------
% GF - for Grid-based filters
%--------------------------------------------------------------------------
GF_px = zeros(size(domain.x));
for i = 1 : prior.n
    GF_px = GF_px + prior.weig(i) .* normpdf(domain.x,prior.mu{i},sqrt(prior.sig{i}));
end;

%% ------------------------------------------------------------------------
% PF - for particle filters
%--------------------------------------------------------------------------
fprintf('  - create initial particles\n');

% generate i.c. samples
u_pf = rand(1,pf.no_particles);
tu_pf = 0;
X_pf = [];
for i = 1 : prior.n
    % for each gaussian component draw samples dictated by its weight
    % magnitude
    ind = find((u_pf >= tu_pf) & (u_pf < tu_pf + prior.weig(i)));
    Xtmp_pf = repmat(prior.mu{i},1,length(ind)) + chol(prior.sig{i})' * randn(model.fn,length(ind));
    
    % collect all the samples and advance
    X_pf = [X_pf Xtmp_pf];
    tu_pf = tu_pf + prior.weig(i);
end
w_pf = ones(1, pf.no_particles) / pf.no_particles;

%% ------------------------------------------------------------------------
% compute initial estimates
%--------------------------------------------------------------------------
F1.w{1} = prior.weig; 
F2.w{1} = prior.weig; 
F3.w{1} = prior.weig; 

[F1.x{1},F1.P{1},F1.ll{1},F1.ISD{1}] =  getGSdata(F1.w{1},prior.mu,prior.sig,X_pf,w_pf,domain.x,GF_px);
[F2.x{1},F2.P{1},F2.ll{1},F2.ISD{1}] =  getGSdata(F2.w{1},prior.mu,prior.sig,X_pf,w_pf,domain.x,GF_px);
[F3.x{1},F3.P{1},F3.ll{1},F3.ISD{1}] =  getGSdata(F3.w{1},prior.mu,prior.sig,X_pf,w_pf,domain.x,GF_px);
[PF.x{1},PF.P{1}] = getPFdata(X_pf,w_pf);

GF.x{1} = trapz(domain.x, domain.x .* GF_px);
GF.P{1} = trapz(domain.x, (domain.x - GF.x{1}).^2 .* GF_px);

%% ------------------------------------------------------------------------
% init filters
%--------------------------------------------------------------------------
mu = prior.mu;
sig = prior.sig;

%% ------------------------------------------------------------------------
% PLOT - Create the subplot
%--------------------------------------------------------------------------
if pf.plot,
    total_rows = 3;
    str_filters{1} = 'Grid Filter';
    str_filters{2} = 'Classic GSF';
    str_filters{3} = 'Adaptive GSF';

    total_columns = (time.nSteps-1)/2+1;
    h=subplot1(total_rows,total_columns,'Gap',[0 0]);
    for i = 1 : total_rows
        for j = 1 : total_columns
            set(h(i,j),'XTick',[]);
            if (j == 1)
                subplot1((i-1)*total_columns+1);
                ylabel(str_filters{i});
            end;
            if (i == total_rows)
                subplot1((total_rows-1)*total_columns + j);
                xlabel(num2str(j-1));
            end;
        end;
    end;
end;

k = 1;

subplot1(k); 

    ind_x = find((domain.x >= domain.plot_lim(1)) & (domain.x <= domain.plot_lim(2)));

     %--- grid based filter --------------------------------------
    patch(GF_px(ind_x), domain.x(ind_x), [.5 .5 .5]);

    %--- gaussian sum filters -------------------------------------
    mmZ_tmp1 = zeros(size(domain.x));
    mmZ_tmp2 = zeros(size(domain.x));
    mmZ_tmp3 = zeros(size(domain.x));
    for ii = 1 : prior.n
        tmp_mu = mu{ii};
        tmp_sig = sig{ii};

        % get likelihood only on the jth direction
        mmZ_tmp1 = mmZ_tmp1 + F1.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
        mmZ_tmp2 = mmZ_tmp2 + F2.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
        mmZ_tmp3 = mmZ_tmp3 + F3.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
    end     
    m = max([max(GF_px) max(mmZ_tmp2) max(mmZ_tmp3)]);
    xlim([0 1.05*m]); ylim(domain.plot_lim);

subplot1(total_columns + k);
    patch(mmZ_tmp2(ind_x), domain.x(ind_x), [.5 .5 .5]); 
    xlim([0 1.05*m]); ylim(domain.plot_lim);

subplot1(2*total_columns + k);
    patch(mmZ_tmp3(ind_x), domain.x(ind_x), [.5 .5 .5]);
    xlim([0 1.05*m]); ylim(domain.plot_lim);

%% ------------------------------------------------------------------------
% RUN FILTERS
%--------------------------------------------------------------------------   
for k = 2 : time.nSteps
    fprintf('  - run filters tStep %d / %d\n', k, time.nSteps);    

    % save data 
    mu_old = mu;
    sig_old = sig;
    
%% ------------------------------------------------------------------------
% GF - propagate pdf using grid based methods
%--------------------------------------------------------------------------    
    fprintf('     - Grid-based Filter\n'); 
    
    px_new = zeros(size(GF_px));
    for i = 1 : domain.no_points
        fx1D = feval(model.fx, time.dt, k, domain.x);
        tmp1 = normpdf(domain.x(i) - fx1D, 0, model.sQ);
        px_new(i) = sum(tmp1 .* GF_px);
        if (measAvail(k) == 1)
            tmp2(i) = normpdf(y_meas{k} - feval(model.hx, domain.x(i)), 0, sqrt(model.R));
        end;
    end;
    GF_px = px_new ./ trapz(domain.x, px_new);  

     % GF - measurement update    
    if (measAvail(k) == 1)
        GF_px = GF_px .* tmp2;
        GF_px = GF_px ./ trapz(domain.x, GF_px);
    end;
    
    % save data
    GF.x{k} = trapz(domain.x, domain.x .* GF_px);
    GF.P{k} = trapz(domain.x, (domain.x - GF.x{k}).^2 .* GF_px);
    
%% ------------------------------------------------------------------------
% PF - propagate particles
%--------------------------------------------------------------------------
    fprintf('     - Particle Filter\n'); 
    
    % PF - time update
    X_pf = PF_time_update(X_pf, model, time, k);
    
    % PF - measurement update    
    if (measAvail(k) == 1)
        [X_pf, w_pf] = PF_meas_update(X_pf, w_pf, model, pf, y_meas{k});
    end; 
                
    % PF - compute estimates
    [PF.x{k},PF.P{k}] = getPFdata(X_pf, w_pf);
     
%% ------------------------------------------------------------------------
% propagate and update each gaussian component using Kalman Filter
%--------------------------------------------------------------------------
    fprintf('     - Kalman Filter\n');
    
    for j = 1 : prior.n          
        % time update
        [mu{j}, sig{j}] = feval([myfilter.name '_time_update'], model, mu{j}, sig{j}, time, k);
        
        % do update
        if (measAvail(k) == 1)               
            [mu_up{j}, sig_up{j}, z_mu{j}, z_sig{j}] = feval([myfilter.name '_measurement_update'], model, mu{j}, sig{j}, y_meas{k});         
        else
            mu_up{j} = mu{j};
            sig_up{j} = sig{j};                
        end             
    end   
    
%% ------------------------------------------------------------------------
% GS1 - no weight update
%--------------------------------------------------------------------------
    fprintf('     - Weights GS1\n');
    F1.w{k} = F1.w{k-1};    

%% ------------------------------------------------------------------------
% GS2 - Alspach - classic weight update
%--------------------------------------------------------------------------
    fprintf('     - Weights GS2\n');
    if (measAvail(k) == 1)
        for j = 1 : prior.n
            F2.w{k}(j) = F2.w{k-1}(j)*getLikelihood(y_meas{k} - z_mu{j}, z_sig{j} + model.R);
        end
        if (sum(F2.w{k}) > 0)
            F2.w{k} = F2.w{k} ./ sum(F2.w{k});
        else
            F2.w{k} = F2.w{k-1};
        end;
    else
        F2.w{k} = F2.w{k-1}; 
    end 

%% ------------------------------------------------------------------------
% GS3 - Discrete Update - using CKE
%--------------------------------------------------------------------------   
    fprintf('     - Weights GS3\n');

    switch myfilter.integration_method
        case 'IUT'
            % compute matrix A & B - Unscented Transformation
        case 'GQ_perComp'
            % compute matrix A & B - using Gaussian Quadrature with quadrature
            % points for each gaussian component
        case 'GQ_all'
            % compute matrix A & B - using Gaussian Quadrature with
            % quadrature points that cover the entire domain
            [A,B] = ComputeCKE_QP_GQ_all(F3.w{k-1}, mu, sig, mu_old, sig_old, model, time, myfilter, k);
        otherwise
            error('Integration method unknown');
    end; 

    % get the weights   
    cc = ones(prior.n,1);
    xx = sdpvar(prior.n,1);
    xx_old = reshape(F3.w{k-1},prior.n,1);
    optyalmip = sdpsettings('solver','sedumi','verbose',0,'usex0',1);
    const = set(cc'*xx == 1) + set(xx >= 0);
    solvesdp(const,1/2*xx'*A*xx - xx'*B*xx_old,optyalmip);    
    F3.w{k} = double(xx);
    saveF3w = F3.w{k};

    % do the classic weight update GS2
    if (measAvail(k) == 1)
        for j = 1 : prior.n
            F3.w{k}(j) = F3.w{k}(j)*getLikelihood(y_meas{k} - z_mu{j}, z_sig{j} + model.R);
        end
        if (sum(F3.w{k}) > 0)
            F3.w{k} = F3.w{k} ./ sum(F3.w{k});
        else
            F3.w{k} = saveF3w;
        end;
    end   
    
%% ------------------------------------------------------------------------
% Store estimates
%--------------------------------------------------------------------------
    mu = mu_up;
    sig = sig_up;

%% ------------------------------------------------------------------------
% Compute estimates
%--------------------------------------------------------------------------
    [F1.x{k},F1.P{k},F1.ll{k},F1.ISD{k}] =  getGSdata(F1.w{k},mu,sig,X_pf,w_pf,domain.x,GF_px);
    [F2.x{k},F2.P{k},F2.ll{k},F2.ISD{k}] =  getGSdata(F2.w{k},mu,sig,X_pf,w_pf,domain.x,GF_px);
    [F3.x{k},F3.P{k},F3.ll{k},F3.ISD{k}] =  getGSdata(F3.w{k},mu,sig,X_pf,w_pf,domain.x,GF_px);    
    
%% ------------------------------------------------------------------------
% Plot 
%--------------------------------------------------------------------------

%     subplot(511); 
%          %--- grid based filter --------------------------------------
%         plot(domain.x, GF_px);
%         title(['tstep ' num2str(k) ' / ' num2str(time.nSteps) ' - GF']);
%         
%     subplot(512);
%         %--- particle filter --------------------------------------
%         I = myresampling(w_pf);
%         I = round(I);
%         tmp_X_pf = zeros(1, pf.no_particles);
%         for jj = 1 : pf.no_particles
%             tmp_X_pf(:,jj) = X_pf(:,I(jj));
%         end;     
% 
%         % get only the jth dimension
%         [mTN,mTX] = hist(tmp_X_pf, pf.no_bins);        
%     
%         [tmpTX,tmpI] = unique(mTX);
%         mmZ_tmp = interp1(tmpTX,mTN(tmpI),domain.x);
%         ind_tmp = find(isnan(mmZ_tmp));
%         mmZ_tmp(ind_tmp) = 0;
%         
%         plot(domain.x, mmZ_tmp);
%         title('PF');
%         
%     subplot(513);
%         %--- gaussian sum filters -------------------------------------
%         mmZ_tmp1 = zeros(size(domain.x));
%         mmZ_tmp2 = zeros(size(domain.x));
%         mmZ_tmp3 = zeros(size(domain.x));
%         for ii = 1 : prior.n
%             tmp_mu = mu{ii};
%             tmp_sig = sig{ii};
% 
%             % get likelihood only on the jth direction
%             mmZ_tmp1 = mmZ_tmp1 + F1.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
%             mmZ_tmp2 = mmZ_tmp2 + F2.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
%             mmZ_tmp3 = mmZ_tmp3 + F3.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
%         end      
%         
%         plot(domain.x, mmZ_tmp1);
%         title('GSF - no update');
%         
%     subplot(514);
%         plot(domain.x, mmZ_tmp2);
%         title('GSF - classic update');        
% 
%     subplot(515);
%         plot(domain.x, mmZ_tmp3);
%         title('AGSF');          
%                 
%     drawnow;
%     
%     mu
% %     pause;
%     F2.w{k}
%     F3.w{k}'

%% ------------------------------------------------------------------------
% Plot 
%--------------------------------------------------------------------------
    if ((pf.plot) & (floor(k/2) ~= k/2)),
        
        kk = (k-1)/2+1;
        
        ind_x = find((domain.x >= domain.plot_lim(1)) & (domain.x <= domain.plot_lim(2)));
        
        subplot1(kk); 
             %--- grid based filter --------------------------------------
            patch(GF_px(ind_x), domain.x(ind_x), [.5 .5 .5]);

            %--- gaussian sum filters -------------------------------------
            mmZ_tmp1 = zeros(size(domain.x));
            mmZ_tmp2 = zeros(size(domain.x));
            mmZ_tmp3 = zeros(size(domain.x));
            for ii = 1 : prior.n
                tmp_mu = mu{ii};
                tmp_sig = sig{ii};

                % get likelihood only on the jth direction
                mmZ_tmp1 = mmZ_tmp1 + F1.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
                mmZ_tmp2 = mmZ_tmp2 + F2.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
                mmZ_tmp3 = mmZ_tmp3 + F3.w{k}(ii)*getLikelihoodMany(tmp_mu-domain.x, tmp_sig);
            end      
            
            m = max([max(GF_px) max(mmZ_tmp2) max(mmZ_tmp3)]);
            xlim([0 1.05*m]); ylim(domain.plot_lim);

        subplot1(total_columns + kk);
            patch(mmZ_tmp2(ind_x), domain.x(ind_x), [.5 .5 .5]); 
            xlim([0 1.05*m]); ylim(domain.plot_lim);

        subplot1(2*total_columns + kk);
            patch(mmZ_tmp3(ind_x), domain.x(ind_x), [.5 .5 .5]);
            xlim([0 1.05*m]); ylim(domain.plot_lim);

        drawnow;
    end;

end