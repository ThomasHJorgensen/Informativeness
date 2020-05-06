% This static class contains all functions used to estimate the model.
% Functions are devided into groups based on their purpose.
% within each group, functions appear in heirarical order. The most high
% level functions appear first and the functions used within a function
% appear below it.
classdef fun
    methods (Static)

        % Model
        function [Opt_h,Opt_w]  = calc_ret_combi_mex(data,par)
            % the discounted sum of utility from all combinations of
            % retirement ages of husband and wife.
            % INPUT:
            %   data:       structure with sub-fields "h" and "w" for husband
            %               and wife, respectively.
            %   par:        structure with all model parameters and settings
            % OUTPUT:
            %   combi:      N-by-S-by-NumRetChoice-by-NumRetChoise array with
            %               weighted discounted sum of utility from the
            %               combination of retirement ages
            
            if par.restricted == 1
                par.joint_w = par.joint_h;
            end
            
            % Pre-calculate parts that are independent of the retirement choise in the instantaneous utility this is now just a vector in par-struct
            par.disc_vec = (par.disc.^(0:(par.age_max-par.age_min+1)));
            
            data.h.Xbeta = data.h.X*par.beta_h';
            data.w.Xbeta = data.w.X*par.beta_w';

            T       = [par.var_h, 0;
                       par.cov  ,  par.var_w]; % par.var_h = 1; 
            Omega   = T*T';                   % ensures positive defenite
    
            sqrtOmega      = sqrtm(Omega); 
            data.h.e_draws = sqrtOmega(1,1)*data.h.draws + sqrtOmega(1,2)*data.w.draws;
            data.w.e_draws = sqrtOmega(2,2)*data.w.draws + sqrtOmega(2,1)*data.h.draws;
            
            % loop through all combinations of retirement ages to calculate
            % the assocoated value of that choice
            par.Nobs        = data.Nobs;
            par.Nsim        = data.S;
            par.NumRetGrid  = numel(par.RetTimeGrid);
            par.NumAge      = numel(par.age_h);

            [Opt_h,Opt_w]   = mexCombi(par,data.h,data.w);
            
        end
        
        % Misc: setup and loading data                              
        function par  = setup(data)
            % Sets up the different parameters and settings related to
            % estimation of the model.
            % INPUT:
            %   data:   structure containing the loaded BHPS data. Mainly
            %           used to determine dimensions of the data (observations and
            %           explanatory variables)
            % OUTPUT:
            %   par:    structure containing all model parameters and
            %           settings. This struct is passed to all subsequent
            %           functions.

            % Estimation settings
            par.IndirectInference = 0; % 1->II, 0->SML 
            par.do_particleswarm= 1;
            par.restricted      = 1; % restric the value of joint leisure to be the same across couples
            par.PRINT           = 1;
            par.NumStart        = 4; % number of times the global search sequence of optimizers are called
            
            options             = optimset('MaxIter',5000,'TolX',1.0e-6,'TolFun',1.0e-6,'Display','iter');
            options_swarm       = optimoptions('particleswarm','Display','iter','SwarmSize',80,'MaxIter',500);
            par.options_swarm   = options_swarm;
            par.options         = options;

            % Potential retirement times/ages
            par.ret_min = 50;
            par.ret_max = 70;
            par.RetTimeGrid = par.ret_min:par.ret_max;
            
            % Age of certain death
            par.age_max = 80;
            par.age_min = min(min(data.h.Age) , min(data.w.Age));
            
            % parameters on explanatory variables
            par.beta_h = zeros(size(data.h.X,2),1)';
            par.beta_w = zeros(size(data.w.X,2),1)';
            
            % Age-dummy variables
            par.RetAge_cons  = 0*ones(1,2);
            par.RetAge_trend = 0.0*ones(1,2);
            par.RetAge50     = 0*ones(1,2);
            par.RetAge55     = 0*[.35 .3];%3*ones(1,2);
            par.RetAge60     = 0*[.21 .3];%.3*ones(1,2);
            par.RetAge65     = 0*[.68 .05];%.2*ones(1,2);
            par.RetAge70     = 0*ones(1,2);

            par.RetAgeDum_fun = @(Age,m,par)   par.RetAge_cons(m) + par.RetAge_trend(m)*(Age-65) ...
                                          + par.RetAge50(m)*(Age>=50) ...
                                          + par.RetAge55(m)*(Age>=55) ...
                                          + par.RetAge60(m)*(Age>=60) ...
                                          + par.RetAge65(m)*(Age>=65) ...
                                          + par.RetAge70(m)*(Age>=70);
            
            par.age_grid{1}     = 50:54;
            par.age_grid{2}     = 55;
            par.age_grid{3}     = 56:59;
            par.age_grid{4}     = 60;
            par.age_grid{5}     = 61:64;
            par.age_grid{6}     = 65;
            %par.age_grid{7}     = 66:70;   % this is the reference group.
            %Otherwise there will be perfect colliniarity in the moments
            
                                      
            par.joint_h  = 0.02;%0.2; % value of joint retirement
            par.joint_w  = par.joint_h; % Value of joint retirement
            par.disc     = .97; % Annual discount factor
            par.SPA_w    = .1; % The value of being retired once older than the SPA age
            par.weight_h = 1; % The relative weight on the husband in the household utility function
            
            % Shock variance and covariance
            par.var_h    = 1; % normalization
            par.var_w    = 0.7;
            par.cov      = 0.27;

        end
        function data = load_data(seed,S,dest)
            % loads BHPS data from .txt files
            % INPUT:
            %   seed:   the seed used to simulate draws of the random taste
            %           shock
            %   S:      the number of random taste shocks to simulate
            %   dest:   the denstination folder of the data
            % OUTPUT:
            %   data:   structure with sub-fields "h" and "w" for husband
            %           and wife, respectively.
            if nargin<3
                dest    = 'data';
            end
            
            members = {'h','w'};
            for m=1:numel(members)
                data.(members{m}).ExpRetAge     = importdata([dest '\ExpRetAge' num2str(m) '.txt']);
                data.(members{m}).Age           = importdata([dest '\Age' num2str(m) '.txt']);
                data.(members{m}).Age2          = data.(members{m}).Age.^2;
                data.(members{m}).Retired       = importdata([dest '\retired' num2str(m) '.txt']);
                data.(members{m}).SPA           = importdata([dest '\SPA' num2str(m) '.txt']);
                data.(members{m}).White         = importdata([dest '\White' num2str(m) '.txt']);
                data.(members{m}).PoorHealth    = importdata([dest '\gp_many' num2str(m) '.txt']);
                data.(members{m}).ExpWorseHealth= importdata([dest '\ExpWorseHealth' num2str(m) '.txt']);
                data.(members{m}).Highskilled   = importdata([dest '\Highskilled' num2str(m) '.txt']);
                % Financial variables
                data.(members{m}).Income        = importdata([dest '\Income' num2str(m) '.txt']);
                data.(members{m}).lIncome       = log(data.(members{m}).Income);
                data.(members{m}).Wealth        = importdata([dest '\lWealth' num2str(m) '.txt']);
                data.(members{m}).pensPPP       = importdata([dest '\pensPPP' num2str(m) '.txt']);
                data.(members{m}).pensEPS       = importdata([dest '\pensEPS' num2str(m) '.txt']);
                
                data.(members{m}).cohort       = importdata([dest '\cohort' num2str(m) '.txt']);

            end
            % remove NaN (=-2)
            vars = {'ExpRetAge','Retired','Highskilled','Age','Age2','PoorHealth','ExpWorseHealth','Income','White','pensPPP','pensEPS','Wealth','Income','pensPPP','pensEPS','cohort'};
            for m=1:numel(members)
                for v=1:numel(vars)
                    var_nan     =  data.(members{m}).(vars{v})<(-2+1.0e-6);
                    data.(members{m}).(vars{v})(var_nan) = NaN;
                end
            end

            data.w.SPA_above60  = data.w.SPA - 60;
            data.w.SPA_5155     = data.w.cohort>=1951 & data.w.cohort<1955 + 0*data.w.cohort; 
            data.w.SPA_56       = data.w.cohort>=1955 + 0*data.w.cohort;
            
            % TJ: June '19: re-scale cohort variable
            cohort_base = 1955;
            for m=1:numel(members)
                data.(members{m}).cohort       = data.(members{m}).cohort - cohort_base;
            end
            
            % add draws
            data.Nobs = size(data.h.Age,1);
            data.S    = S;
            
            rng(seed);
            data.h.draws    = randn(data.Nobs,data.S);
            data.w.draws    = randn(data.Nobs,data.S); 

        end
        function data = explanatory_vars(data,ExplVar,ExplVar_excl,EstPar)
            % construct explanatory variables in matrix
            % INPUT:
            %   data:           structure with sub-fields "h" and "w" for husband
            %                   and wife, respectively.
            %   ExplVar:        cell-array with names of parameters to be included
            %                   as explanatiry variables in the utility
            %   ExplVar_excl:   cell-array with names of parameters to be included
            %                   exclusively in the auxiliary model. Only
            %                   used in II
            % OUTPUT:
            %   data:           updated data-structure with the field X for
            %                   both sub-fields h and w. Potentially also
            %                   another X containing more variables in
            %                   data."member".aux.X.
            
            % Loop over household members
            members = {'h','w'};
            for m=1:numel(members)
                X = NaN(data.Nobs,numel(ExplVar.(members{m}).Own)+numel(ExplVar.(members{m}).Spouse));
                for i=1:numel(ExplVar.(members{m}).Own)
                    % change "-2" to NaN (missings). This is done in STATA because
                    % missings are not exported into txt-files....
                    VAR          = data.(members{m}).(ExplVar.(members{m}).Own{i});
%                     VAR(VAR<0) = NaN; % NaN are set to -2 in STATA. All variables are supposted to be positive.

                    % insert the explanatory variable
                    X(:,i) = VAR;
                    
                    names{i,m} =  [ExplVar.(members{m}).Own{i} '_' members{m}];
                end
                for j=1:numel(ExplVar.(members{m}).Spouse)
                    if strcmp(members(m),'h')
                        ms = 2;
                    else
                        ms = 1;
                    end
                    VAR          = data.(members{ms}).(ExplVar.(members{m}).Spouse{j});
%                     VAR(VAR<0) = NaN; % NaN are set to -2 in STATA. All variables are supposted to be positive.
                    X(:,j+i) = VAR;
                    
                    names{j+i,m} = [ ExplVar.(members{m}).Spouse{j} '_' members{m} '_s'];
                end
                data.(members{m}).X = X;
                
                % add variables only entering the auxiliary model
                Z = X(:,1:numel(ExplVar.(members{m}).Own)); % but remove spousal variables: all are included directly
                for i=1:numel(ExplVar_excl.(members{m}))
                    VAR          = data.(members{m}).(ExplVar_excl.(members{m}){i});
%                     VAR(VAR<0) = NaN;
                    Z = [Z, VAR];

                end
                data.(members{m}).aux.X = Z;
                
                % Construct missings in response variable
                data.(members{m}).ExpRetAge(data.(members{m}).ExpRetAge<0) = NaN;
            end
            data.X.names = [names(:,1) ; names(:,2)];

            % Names:
            clear names;
            j=1;
            for i1=1:numel(EstPar)
                if strcmp(EstPar{i1},'beta_h')
                    beta_names=data.X.names;
                    num_now = numel(beta_names);
                    for i2 = 1:num_now
                        names{j+i2-1} = beta_names{i2};
                    end
                elseif strcmp(EstPar{i1},'beta_w')
                    % do nothing because all are included above
                    num_now = 0;
                elseif strcmp(EstPar{i1},'joint_h')  || strcmp(EstPar{i1},'joint_w') || strcmp(EstPar{i1},'SPA_w') || strcmp(EstPar{i1},'var_w') || strcmp(EstPar{i1},'cov')
                    num_now = 1;
                    names{j} = EstPar{i1};
                else
                    num_now = 2; % are parameters
                    for i2 = 1:num_now
                        names{j+i2-1} = [EstPar{i1} '_' members{i2}];
                    end
                end
                j=j+num_now;
            end
            
            data.names = names;
        end
        function [bsCov, bsDiagVar] = bootstrap_cov(b_aux,data,par)
            
            N = data.Nobs;
%             NumMoments  = numel(fun.calc_moments(data,par));
%             bsMoments   = NaN(NumMoments , par.bsNum);
            
            for bs_i=1:200
                rng(9210+bs_i) % Set seed
                
                % Draw random sample from the data. With replacement
                [~,id]  = datasample(ones(N,1),N);
                id      = id';
                
                % Draw the same households used to calculate moments from data
                 members = {'h','w'};
                 for m=1:numel(members)
                     names = fieldnames(data.(members{m}));
                     for n=1:numel(names)
                         
                         if strcmp(names{n},'aux')
                             rdata.(members{m}).aux.X    = data.(members{m}).aux.X(id,:);
                         else
                            rdata.(members{m}).(names{n})   = data.(members{m}).(names{n})(id,:);
                         end
                         
                     end
                 end
                rdata.Nobs = data.Nobs;
                rdata.S    = data.S;
                bsMoments(:,bs_i) = fun.aux_moments(b_aux,rdata,par,'data');
            end
            
            bsCov = cov(bsMoments');
            bsVar = diag(bsCov);
            bsDiagVar = diag(bsVar);

        end
        
        % Estimation
        function par = estimate(EstPar,data,par,W,b_aux)
            % Global estimation routine applying succesive numerical
            % optimizers to estimate the model parameters.
            % INPUT:
            %   EstPar: cell-array with names of parameters to be estimated
            %   data:   structure with sub-fields "h" and "w" for husband
            %           and wife, respectively.
            %   par:    structure with all model parameters and settings
            %   W:      weight matrix, only used for indirect inference (II)
            %   b_aux:  auxiliary parameter estimates, only used for II
            % OUTPUT:
            %   par:    updated par-structure with estimates and updated
            %           values of estimated parameters
            
            % 0. store the number of simulations for later use
            S_orig = data.S;
            data.S = 100;
            par.S = data.S;
            
            % 1. Construct initial values
            orig_PRINT = par.PRINT;
            par.PRINT  = 0;

            % initial starting values: ordered probit
            EstAgeDum = {'RetAge_cons','RetAge_trend','RetAge55','RetAge60','RetAge65'};
            b_init    = fun.oprobit_est(data,par,EstAgeDum);
            
            % insert estimates into par-struct
            num_beta.h = size(b_init.h,1) - numel(EstAgeDum);
            num_beta.w = size(b_init.w,1) - numel(EstAgeDum);
            par_init = par;
            
            % explanatory variables
            par_init.beta_h(1:num_beta.h) = b_init.h(1:num_beta.h);
            par_init.beta_w(1:num_beta.w) = b_init.w(1:num_beta.w);
            
            % age-structure
            members = {'h','w'};
            for m=1:numel(members)
                for p=1:numel(EstAgeDum)
                    par_init.(EstAgeDum{p})(:,m) = b_init.(members{m})(num_beta.(members{m})+p);
                end
            end
            % set the start value for the RetAge60.
            par_init.RetAge60(2) = 0.6;
            
            par       = par_init;
            par.PRINT = orig_PRINT;
            
            theta0 = NaN(1,numel(EstPar));
            i_par  = 1;
            for p=1:numel(EstPar)
                num_now                        = numel(par.(EstPar{p}));
                theta0(i_par:i_par+num_now-1)  = par.(EstPar{p});           % set the initial values of the values in the par-struct as initial guess
                i_par                          = i_par + num_now;
            end
            
            % 2. bounds
            UB = theta0 + abs(theta0)*1.2 + (theta0==0)*.8;
            LB = theta0 - abs(theta0)*1.2 - (theta0==0)*.5;
            
            i_par  = 1;
            for p=1:numel(EstPar)
                num_now                        = numel(par.(EstPar{p}));
                if strcmp(EstPar{p},'RetAge_trend')
                    LB(i_par:i_par+num_now-1) = 0;
                elseif strcmp(EstPar{p},'var_w')
                    LB(i_par:i_par+num_now-1) = 0.1;
                    UB(i_par:i_par+num_now-1) = 1.5;
                elseif strcmp(EstPar{p},'cov')
                    LB(i_par:i_par+num_now-1) = 0.0;
                    UB(i_par:i_par+num_now-1) = 1;
                elseif strcmp(EstPar{p},'joint_h')
                    LB(i_par:i_par+num_now-1) = 0.0;
                    UB(i_par:i_par+num_now-1) = .3;  
                elseif strcmp(EstPar{p},'SPA_w')
                    LB(i_par:i_par+num_now-1) = 0.0;
                    UB(i_par:i_par+num_now-1) = 0.7; 
                elseif strcmp(EstPar{p},'RetAge55') || strcmp(EstPar{p},'RetAge60') || strcmp(EstPar{p},'RetAge65')
                    LB(i_par:i_par+num_now-1) = 0.1;     
                end
                
                i_par                          = i_par + num_now;
            end
            
            % 3. Run a series of estimators to get at the global minimum
            obj     = @(theta) fun.ObjFun(theta,EstPar,data,par,W,b_aux);

            % Initialize objects used across the multistarts
            theta0_orig = theta0;
            UB_orig     = UB;
            LB_orig     = LB;
            MIN         = Inf;
            est_all 	= [];
            
            rng(9211); %  set seed
            
            for i=1:par.NumStart % loop over multistarts
                
                
                all(i).UB = UB;
                all(i).LB = LB;
                
                % a. particle swarm
                est_swarm = theta0;
                if par.do_particleswarm==1
                    if i>3
                        par.options_swarm.InitialSwarmMatrix(i-3,:) = theta0; % Include the estimates in the next particle swarm
                    else
                        % start from the original values again
                        par.options_swarm.InitialSwarmMatrix = theta0_orig;
                        LB = LB_orig;
                        UB = UB_orig;
                    end
                    est_swarm = particleswarm(obj,numel(theta0),LB,UB,par.options_swarm);
                end
                
                % b. use those estimates to initialize fminsearch
                [all(i).est,all(i).fval,all(i).exitflag] = fminsearch(obj,est_swarm,par.options);
                est_all  = [est_all;all(i).est];
                
                % c. update the estimate to the current if that improves
                % the objective function over the last run
                if all(i).fval<MIN
                    MIN  = all(i).fval;
                    fval = all(i).fval;
                    est  = all(i).est;
                end
                
                % d. update theta0 and bounds
                theta0   = est;

                    % d.i adjust the tuning parameters based on the past variance of estimates
                    % only efter two estimations using the original bounds
                    if i>3
                        past_est = est_all(max(1,i-5):i,:);
                        var_past = var(past_est);

                    damp = 1.0;
                    d = damp*sqrt(var_past) + .1;

                    % d.ii calculate/update the bounds: based on the best
                    % (est) but restricted to be within the initial supplied bounds
                    UB =  (est+d) + abs(est).*d ;
                    LB =  (est-d) - abs(est).*d ;
                    end
            end
            
            % 3. run estimation using best estimates as start with more simulations
            data.S = S_orig;
            par.S = data.S;
            est_now = est;
            obj     = @(theta) fun.ObjFun(theta,EstPar,data,par,W,b_aux);
            est = fminsearch(obj,est_now,par.options);
            
            % 3. Output the updated par struct
            par         = fun.update_par(est,EstPar,par);
            par.est     = est;
            par.EstPar  = EstPar;
            par.fval    = fval;
            par.all     = all;
            [~,par.momSim] = fun.ObjFun(est,EstPar,data,par,W,b_aux);
        end

        function par = update_par(theta,EstPar,par)
            % This function takes a vector (theta) of values and inserts
            % them into the par-struct with the associated names contained
            % in "EstPar". 
            i_par = 1;
            for i=1:numel(EstPar)
                num_now         = numel(par.(EstPar{i}));                   % allows for vectors in the par-struct
                par.(EstPar{i}) = theta(i_par:i_par+num_now-1);
                i_par           = i_par + num_now;

                % Print current values:
                if par.PRINT==1
                    if num_now>1
                        for j=1:num_now
                            var_now = par.(EstPar{i});
                            fprintf('%s[%d]=%2.4f, ',EstPar{i},j,var_now(j));
                        end
                    else
                        fprintf('%s=%2.4f, ',EstPar{i},par.(EstPar{i}));
                    end
                end
            end
        end
        function [SE,grad,sens,Avar] = sensitivity(theta,EstPar,data,par,W,b_aux)
            
            % Calculate numerical gradient for all parameters
            num_mom  = size(par.CoVar,1);
            num_par  = numel(theta);
            h    = 1.0e-4;
            grad = NaN(num_mom,num_par);
            for i=1:numel(theta)
                var_now      = zeros(size(theta));
                var_now(i)   = 1;
                
                [~,forward]  = fun.ObjFun(theta+h*var_now,EstPar,data,par,W,b_aux);
                [~,backward] = fun.ObjFun(theta-h*var_now,EstPar,data,par,W,b_aux);
                
                grad(:,i)   = (forward-backward)./(2*h); 
            end
            
            % calculate objects re-used below
            GW       = grad'*W;
            GWG      = GW*grad;
            S        = par.CoVar;
            
            % Calculate asymptotic variance and (adjusted) standard errors
            Avar    = (GWG\(GW*S*GW'))/GWG;
            AvarOpt = inv((grad'/S)*grad);
            SE      = sqrt( diag((1+1/data.S)*Avar)./data.Nobs );
            
            % 4. calculate sensitivity measure (NumPar x NumMom) matrix
            sens.M1       = -GWG\GW;
            std_moms      = sqrt(diag(S));  
            sens.M1e      = sens.M1.*repmat(std_moms',num_par,1);
            
            % 5. alternatives from Honoré, Jørgensen and de Paula
            GSi  = grad'/S;
            GSiG = GSi*grad;
            
            sens.M2 = NaN(numel(theta),num_mom);
            sens.M3 = NaN(numel(theta),num_mom);
            sens.M4 = NaN(numel(theta),num_mom);
            sens.M5 = NaN(numel(theta),num_mom);
            sens.M6 = NaN(numel(theta),num_mom);
            
            sens.M2e = NaN(numel(theta),num_mom);
            sens.M3e = NaN(numel(theta),num_mom);
            sens.M4e = NaN(numel(theta),num_mom);
            sens.M5e = NaN(numel(theta),num_mom);
            sens.M6e = NaN(numel(theta),num_mom);
            
            for k = 1:num_mom
                % pick out the kk'th element: Okk
                O      = zeros(num_mom);
                O(k,k) = 1;
                
                M2kk     = (GSiG\(GSi*O*GSi'))/GSiG;          % NumPar-by-NumPar
                M3kk     = (GWG)\(GW*O*GW')/(GWG);%sens.M1*O*sens.M1'; %(Avar/GWG)*(GW*O*GW')*(GWG\Avar);  % NumPar-by-NumPar
                M6kk     =  - GWG\(grad'*O*grad)*Avar ...
                            + GWG\(grad'*O*S*W*grad)/GWG ...
                            + GWG\(grad'*W*S*O*grad)/GWG ...
                            - Avar*(grad'*O*grad)/GWG;  % NumPar-by-NumPar
                
                sens.M2(:,k)  = diag(M2kk); % store only the diagonal: the effect on the variance of a given parameter from a slight change in the variance of the kth moment
                sens.M3(:,k)  = diag(M3kk); % store only the diagonal: the effect on the variance of a given parameter from a slight change in the variance of the kth moment
                sens.M6(:,k)  = diag(M6kk); % store only the diagonal: the effect on the variance of a given parameter from a slight change in the variance of the kth moment
                
                sens.M2e(:,k)  = diag(M2kk)./diag(AvarOpt) * S(k,k); % store only the diagonal: the effect on the variance of a given parameter from a slight change in the variance of the kth moment
                sens.M3e(:,k)  = diag(M3kk)./diag(Avar) * S(k,k); % store only the diagonal: the effect on the variance of a given parameter from a slight change in the variance of the kth moment
                sens.M6e(:,k)  = diag(M6kk)./diag(Avar) * W(k,k); % store only the diagonal: the effect on the variance of a given parameter from a slight change in the variance of the kth moment
                
                
                % remove the kth moment from the weight matrix and
                % calculate the asymptotic variance without this moment
                W_now      = W;
                W_now(k,:) = 0;
                W_now(:,k) = 0;
                
                GW_now   = grad'*W_now;
                GWG_now  = GW_now*grad;
                Avar_now = (GWG_now\(GW_now*S*GW_now'))/GWG_now;
                
                sens.M4(:,k)  = diag(Avar_now) - diag(Avar);
                sens.M4e(:,k) = (diag(Avar_now) - diag(Avar))./diag(Avar);
                
                S_now=S;
                S_now(k,:)=[];
                S_now(:,k)=[];
                grad_now=grad;
                grad_now(k,:)=[];
                AvarOpt_now = inv((grad_now'/S_now)*grad_now);
                sens.M5(:,k)  = diag(AvarOpt_now) - diag(AvarOpt);
                sens.M5e(:,k) =  (diag(AvarOpt_now) - diag(AvarOpt))./ diag(AvarOpt);
            end
            
            
        end

        function [Obj,mom] = ObjFun(theta,EstPar,data,par,W,b_aux)
            % 1. Update par-struct with the values of theta
            par     = fun.update_par(theta,EstPar,par);
            if par.restricted == 1
                par.joint_w = par.joint_h;
            end
            
            % age-dummies are constructed from the RetAgeDum_fun function
            % that is constructed in "fun.setup"
            par.age_h = par.RetAgeDum_fun(par.age_min:par.age_max,1,par); %par.age_h = par.RetAgeDum_fun(par.RetTimeGrid,1);
            par.age_w = par.RetAgeDum_fun(par.age_min:par.age_max,2,par);
            
            % 2. Calculate all possible combinations of retirement ages and
            % optimal planned retirement age for each spouse 
            [ExpRetAge_h,ExpRetAge_w]   = fun.calc_ret_combi_mex(data,par);
            
            % 4. Evaluate the scores at the estimated parameters using data for the simulated data
            mom_s = NaN(size(W,1),data.S);
            sim = data;
            for s=1:data.S
                sim.h.ExpRetAge = ExpRetAge_h(:,s);
                sim.w.ExpRetAge = ExpRetAge_w(:,s);
                
                mom_s(:,s) = fun.aux_moments(b_aux,sim,par,'sim');%fun.aux_scores(b_aux,sim,par);
            end
            mom = nanmean(mom_s,2); % J-by-1

            % 5. calculate objective function
            Obj = mom'*W*mom;
            
            if par.PRINT==1,          fprintf(' -> Obj=%2.4f\n',Obj);            end

        end
        function b_aux = aux_est_data(data,par)
            
            % 2. Estimate residual correlation from regression models for both spouses
            [b_cov_h,b_cov_w,covi,mom_h,mom_w] = fun.aux_res_cov(data);
            
            % 3. Estimate the share expecting retiring at the same calender year
            ret_year_h = data.h.ExpRetAge - data.h.Age; 
            ret_year_w = data.w.ExpRetAge - data.w.Age;
            ret_diff   = ret_year_h-ret_year_w;
            joint_ret  = (abs(ret_diff)<=0) + 0*ret_diff;
            
            % 3. Return the estimates in a struct
            b_aux.h.cov     = b_cov_h;
            b_aux.w.cov     = b_cov_w;
            b_aux.h.mom     = mom_h;
            b_aux.w.mom     = mom_w;
            b_aux.covi      = covi;
            b_aux.joint_ret = joint_ret;
            
            b_aux.within2minus = (ret_diff<0).*(ret_diff>=-2) + 0*ret_diff;
            b_aux.within2plus  = (ret_diff>0).*(ret_diff<=2) + 0*ret_diff;
            

            % include the coefficients from the histogram
            members = {'h','w'};
            for m=1:1:numel(members)
                for j=1:numel(par.age_grid)
                    h                          = sum(data.(members{m}).ExpRetAge==repmat(par.age_grid{j} , size(data.(members{m}).ExpRetAge,1),1) , 2);
                    b_aux.(members{m}).hist(:,j) = h;
                end
            end
            
            
        end 
        function [s_aux,s_aux_var,names] = aux_moments(b_aux,data,par,type)
             
            % This function returns the individual-level scores evaluated
            % at the estimated parameters in b_aux.

            % Covariance from regression-residuals
            [b_h,b_w,covi,mom_h,mom_w] = fun.aux_res_cov(data,b_aux);
            if strcmp(type,'data')
                gradi_cov  = covi; 
                
                ols_h = mom_h(:,1:end);
                ols_w = mom_w(:,1:end);
            else
                gradi_cov  = b_aux.covi - covi; % difference between the data and the simulation moment
                
                ols_h = b_aux.h.mom(:,1:end) - mom_h(:,1:end);
                ols_w = b_aux.w.mom(:,1:end) - mom_w(:,1:end);
            end
            
      
            % the histograms 
            members = {'h','w'};
            for m=1:1:numel(members)

                for j=1:numel(par.age_grid)
                    h                         = sum(data.(members{m}).ExpRetAge==repmat(par.age_grid{j} , size(data.(members{m}).ExpRetAge,1),1) , 2);
                    if strcmp(type,'data')
                    
                        HIST.(members{m}).hist(:,j) = h;
                    else
                        HIST.(members{m}).hist(:,j) = b_aux.(members{m}).hist(:,j) - h  +  0*data.(members{m}).ExpRetAge;
                        
                    end
                end
            end
            
            
            % Share retiring at the same time.
            ret_year_h = data.h.ExpRetAge - data.h.Age;
            ret_year_w = data.w.ExpRetAge - data.w.Age;
            ret_diff   = ret_year_h-ret_year_w;
            joint_ret  = (abs(ret_diff)<=0) + 0*ret_diff;
            
            within2minus = (ret_diff<0).*(ret_diff>=-2) + 0*ret_diff;
            within2plus  = (ret_diff>0).*(ret_diff<=2) + 0*ret_diff;
                
            if strcmp(type,'data')
                gradi_ret  = joint_ret; % differnece between the data and simulation. a vector but like the other moments we could take means already here and not carry around the vectors
                
                gradi_within2minus = within2minus;
                gradi_within2plus  = within2plus;
            else
                
                gradi_ret  = b_aux.joint_ret - joint_ret; % differnece between the data and simulation. a vector but like the other moments we could take means already here and not carry around the vectors  
                
                gradi_within2minus = b_aux.within2minus - within2minus;
                gradi_within2plus  = b_aux.within2plus - within2plus;
            end
            gradi_ret(isnan(joint_ret)) = NaN; % add nans where the model has NaNs due to explanatory variables
             
            gradi_within2minus(isnan(within2minus)) = NaN; % add nans where the model has NaNs due to explanatory variables
            gradi_within2plus(isnan(within2plus)) = NaN; % add nans where the model has NaNs due to explanatory variables
         
            % Stack the moments
            s_aux           = [ols_h , ols_w , HIST.h.hist, HIST.w.hist, gradi_cov ,gradi_within2minus, gradi_within2plus, gradi_ret ]; % N-by-J (J-3)
            s_aux_var       = nancov(s_aux);
            s_aux           = nanmean(s_aux,1)'; %J-by-1
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % return the names of the moments
            if nargout>2
                members = {'h','w'};
                clear nam OP;
                
                for ms=1:numel(members)
                    
                    % 1. Regression
                    for m=1:numel(members)                        
                        j=0;
                        
                        % 1. regression, own
                        if m==1 % constant
                        j=j+1;
                        nam(m).names{j} = ['reg' members{ms} '_cons'];
                        end
                        
                        num = numel(par.ExplVar.(members{m}).Own);
                        for i = 1:num
                            j=j+1;
                            nam(m).names{j} = ['reg' members{ms} '_' par.ExplVar.(members{m}).Own{i} '_' members{m}];
                        end
                        
                        num = numel(par.ExplVar_excl.(members{m}));
                        num = num*(num>0);
                        for i = 1:num
                            j=j+1;
                            nam(m).names{j} = ['reg' members{ms} '_' par.ExplVar_excl.(members{m}){i} '_' m];
                        end
                        
                    end
                    OP(ms).names = [nam(1).names nam(2).names];
                    
                    % 2. Histogram
                    num = numel(par.age_grid);
                    for i=1:num
                        HIST(ms).names{i} = ['H_' members{ms} '_' num2str(min(par.age_grid{i})) '_' num2str(max(par.age_grid{i}))];
                    end
                    
                end
                
                names = [OP(1).names,OP(2).names,HIST(1).names,HIST(2).names,{'Reg_var_h','Reg_var_w','Reg_cov_hw','diff [-2,-1]','diff [1,2]','Joint retirement'}];

            end
        end
        function [b_h,b_w,gradi,mom_h,mom_w]  = aux_res_cov(data,vararg_b_aux)
            % setup data
            xi        = [ones(data.Nobs,1), data.h.aux.X, data.w.aux.X];  % husband and wife information included. Add a constant here because that is not included in the X's
            ret_h     = data.h.ExpRetAge;
            ret_w     = data.w.ExpRetAge;
                
            if nargin<2
                b_h = regress(ret_h,xi);
                b_w = regress(ret_w,xi);
            else
                % use the input'et estimated values
                b_h = vararg_b_aux.h.cov;
                b_w = vararg_b_aux.w.cov;
            end
            
            e_h = ret_h  - xi*b_h;
            e_w = ret_w  - xi*b_w;
            
            % return the covariance matrix
            gradi = [e_h.^2 e_w.^2 e_h.*e_w]; % N-by-3
           
            % and the moment restriction
            mom_h = xi.*repmat(e_h,1,size(xi,2));
            mom_w = xi.*repmat(e_w,1,size(xi,2));
            
        end
        
        % Initial values: Ordered Probit
        function [obj]      = oprobit_obj(theta,data,par,member_str,EstAgeDum)
            LARGE_NUMBER = 10000000000;
            ExpRetAge = data.(member_str{1}).ExpRetAge;

            xi        = data.(member_str{1}).X;%[data.h.aux.X , data.w.aux.X];  % husband and wife information included. No constant because that is in the cut-offs
%             xi        = -xi;                    % this is how STATA does it...
            num_beta  = size(xi,2);       
            beta      = theta(1:num_beta);
            
            age_i     = ExpRetAge - par.ret_min + 1; 

            % handle NaNs
            nan_y        = isnan(age_i);
            age_i(nan_y) = 1; % going to ignore these observations below so the used index does not matter.

            
            % Now there is a functional form for the age-dummies
            par         = fun.update_par(theta(num_beta+1:end),EstAgeDum,par); % update age dum parameters
            age_grid    = (min(age_i):max(age_i))' + par.ret_min-1;
            sex         = 1;%1*strcmp(member_str,'h') + 2*strcmp(member_str,'w'); 
            cut_offs    = par.RetAgeDum_fun(age_grid,sex,par);
            cut_offs    = [-LARGE_NUMBER;cut_offs;LARGE_NUMBER]; % add the end-points
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            age_i     = age_i + 1; % add one because there is now added an element in the bottom of cut_offs.
            age_dum   = cut_offs(age_i);

            % lagged age-dummies
            age_i_l   = age_i - 1;
            age_dum_l = cut_offs(age_i_l);
            
            % calculate objective
            age_dum(nan_y)   = NaN; % construct NaN again
            age_dum_l(nan_y) = NaN;
            index   = xi*beta + age_dum; 
            index_l = xi*beta + age_dum_l; 
            
            Lik    = normcdf(index) - normcdf(index_l);
            obj    = -nanmean(log(Lik));
            
        end
        function b_oprobit  = oprobit_est(data,par,EstAgeDum)
            
            % reduce the dimension to one.
            par.RetAge_cons  = 0;
            par.RetAge_trend = 0.01;
            par.RetAge50     = 0;
            par.RetAge55     = 0;
            par.RetAge60     = 0;
            par.RetAge65     = 0;
            par.RetAge70     = 0;
            
            % 1. Estimate ordered probit models for each household member
            options_unc = optimoptions('fminunc','Display','off','Algorithm','quasi-newton','SpecifyObjectiveGradient',false,'CheckGradients',false,'MaxFunEvals',1000*6,'MaxIter',1000*6);
            members = {'h','w'};
            for m=1:numel(members)
                num_theta   = size(data.(members{m}).X,2);
                theta_beta  = 0*ones(num_theta,1);
                % change this to functional form
                num_age     = numel(EstAgeDum);
                theta_age   = [zeros(1,num_age)+.01]';

                theta0      = [theta_beta;theta_age];
                
                b_oprobit.(members{m}) = fminunc(   @(theta) fun.oprobit_obj(theta,data,par,members(m),EstAgeDum) , theta0, options_unc);

            end
        end
    end
end