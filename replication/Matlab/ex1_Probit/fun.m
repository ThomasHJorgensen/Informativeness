classdef fun
    methods (Static)
        
        function [mm,S]=MoM(b,y,x)
            
            a=y-normcdf(x*b);
            [n,k]=size(x);
            m=zeros(n,k*(k+1)/2);
            
            ii=1;
            for i1=1:k
                for i2=i1:k
                    m(:,ii)=a.*x(:,i1).*x(:,i2);
                    ii=ii+1;
                end
            end
            mm=mean(m,1)';
            if nargout>1
                S=cov(m);
            end
        end
        
        function [SE,grad,sens,Avar] = sensitivity(theta,y,x,S,W)
            
            % Calculate numerical gradient for all parameters
            num_mom  = size(W,1);
            num_par  = length(theta);
            h    = 1.0e-4;
            grad = NaN(num_mom,num_par);
            for i=1:numel(theta)
                var_now      = zeros(size(theta));
                var_now(i)   = 1;
                
                [forward]  = fun.MoM(theta+h*var_now,y,x);
                [backward] = fun.MoM(theta-h*var_now,y,x);
                
                grad(:,i)   = (forward-backward)./(2*h);
            end
            
            
            % calculate objects re-used below
            GW       = grad'*W;
            GWG      = GW*grad;
            
            
            % Calculate asymptotic variance and (adjusted) standard errors
            Avar    = (GWG\(GW*S*GW'))/GWG;
            AvarOpt = inv((grad'/S)*grad);
            
            %Avar_opt = (grad'/S)*grad;
            SE      = sqrt( diag(Avar) );
            
            % 4. calculate sensitivity measure (NumPar x NumMom) matrix
            sens.M1       = -GWG\GW;
            [val0]      = fun.MoM(theta,y,x);
            
            sens.M1e      = sens.M1.*( repmat(val0' , num_par ,1)./repmat(theta , 1, num_mom ) );
            
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
                M3kk     = (GWG)\(GW*O*GW')/(GWG);
                M6kk     = - GWG\(grad'*O*grad)*Avar ...
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
                sens.M4e(:,k) = (diag(Avar_now) - diag(Avar))./ diag(Avar);
                
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
        
        function write_res(res,titl,lbl,hdr)
            
            [k1, k2]=size(res);
            
            fprintf('\\begin{center}\n');
            fprintf('\\begin{tabular}{l');
            for i2=1:k2
                fprintf('r');
            end
            fprintf('}\n');
            k3=k2+1;
            fprintf('\\multicolumn{%s}',num2str(k3))
            fprintf('{c}{\\textbf{ %s}} \\\\ \n',titl)
            
            for i2=1:k2
                fprintf('& %s   ',char(hdr(i2)));
            end
            fprintf('\\\\ \n ');
            
            for i1=1:k1
                fprintf('%s',char(lbl(i1)));
                for i2=1:k2
                    fprintf('& $ %8.3f $ ',res(i1,i2));
                end
                fprintf('\\\\ \n ');
            end
            fprintf('\\end{tabular}\n');
            fprintf('\\end{center}\n\n\n\n');
        end
    end
end


