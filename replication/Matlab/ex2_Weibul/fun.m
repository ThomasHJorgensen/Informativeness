classdef fun
    methods (Static)
        function [mm,S]=MoM(theta,T,x1,x2,x3)
            
            [n,k]=size(x1);
            m=zeros(n,k+2);
            theta1=theta*theta(length(theta));
            theta1(length(theta))=theta(length(theta));
            [ZZ]=fun.Z(T,x1,x2,x3,ones(n,1),theta1);
            a=log(ZZ)+double(eulergamma);
            
            ii=1;
            m(:,ii)=a;
            ii=ii+1;
            m(:,ii)=a.*x1(:,2);
            ii=ii+1;
            for i1=3:k
                m(:,ii)=a.*x1(:,i1);
                ii=ii+1;
                m(:,ii)=a.*x2(:,i1);
                ii=ii+1;
                m(:,ii)=a.*x3(:,i1);
                ii=ii+1;
            end
            mm=mean(m,1)';
            if nargout>1
                S=cov(m);
            end
            
        end
        
        function [Zout]=Z(t,x1,x2,x3,eta,theta)
            
            k=length(theta)-1;
            beta=theta(1:k);
            alpha=theta(k+1);
            
            xb1=exp(x1*beta);
            xb2=exp(x2*beta);
            xb3=exp(x3*beta);
            ta=t.^alpha;
            
            i1=(t <= 1.0);
            i2=(t <= 2.0)-i1;
            i3=1-i1-i2;
            
            Zout=i1.*ta.*xb1+i2.*((ta-1).*xb2+xb1)+i3.*((ta-2^alpha).*xb3+(2^alpha-1).*xb2+xb1);
            Zout=Zout.*eta;

        end

        function [SE,grad,sens,Avar] = sensitivity(theta,T,x1,x2,x3,S,W)
            
            % Calculate numerical gradient for all parameters
            num_mom  = size(W,1);
            num_par  = length(theta);
            h    = 1.0e-4;
            grad = NaN(num_mom,num_par);
            for i=1:numel(theta)
                var_now      = zeros(size(theta));
                var_now(i)   = 1;
                
                [forward]  = fun.MoM(theta+h*var_now,T,x1,x2,x3);
                [backward] = fun.MoM(theta-h*var_now,T,x1,x2,x3);
                
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
            [val0]      = fun.MoM(theta,T,x1,x2,x3);
            
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
        
        function write_res(meas,meas_str,lbl)
            num_par = numel(lbl);
            num_mom = size(meas{1},2);
            cut_off = 100;
            
            for m=1:numel(meas_str)
                meas_now = meas{m};
                too_large = meas_now>cut_off;
                star_col = max(too_large);
                
                fprintf('& \\multicolumn{%d}{c}{$%s$} \\\\ \\cmidrule(lr){2-%d}  \n',num_mom,meas_str{m},num_mom+1);
                for p=1:num_par
                    fprintf(' %s',lbl{p});
                    for j=1:num_mom
                        if too_large(p,j)
                            fprintf('& $> \\negthinspace 100^*$');
                        elseif star_col(j)
                            fprintf('& $%2.3f^*$',meas_now(p,j));
                        else
                            fprintf('& $%2.3f$',meas_now(p,j));
                        end
                    end
                    fprintf('\\\\ \n');
                end
                if m<numel(meas_str)
                    fprintf('\\cmidrule(lr){2-%d} \n',num_mom+1);
                end
            end
        end
        
        function print(x,s)
            
            if nargin<2
                a=floor(log10(max(max(abs(x)))+0.1))+1;
                if a>0
                    b=a+3;
                    s=[num2str(b+3) '.3'];
                else
                    b=2-a;
                    s=[num2str(b+6) '.' num2str(b+3)];
                end
            end
            
            [k1,k2]=size(x);
            
            formatSpec = ['%' s 'f '];
            for i1=1:k1
                for i2=1:k2
                    fprintf(formatSpec,x(i1,i2))
                end
                fprintf('\n')
            end
        end
        
    end
end


