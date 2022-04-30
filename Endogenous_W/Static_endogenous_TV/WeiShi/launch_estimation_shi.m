function [outputs,allresults]= launch_estimation_shi(depend,exo,endo_var,...
                                   endo_varx,weight,info)
% endo_var = depend_z;
% endo_varx = {};
% sdata = estim_data;
% weight = weights;
     
    fieldsinf = fieldnames(info);
    nf = length(fieldsinf);
    for i=1:nf
        if strcmp(fieldsinf{i},'ind_id')
            ind_id=info.ind_id;
            % Individual Id's
            id = sdata.(sprintf('%s',ind_id));
        elseif strcmp(fieldsinf{i},'time_id')
            time_id=info.time_id;
            % Time variable (think about the case where the time period is not
             % properly ordered)
            period = sdata.(sprintf('%s',time_id));
        elseif strcmp(fieldsinf{i},'print_res')
            print_res = info.print_res;
        elseif strcmp(fieldsinf{i},'n')
            n=info.n;
        elseif strcmp(fieldsinf{i},'T')
            T= info.T; 
        elseif strcmp(fieldsinf{i},'optimtyp')
            optimtyp = info.optimtyp; 
        elseif strcmp(fieldsinf{i},'boundtyp')
            boundtyp = info.boundtyp; 
        elseif strcmp(fieldsinf{i},'modeltyp')
            modeltyp = info.modeltyp;
        elseif strcmp(fieldsinf{i},'lagexo')
            lagexo = info.lagexo;
        elseif strcmp(fieldsinf{i},'varset')
            varset = info.varset;
        end
    end

    % dependent variable 
    y = depend;
    yori = y;
    y = reshape(y,n,T);
%     if strcmp(varset,'VarAOut')==1
%          y = reshape(y,n,T);
%     else
%         y = reshape(y,n,T+1);
%         y = y(:,2:end);
%         yori = reshape(y,n*T,1);
%     end
    
    % independent variables
    ky = size(exo,2);
    xtemp = [];
    xtempori = [];
    for li = 1:ky
       
            xini = exo(:,li);
    
        xtempori = [xtempori,xini];
        xini = reshape(xini,n,T);
        xtemp = [xtemp, xini];
%         if strcmp(varset,'VarAOut')==1
%            xini = reshape(xini,n,T);
%            xtemp = [xtemp, xini];
%         else
%            xini = reshape(xini,n,T+1);
%            xtemp = [xtemp, xini(:,1:end-1)];
%            xtempori = [xtempori,reshape(xini(:,1:end-1),n*T,1)];
%         end        
    end 
    exo_temp = exo;
    x = xtemp;
    xori = xtempori;
    %Xy = xesfact;
    
    
        kz = size(endo_varx,2);
        x_ztemp = [];
        xztempori = [];
        for li = 1:kz
            xzini   = endo_varx(:,li);
            xztempori = [xztempori,xzini];
            xzini   =  reshape(xzini,n,T);
            x_ztemp = [x_ztemp,xzini];
            %x_ztemp = [x_ztemp, xzini(:,1:end-1)];
        end 
        x_z = x_ztemp;
        xzori = xztempori;
     
    % Compute correlation between variables
    % correlation
    %[corr_xtemp, ~] = corr(xtemp);
    %[corr_x, ~] = corr(x);
    
     n_beta = exo;
          
  
 
   
    %initial_values = zeros(ky+2,1);

   % Estimation for each weight matrices
   
    
        %iw = 10;      
        % Specific weight matrice
        
        W = weight;
       
        
        % Variables to control for the endogeneity of interaction matrix
        z = endo_var;
        zori = z;
        %if strcmp(varset,'VarAOut')==1
            z = reshape(z,n,T);
        %else
         %   z = reshape(z,n,T+1);
          %  z = z(:,2:end);
           % zori = reshape(z,n*T,1);
        %end
        
         if strcmp(optimtyp,'Shi_con')==1
           if boundtyp == 2 
              % kz=1;
              initial_values = [zeros(kz+ky+4,1),...
                  [randn(kz+ky,9);-1 + 2.*rand(1,9); 10.*rand(2,9);randn(1,9)]];
              initial_values(end-2:end-1) = exp(initial_values(end-2:end-1));
           elseif boundtyp == 1
              initial_values = [zeros(kz+ky+4,1),...
                      [randn(kz+ky,9);-1 + 2.*rand(1,9); 10.*rand(2,9);randn(1,9)]];
           end
       elseif strcmp(optimtyp,'Shi_unc')==1
           olsy = ols(yori,xori);
           olsz = ols(zori,xzori);
          initial_values = [zeros(kz+ky+4,1),...
           [olsz.beta;olsy.beta;randn(1,1);olsy.sige;olsz.sige;randn(1,1)],...
           randn(kz+ky+4,8)];
       end
%         foldnam = strcat('C:\Users\Sessi\Dropbox\Phd_thesis\shock_propagation_financial_inst\Output\Estimation\',...
%         modeltyp,'\Z_variables_as_ratio_to_tot_asset');  
%         cd(foldnam)
        
%         filnam = strcat('sensitivity_',lower(strrep(optimtyp,'_','')),'_',...
%                      modeltyp,'_',varset,'.xlsx');
%         xlswrite(filnam,initial_values,W_nam,'B3')

        [estim, selectstat,goodfit,nfactors,exitflag,results] = ...
              estimodel_shi(y,x,z,x_z,W,initial_values,info,exo,endo_varx,W);
          
%         cd(foldnam)
%         xlswrite(filnam,n_param,W_nam,'A2')     
        
        
            allestimate=estim(:,1);
            allstds=estim(:,2);
            alltstat=estim(:,3);
            allselect=selectstat; 
            allnfactors = nfactors;
            allexitflag = exitflag;
            allini_val = initial_values;
            allgoodfit = goodfit;
       
        
        allresults.(sprintf('%s','W')) = results;

    
    outputs = struct('allestimate',allestimate,'allstds',allstds,...
        'alltstat',alltstat,'allselect',allselect,'allnfactors',allnfactors,...
        'allexitflag',allexitflag);

%     if print_res == 1
%         %%%%%% Display results %%%%%%%
%          fieldm = field'
%          if length(char(fieldm))<=15
%               kfac=13;
%           else
%               kfac=17;
%          end
% 
%         % Correlation between variables
% %         display('Correlation between variables')
% %         display('------------------------------')
% %         display('Case 1: variables with no lags')
% % 
% %         fprintf(strcat('\n                     \t                ',repmat('%13s\t',1,size(exo_temp,1)), '\n'),exo_temp{:});
% %         for i=1:size(exo_temp,1)
% %             fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(exo_temp,1)), '\n'),exo_temp{i},corr_xtemp(i,:));
% %         end
% %         display('Case 2: variables with lags')
% %         fprintf(strcat('\n                     \t                ',repmat('%13s\t',1,size(exo,1)), '\n'),exo{:});
% %         for i=1:size(n_beta,1)
% %             fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(n_beta,1)), '\n\n'),n_beta{i},corr_x(i,:));
% %         end
% 
%         fprintf('\n')
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         if size(field)<4
%             fprintf('\n                              Weight matrices                                   \n')
%         else
%             fprintf('\n                                                       Weight matrices                                                      \n')
%         end
% 
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         fprintf(strcat('\n          Parameters | \t                ',repmat('%13s\t',1,size(fieldm,2)), '\n'),fieldm{:});
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         fprintf('\n%20s\n','Parameters in the main equation y');
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         
%         for i=kz+1:kz+ky
%                 fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),n_param{i},allestimate(i,:));
%                 fprintf(strcat('                    \t                ',repmat('      (%4.3f)  \t',1,size(fieldm,2)), '\n'),allstds(i,:));
%                 fprintf(strcat('                    \t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),alltstat(i,:));
%         end
%         fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),n_param{kz+ky+1},allestimate(kz+ky+1,:));
%         fprintf(strcat('                    \t                ',repmat('      (%4.3f)  \t',1,size(fieldm,2)), '\n'),allstds(kz+ky+1,:));
%         fprintf(strcat('                    \t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),alltstat(kz+ky+1,:));
%         fprintf(strcat('\n%20s\t                ',repmat('%12.6f\t',1,size(fieldm,2)), '\n'),n_param{end},allestimate(end,:));
%         fprintf(strcat('                    \t                ',repmat('      (%4.6f)  \t',1,size(fieldm,2)), '\n'),allstds(end,:));
%         fprintf(strcat('                    \t                ',repmat('%12.6f\t',1,size(fieldm,2)), '\n'),alltstat(end,:));
%         
% 
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         fprintf('\n%20s\n','Parameters in the endogeneity control equation z');
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         for i= 1:kz
%                 fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),n_param{i},allestimate(i,:));
%                 fprintf(strcat('                    \t                ',repmat('      (%4.3f)  \t',1,size(fieldm,2)), '\n'),allstds(i,:));
%                 fprintf(strcat('                    \t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),alltstat(i,:));
%         end
% 
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         fprintf('\n%20s\n','Components of errors variance');
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
% 
%         for i= kz+ky+2:size(n_param,1)-1
%                 fprintf(strcat('\n%20s\t                ',repmat('%12.6f\t',1,size(fieldm,2)), '\n'),n_param{i},allestimate(i,:));
%                 fprintf(strcat('                    \t                ',repmat('      (%4.6f)  \t',1,size(fieldm,2)), '\n'),allstds(i,:));
%                 fprintf(strcat('                    \t                ',repmat('%12.6f\t',1,size(fieldm,2)), '\n'),alltstat(i,:));
%         end
%         
% 
% %         fprintf(repmat('-',1,kfac*length(char(fieldm))))
% %         fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),n_param{end},allestimate(end,:));
% %         fprintf(strcat('                    \t                ',repmat('      (%4.3f)  \t',1,size(fieldm,2)), '\n'),allstds(end,:)); 
% %         fprintf(strcat('                    \t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),alltstat(end,:)); 
%         
%  fit_nam = {'Rsquare1'};
%         for i=1:size(fit_nam,1)
%             fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),fit_nam{i},allgoodfit(i,:));
%         end
%         
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         select_nam={'Likhood';'YLikhood';'ZLikhood'; 'AIC'; 'BIC'};
% 
%         for i=1:size(select_nam,1)
%             fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),select_nam{i},allselect(i,:));
%         end
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         fprintf('\n%20s\n','Additional information');
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         inform ={'Nfactors Y';'Nfactors Z';'Exitflag'};
% 
%         for i=1:size(inform ,1)-1
%             fprintf(strcat('\n%20s\t                ',repmat('%3.0f\t',1,size(fieldm,2)), '\n'),inform{i},allnfactors(i,:));
%         end
%         fprintf(strcat('\n%20s\t                ',repmat('%3.0f\t',1,size(fieldm,2)), '\n'),inform{end},allexitflag(1,:));
%         
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         fprintf('\n%20s\n','Initial conditions');
%         fprintf(repmat('-',1,kfac*length(char(fieldm))))
%         for i=kz+1:kz+ky
%                 fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),n_param{i},allini_val(i,:));
%         end
%         fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),n_param{kz+ky+1},allini_val(kz+ky+1,:));
%         fprintf(strcat('\n%20s\t                ',repmat('%12.6f\t',1,size(fieldm,2)), '\n'),n_param{end},allini_val(end,:));
%         for i= 1:kz
%                 fprintf(strcat('\n%20s\t                ',repmat('%12.3f\t',1,size(fieldm,2)), '\n'),n_param{i},allini_val(i,:));
%         end
%         for i= kz+ky+2:size(n_param,1)-1
%                 fprintf(strcat('\n%20s\t                ',repmat('%12.6f\t',1,size(fieldm,2)), '\n'),n_param{i},allini_val(i,:));
%         end
%           
%        
% 
%        fprintf('\n-------------------------------------------------------------------------------------------------------------------------------------------------------\n \n');
% 
%    end
    
end