%   ********************************************************************
%   * Copyright(c) M. Ryzhii, University of Aizu, Japan                *
%   *              E. Ryzhii, Fukushima Medical University, Japan      *
%   *              05/12/2021                                          *
%   * (2022) "Pacemaking function of two simplified cell models",      *
%   *  PLoS ONE 17(4): e0257935. doi.org/10.1371/journal.pone.0257935  *
%   ********************************************************************
% Pacemaker (oscillator) variant of the Aliev & Panfilov model (1996)
% Comparison of Backward Euler and Forward Euler methods with different time step dt
%
clear all
total_time = 60000; % Time in [ms]
% Constants of the model
k    = 8.0;
a    = 0.15; 
mu1  = 0.20;
mu2  = 0.30;
ct   = 1.0/12.9; % Time scaling coefficient
eps0 = 0.002;
bAP  = 0.02;

dt_backward = 0.0001;   %0.01;  % Time step for backward Euler [ms]
dt_forward  = 0.01;     % Time step for forward Euler [ms]
delta  = 0.1;
si1 = delta/dt_backward; % Output intervals
si2 = delta/dt_forward;
        
nruns = 1; % Set nruns >1 for precise calculation of the simulation time
sim_time_BE = zeros(nruns,1);
sim_time_FE = zeros(nruns,1);

spat = 0.0; % 0.001 % Spatial coupling term (external current)

fprintf('Pacemaking Aliev-Panfilov model: BackwardEuler / ForwardEuler:\n');  
for run = 1:nruns
    if nruns>1, fprintf("Run %d\n",run); end
    
for nt = 1:2
    u  =  0.01;
    u1 =  0.01;
    v  =  0.01;
    v1 =  0.01;
    switch nt
    case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Backward Euler  
      tic
      T = total_time/dt_backward;
      ts_uvt_BE = zeros(3,T/si1); % Matrix with u, v, and t
      for t = 1:T 
        iter = 0;
        residual0 = 1.0;
        while residual0 > 1.0e-7 && iter <= 20
            dtt = dt_backward*ct;
            tmp  =  (1.0 + dtt*(eps0 + mu1*k*u1*(u1-a-1.0)/(u1+mu2)) );
            v1   =  (-tmp + sqrt(tmp^2 - 4.0*dtt*mu1/(u1+mu2)*(dtt*eps0*k*u1*(u1-a-1.0)- v)))/(2.0*dtt*mu1/(u1+mu2));
            dvdu  =  mu1*v1/(u1+mu2)*(v1 + k*u1*(u1-a-1.0)) - k*(1.0*u1-a-1.0)*(eps0*(u1+mu2)+mu1*v1)/(mu1*(v1+k*u1*(u1-a-1.0)) + (2.0/dtt+eps0)*(u1+mu2));
            itotal =  k*u1*(u1+bAP)*(u1-1.0) + u1*v1 - spat/ct;
            ditotal =  k*(u1*(3.0*u1-2.0) + bAP*(2.0*u1-1.0)) + v1 + u1*dvdu;
            residual =  (u1 - u)/dtt + itotal; 
            u1    =  u1 - residual/(1.0/dtt + ditotal);
            iter  =  iter + 1;
            %residual0 = max(abs(residual), [], 'all'); % In the case of tissue
            residual0 = abs(residual);
        end % while
        if iter > 20, fprintf("Warning: iter > 20, Residual = %e\n", residual0); end
        u = u1;
        v = v1;
%  Downsample to create output matrix
        if rem(t,si1) == 0
              j = floor(t/si1);
              ts_uvt_BE(1,j) = u;
              ts_uvt_BE(2,j) = v;
              ts_uvt_BE(3,j) = t*dt_backward;
        end
      end  % t
      sim_time_BE(run) = toc;
      [peaks,locs,widths,proms] = findpeaks(ts_uvt_BE(1,:),ts_uvt_BE(3,:),...
        'MinPeakHeight',0.1,'MinPeakDistance',0.10);
        Period_BE = 1.e-3*mean(diff(locs));   % In [s]
      if isnan(Period_BE) 
        nloc = 0;
        Freq_BE = NaN;
        Ampl_BE = 0;
        Maxp_BE = 0;  
      else
        nloc = length(locs);
        Freq_BE = 1.0/Period_BE;
        Ampl_BE = max(proms(floor(end/2):end));
        Maxp_BE = max(peaks(floor(end/2):end));
      end
      fprintf('BE: dt_backward = %0.5f  Period_BE = %0.5f  Freq_BE = %0.4f  Ampl_BE = %0.4f \n',...
                    dt_backward,Period_BE,Freq_BE,Ampl_BE);

    case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Forward Euler
      tic  
      T = total_time/dt_forward;
      ts_uvt_FE = zeros(3,T/si2); % Matrix with u, v, and t
      for t = 1:T
         dudt = ct*(k*u*(u + bAP)*(1.0-u) - u*v) + spat;
         dvdt = ct*(eps0 + mu1*v/(u+mu2))*(-v-k*u*(u-a-1.0));
         u =  u + dt_forward*dudt; 
         v =  v + dt_forward*dvdt; 
 %  Downsample to create output matrix         
         if rem(t,si2) == 0
            j = floor(t/si2);
            ts_uvt_FE(1,j) = u;
            ts_uvt_FE(2,j) = v;
            ts_uvt_FE(3,j) = t*dt_forward;
         end      
      end % t
      sim_time_FE(run) = toc;  
      [peaks2,locs2,widths2,proms2] = findpeaks(ts_uvt_FE(1,:),ts_uvt_FE(3,:),...
      'MinPeakHeight',0.1,'MinPeakDistance',0.10);
      Period_FE = 1.e-3*mean(diff(locs2));   % In [s]
      if isnan(Period_FE) 
         nloc2 = 0;
         Freq_FE = NaN;
         Ampl_FE = 0;
         Maxp_FE = 0;  
      else
         nloc2 = length(locs2);
         Freq_FE = 1.0/Period_FE;
         Ampl_FE = max(proms2(floor(end/2):end));
         Maxp_FE = max(peaks2(floor(end/2):end));
      end
 
      fprintf('FE: dt_forward  = %0.5f  Period_FE = %0.5f  Freq_FE = %0.4f  Ampl_FE = %0.4f\n',...
            dt_forward,Period_FE,Freq_FE,Ampl_FE);
    end % switch
  end % nt
end % runs

if (nloc == 0 || nloc2 == 0)
    fprintf(' - No oscillations -\n'); 
else
  nn1 = floor(abs(locs2(1)/delta-locs(1)/delta)); % Correction shift for BE-FE synchronization
  nend = 3500;
  nbegin = 1;

% Relative matrix norm^2 
  udiff_matL_2_rel   = norm(ts_uvt_BE(1,nbegin:nend)-ts_uvt_FE(1,nbegin+nn1:nend+nn1),2)/norm(ts_uvt_BE(1,nbegin:nend),2);
% Relative matrix norm^Inf 
  udiff_matL_inf_rel = norm(ts_uvt_BE(1,nbegin:nend)-ts_uvt_FE(1,nbegin+nn1:nend+nn1),Inf)/norm(ts_uvt_BE(1,nbegin:nend),Inf);

% Comparison
  fprintf(' L_2 relative norm   = %0.8f / %0.4f%% \n',  udiff_matL_2_rel,   udiff_matL_2_rel*100);
  fprintf(' L_inf relative norm = %0.8f / %0.4f%% \n',  udiff_matL_inf_rel, udiff_matL_inf_rel*100);
  fprintf(' Frequency_BE-Frequency_FE = %0.5e  d_Frequency_rel = %0.5e / %0.3f%%\n',...
    Freq_BE-Freq_FE,abs(Freq_BE-Freq_FE)/Freq_BE,abs(Freq_BE-Freq_FE)/Freq_BE*100);
  fprintf(' sim_time_BE/sim_time_FE = %0.2f \n',sum(sim_time_BE/sim_time_FE)/nruns );

%%%%%%%%%%%%%%%%%%%%%%%%%% Plot of action potentials
  Fig = figure();
  clf
  set(gcf,'Position',[100 450 700 250]);
  subplot(1,2,1)
  title('pAP: Action potentials'); 
  box on
  hold on; grid on
  plot(ts_uvt_FE(3,end-floor(end/20):end).*1.e-3,ts_uvt_FE(1,end-floor(end/20):end),'Color',[0.2 0.2 0.7],'LineWidth',2.0)
  plot(ts_uvt_BE(3,end-floor(end/20):end).*1.e-3,ts_uvt_BE(1,end-floor(end/20):end),'-r','LineWidth',1.0)
  xlabel('Time (s)','FontSize',10);
  ylabel('u','FontSize',10);
  set(gca,'FontSize',10);
%%%%%%%%%%%%%%%%%%%%%%%%%% Plot of phase portraits
  subplot(1,2,2)
  title('pAP: Phase portraits'); 
  box on
  hold on
  plot(ts_uvt_FE(1,floor(end*7/8):end),ts_uvt_FE(2,floor(end*7/8):end),'Color',[0.2 0.2 0.7],'LineWidth',4.0)
  plot(ts_uvt_BE(1,floor(end*7/8):end),ts_uvt_BE(2,floor(end*7/8):end),'-r','LineWidth',1.0)
  str2 = sprintf('BE dt=%0.1e ms ',dt_backward); 
  str1 = sprintf('FE dt=%0.1e ms ',dt_forward); 
  legend(str1,str2,'Location',[0.40,0.016,0.231,0.17],'FontSize',10);
  xlabel('u','FontSize',10);
  ylabel('v','FontSize',10);
  set(gca,'FontSize',10);
  
  exportgraphics(Fig,'pAPmodel_BE_FE.png');
end
