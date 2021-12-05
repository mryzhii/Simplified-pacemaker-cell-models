%%%% CORRADO Oscillator model  BACKWARD EULER
%%%% 
%%%% COMPARISON 60s FREQUENCY of Backward Euler and Forward Euler methods with different dt
% 
clear all
time_units = 30000; %60000;
tau_in = 0.3;
tau_out = 6.0;
tau_open = 120.0;
tau_close = 150.0;
u_s =  0.15; 
u_gate = -0.05; 
bCN = 0.2;

dt_backward = 0.0001;  % time step for backward E [ms]
dt_forward = 0.01;     % time step for forward E [ms]
delta = 0.1; % 1.0

nruns = 1;
time_BE = zeros(nruns,1);
time_FE = zeros(nruns,1);

for run = 1:nruns
 fprintf("Run %d\n",run);

for nt = 1:2
    u  = 0.01;
    u1 = 0.01; 
    h  = 0.5;
    h1 = 0.5;
    switch nt
    case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Backward Euler
        tic
        T = time_units/dt_backward;
        si1 =delta/dt_backward;
        ts_U_BE = zeros(1,T/si1);
        ts_T_BE = zeros(1,T/si1);
        ts_UH_BE = zeros(2,floor(numel(T)/si1));
        for t = 2:T 
          iter   =  0;
          Residual    =  1.0e6; 
% Reaction terms for the fast and slow variables 
          while abs(Residual)>1.0e-7 && iter <= 20
            h_inf = 0.5*(1.0-tanh((u-u_gate)/u_s)   );
            tau = tau_open*tau_close/(tau_open+h_inf*(tau_close-tau_open));
%  Backward Euler
            h1 = (h*tau+h_inf*dt_backward)/(tau+dt_backward);
            I = h1*u1*(u1+bCN)*(1.0-u1)/tau_in - (1.0-h1)*u1/tau_out;
            dIt = h1/tau_in*( (u1+bCN)*(1.0-u1)+ u1*(1.0-u1) - u1*(u1+bCN) ) - (1.0-h1)/tau_out;
            Residual = (u1-u)/dt_backward - I;
            u1 = u1 - Residual/(1.0/dt_backward + dIt);
            iter = iter + 1;
          end
          if iter > 19, fprintf("Warning: it = 20\n"); end
          u = u1;
          h = h1;     
%  Downsample to create output matrix         
          if rem(t,si1) == 0
            l = floor(t/si1);
            ts_U_BE(1,l) = u;
            ts_UH_BE(1,l) = u;
            ts_UH_BE(2,l) = h;
            ts_T_BE(1,l) = t*dt_backward;
          end
       end % t
       time_BE(run) = toc;
       [peaks,locs,widths,proms]=findpeaks(ts_U_BE(1,1:end),ts_T_BE(1:end),...
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
        Ampl_BE = max(proms(1:end));
        Maxp_BE = max(peaks(1:end));
       end
       fprintf('BE: dt=%0.5f Period_BE=%0.4f Freq_BE=%0.4f Ampl_BE=%0.2f \n',...
                    dt_backward,Period_BE,Freq_BE,Ampl_BE);

    case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Forward Euler
        tic
        T = time_units/dt_forward;
        si2 = delta/dt_forward;
        ts_U_FE = zeros(1,T/si2);
        ts_T_FE = zeros(1,T/si2);
        ts_UH_FE = zeros(2,T/si2);
        
       for t = 1:T 
          h_inf = 0.5*(1.0-tanh((u-u_gate)/u_s)   );
          tau = tau_open*tau_close/(tau_open+h_inf*(tau_close-tau_open));
          del_u = h*u*(u + bCN)*(1.0-u)/tau_in - (1.0-h)*u/tau_out; 
          del_h = (h_inf-h)/tau; 
%  Forward Euler
          u1 = u + dt_forward*del_u;
          h1 = h + dt_forward*del_h; 
          u = u1;
          h = h1;
%  Downsample to create output matrix          
          if rem(t,si2) == 0
            l = floor(t/si2);
            ts_U_FE(1,l) = u;
            ts_UH_FE(1,l) = u;
            ts_UH_FE(2,l) = h;
            ts_T_FE(1,l) = t*dt_forward;
          end %
       end % t
       time_FE(run) = toc;
       
       [peaks2,locs2,widths2,proms2]=findpeaks(ts_U_FE(1,1:end),ts_T_FE(1:end),...
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
         Ampl_FE = max(proms2(1:end));
         Maxp_FE = max(peaks2(1:end));
       end
       fprintf('FE: dt=%0.5f Period_FE=%0.4f Freq_FE=%0.4f Ampl_FE=%0.2f \n',...
           dt_forward,Period_FE,Freq_FE,Ampl_FE);
     
    end % select
  end % nt
end % runs

%%% Vector Norm 
nn1 = locs2(1)/delta-locs(1)/delta;  
%nstop = time_units/delta-nn1; 
%nstart = nstop -3500; 
nstop = 3500;
nstart = 1;
Udiff_matL_2_rel   = norm(ts_U_BE(1,nstart:nstop)-ts_U_FE(1,nstart+nn1:nstop+nn1),2)/norm(ts_U_BE(1,nstart:nstop),2);
Udiff_matL_inf_rel = norm(ts_U_BE(1,nstart:nstop)-ts_U_FE(1,nstart+nn1:nstop+nn1),Inf)/norm(ts_U_BE(1,nstart:nstop),Inf);

%%% Comparison
fprintf('CN BackwardEuler / ForwardEuler:\n     dt1=%0.4f dt2=%0.4f | delta=%0.2f\n',...
    dt_backward,dt_forward,delta); 
fprintf(' Udiff_L2_rel  = %0.8f %0.4f%% \n',  Udiff_matL_2_rel,   Udiff_matL_2_rel*100);
fprintf(' Udiff_Linf_rel= %0.8f %0.4f%% \n',  Udiff_matL_inf_rel, Udiff_matL_inf_rel*100);

fprintf(' Freq_bE-Freq_fE=%0.8e  d_Freq_rel=%0.8e %0.8f%% \n',...
    Freq_BE-Freq_FE,abs(Freq_BE-Freq_FE)/Freq_BE,abs(Freq_BE-Freq_FE)/Freq_BE*100);
fprintf(' time1/time2=%0.2f \n',sum(time_BE./time_FE)/nruns );

fname = 'CN_BE-FE_test.txt';
fileID = fopen(fname,'w');

txt = strcat('# t_in=',num2str(tau_in,'%0.2f'),' t_out=',num2str(tau_out,'%0.2f'),' t_open=',num2str(tau_open,'%0.2f'),...
        ' t_close=',num2str(tau_close,'%0.2f'),' u_s=',num2str(u_s,'%0.3f'),' u_gate=',num2str(u_gate,'%0.3f'),' bCN=',num2str(bCN,'%0.3f'),'\n');
fprintf(fileID,txt);
fprintf(fileID,'CN BackwardEuler / ForwardEuler:\n     dt1=%0.4f dt2=%0.4f | delta=%0.2f\n',dt_backward,dt_forward,delta); 
fprintf(fileID,' Udiff_L2_rel  = %0.8f %0.4f%% \n',  Udiff_matL_2_rel,   Udiff_matL_2_rel*100);
fprintf(fileID,' Udiff_Linf_rel= %0.8f %0.4f%% \n',  Udiff_matL_inf_rel, Udiff_matL_inf_rel*100);

fprintf(fileID,' Freq_bE-Freq_fE=%0.8e  d_Freq_rel=%0.8e %0.8f%% \n',...
    Freq_BE-Freq_FE,abs(Freq_BE-Freq_FE)/Freq_BE,abs(Freq_BE-Freq_FE)/Freq_BE*100);


fprintf(fileID,' time_BE/timeFE = %0.2f\n',sum(time_BE./time_FE)/nruns );
fclose(fileID);


%%%%%%%%%%%%%%%%%%%%%%%%%% Action Potentials
figure(1);
clf
subplot(1,2,1)
set(gcf,'Position',[100 450 600 250]);
%box on
hold on
plot(ts_T_FE(floor(end*9/10):end).*1.e-3,ts_U_FE(floor(end*9/10):end),'Color',[0.7 0.7 0.7],'LineWidth',2.0)
plot(ts_T_BE(end*9/10:end).*1.e-3,ts_U_BE(end*9/10:end),'--r','LineWidth',2.0)
%legend('FE','BE');
xlabel('Time (s)','fontsize',14);
ylabel('u','fontsize',14);
set(gca,'FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%% Phase portrait
subplot(1,2,2)
%box on
hold on
plot(ts_UH_FE(1,floor(end*3/4):end),ts_UH_FE(2,floor(end*3/4):end),'Color',[0.7 0.7 0.7],'LineWidth',4.0)
plot(ts_UH_BE(1,floor(end*3/4):end),ts_UH_BE(2,floor(end*3/4):end),'--r','LineWidth',1.0)
legend('FE','BE');
xlabel('u','fontsize',14);
ylabel('h','fontsize',14);
set(gca,'FontSize',12);