use_spiking = 1;
use_theta = 1;
use_feedforward = 1;
use_envelope = 0;
use_periodic = 1;
use_facilitation = 1;

spatialDim = 1;

use_plot = 0;

% V = VideoWriter('thetaSequences_2.avi','Motion JPEG AVI'); V.Quality = 75;
% open(V);

for kk = 1:2

%parameters
    %simulation duration
    if spatialDim == 1
        T_transition = 4;%10;%2
        T_run = 6;%0;%8; %duration of run between stopping periods (sec)
        T_stop = 12; %duration of stopping periods
    else
        T_transition = 1;
        T_run = 6; %duration of run between stopping periods (sec)
        T_stop = 1; %duration of stopping periods
    end
    if use_plot == 1
        numStops = 4; %number of stops
    else
        if spatialDim == 1
            %numStops = 400; % Full model.
            numStops = 10;
        else
            numStops = 200;
        end
    end
    T = numStops*(T_run + T_stop + 2*T_transition); %total simulation time

    %network parameters
        if spatialDim==2, n = 64; %num of neurons per spatialDim
        elseif spatialDim==1, n = 64;
        end
        big = 2*n; %padding for convolutions
        N = n^spatialDim;
        dt = 0.5/1000; %step size
        if spatialDim==1
            m = 1; %CV = 1/sqrt(m)
        else
            m = 4;
        end
    
        if spatialDim==1
            tau_s = 30/1000; %synaptic time constant
            beta_0 = 400; %uniform excitation    
            w_rec = 60; %amplitude of recurrent exc.
            sig_rec = 0.1; %width of recurrent exc.
            w_inh = 0.07; %global inhibition
            w_adapt = 4; %amplitude of adaptation
            tau_adapt = 3; %adaptation time constant !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            w_FF = 300; %amplitude of feedforward exc.
            sig_FF = 0.01; %width of feedforward exc.
            w_facil = 0; %amplitude of facilitation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tau_facil = 16;
            v_max = 0.1;
            T_theta_max = 0.4; 
            T_replay_min = 0.8;

            if kk == 2, w_facil = 2; end
        else
            % tau_s = 20/1000; %synaptic time constant
            % beta_0 = 600; %uniform excitation    
            % w_rec = 32; %strength of recurrent exc.
            % sig_rec = 0.04; %width of recurrent exc.
            % w_inh = 0.004; %global inhibition
            % w_adapt = 40; %amplitude of adaptation %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % tau_adapt = 1; %adaptation time constant %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % w_FF = 600; %amplitude of feedforward exc.
            % sig_FF = 0.002; %width of feedforward exc.
            % v_max = 0.5;
            % T_theta_max = 0.12;

            tau_s = 20/1000; %synaptic time constant
            beta_0 = 600; %uniform excitation    
            w_rec = 24; %strength of recurrent exc.
            sig_rec = 0.02; %width of recurrent exc.
            w_inh = 0.002; %global inhibition
            w_adapt = 40; %amplitude of adaptation %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            tau_adapt = 0.8; %adaptation time constant %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            w_FF = 400; %amplitude of feedforward exc.
            sig_FF = 0.001; %width of feedforward exc.
            v_max = 0.35;
            T_theta_max = 0.12;

            if kk == 1, v_max = v_max_vec(kkk);
            elseif kk == 2, T_theta_max = T_theta_max_vec(kkk);
            elseif kk == 3, w_adapt = w_adapt_vec(kkk);
            else, tau_adapt = tau_adapt_vec(kkk);
            end
        end

        T_plot = T_theta_max/12;

        %place field centers
        if spatialDim == 2
            [X,Y] = meshgrid((1:n)/n,(1:n)/n); 
        else
            X = (1:n)'/n;
        end

        %activity envelope
        if use_envelope == 1
            kappa = 0.4; %controls width of main body of envelope
            a0 = 40;    %contrls steepness of envelope
            A = zeros(n,1);
            for i = 1:n
                r = abs(i-n/2);
                if r<kappa*n 
                    A(i) = 1;
                else
                    A(i) = exp(-a0*((r-kappa*n)/((1-kappa)*n))^2);
                end
            end
            if spatialDim == 2
                A = A*A';
            end
        else
            if spatialDim == 1
                A = ones(n,1);
            else
                A = ones(n,n);
            end
        end

%rat trajectory
    if spatialDim == 1
        %rat max speed
        % v_max = 0.1;

        v_transition = (1+sin(2*pi*(dt:dt:T_transition)'/2/(T_transition) - pi/2))/2;
        v_lap = v_max*[zeros(T_stop/dt,1); v_transition; ones(T_run/dt,1); v_transition(end:-1:1)];
        v = repmat(v_lap,1,numStops); v = v(:);
            if use_plot == 1, v = circshift(v,-(T_stop + T_transition + T_run/2)/dt); end %uncomment this if sim starts with rat running full speed
        x = mod(cumsum(v)*dt,1);
        acc = [abs(diff(v)); 0];

        % g_transition = 1-acc/max(acc);

        v_thr = 0.03;
        d = 0*(v>v_thr) + 1*(v<v_thr);
        [~,locs] = findpeaks(abs(diff(d)));
        g_transition = zeros(1,T/dt);
        g_transition(locs) = 1;
        g_transition = conv(g_transition,setUp_gaussFilt([1,20/dt],4000/dt),'same');
        g_transition = 1-g_transition/max(g_transition);
        % g_transition = ones(1,T/dt);

        % subplot(211)
        % plot(dt:dt:T,v,'k.','linewidth',2)
        % 
        % subplot(212)
        % plot(dt:dt:T,g_transition,'b','linewidth',2)
        % vline(dt*locs,'g--')
        % keyboard

        % v = v_max*ones(T/dt,1);
        % x = mod(cumsum(v)*dt,1);

        % v = zeros(T/dt,1);
        % x = 0.5*ones(T/dt,1);
    else
        %rat max speed
        % v_max = 0.4; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        % v_lap = v_max*(1+sin(2*pi*(dt:dt:T_run)'/(2*T_run) + pi/2))/2;
        % v_lap = [v_lap; zeros(T_stop/dt,1); v_lap(end:-1:1)];
        % v = repmat(v_lap,1,numStops); v = v(:);
        % x = [0.5*ones(T/dt,1) mod(0.2+cumsum(v)*dt,1)];

        v = v_max*ones(T/dt,1);
        % v = v_max*(1+sin(2*pi*(dt:dt:T)/40))/2;
        x = [0.5*ones(T/dt,1) mod(0.1+cumsum(v)*dt,1)];
        % x = [x(:,2) x(:,1)];

        %correlated random walk
            L = [2+8 98-8]; %size of box (cm)
            % L = [2 98]; %size of box (cm)
            dt_randomWalk = 0.1; %otherwise, takes too long to generate random walk
            speed_sigma = 0; %step size variance
            theta_sigma = 0.02*2*pi; %orientation variance (each step)

            v_transition = (1+sin(2*pi*(dt_randomWalk:dt_randomWalk:T_transition)'/2/(T_transition) - pi/2))/2;
            v_lap = v_max*[zeros(T_stop/dt_randomWalk,1); v_transition; ones(T_run/dt_randomWalk,1); v_transition(end:-1:1)];
            v = repmat(v_lap,1,numStops); v = v(:);
                % v_max_repmat = repmat(normrnd(v_max,0.1,1,numStops),length(v_lap),1);
                % v = repmat(v_lap,1,numStops).*v_max_repmat; 
                % v = v(:);
            speed_mean = 100*v;%30; %step size mean

            x = load_randomWalk(L,T,speed_mean,speed_sigma,theta_sigma,dt_randomWalk,[],spatialDim);
            x(:,2:3) = x(:,2:3)/100;
            v = [diff(x(:,2)), diff(x(:,3))]/mode(diff(x(:,1))); v = [v(1,:); v];
            speed = sqrt(v(:,1).^2 + v(:,2).^2);
            angle = atan2(v(:,2),v(:,1));
            x = [x speed angle]; x = [0 x(1,2:end); x];
            x = compute_dataInterpolation(x,(dt:dt:T)');
            v = x(:,4);
            x = x(:,2:3);

        %plot
            % subplot(211)
            % plot(dt:dt:80,v(1:80/dt),'k','linewidth',2)
            % xlabel('time (sec)'), ylabel('rat speed (cm/s)'), set(gca,'fontsize',14), set(gcf,'color','w'), pbaspect([4 1 1])
            % 
            % subplot(212)
            % cplot(x(:,1),x(:,2),v,'linewidth',2), colormap(flipud(parula)),
            % axis square, set(gca,'fontsize',14), xlabel('x (m)'), ylabel('y (m)')
            % cb = colorbar(); ylabel(cb,'rat speed','Rotation',270)
            % 
            % keyboard
    end

    % plot(x)
    % hold on
    % plot(abs(diff(v)/dt))
    % hold off
    % keyboard


%frequency of sequence generation
    %constant theta freq 
        % g_theta = (sign(sin(2*pi*(dt:dt:T)'*(1/T_theta_max)-pi/8-pi/4)+0.4)+1)/2; 
        % g_theta_FF = (sign(sin(2*pi*(dt:dt:T)'*(1/T_theta_max)+pi/2-pi/4)-0.9)+1)/2; 

    % %smooth changes to theta frequence
    %     f_theta = v/v_max*(1/T_theta_max - 1) + T_replay_min; 
    %     f_theta_cumsum = cumsum(f_theta);
        % g_theta = (square(2*pi*dt*f_theta_cumsum,70)+1)/2; 
        % g_theta_FF = (square(2*pi*dt*f_theta_cumsum+0.5,15)+1)/2; 

    %ITI sampling
        %nonstationary distribution: sample theta ITI's based on running speed
        if spatialDim == 1
            sigma_theta = 0.5;
            f_theta = v/max(v)*(1/T_theta_max) + (1-v/max(v))*(1/T_replay_min);
            isi = 1./normrnd(f_theta(1),sigma_theta);
            while isi<0
                isi = 1./normrnd(f_theta(1),sigma_theta);
            end
            t = 1+isi;
            Theta = t;
            while t<T
                isi = 1./normrnd(interp1(dt:dt:T,f_theta,t),sigma_theta);
                while isi<0
                    isi = 1./normrnd(interp1(dt:dt:T,f_theta,t),sigma_theta);
                end
                t = t+isi;
                Theta = [Theta; t];
            end
            Theta(isnan(Theta)) = [];

            % keyboard
        end

        %stationary distribution
        if spatialDim == 2
            sigma_theta = 1.2;%0.6;
            ITI = 1./normrnd((1/T_theta_max),sigma_theta,ceil(T/T_theta_max),1);
                ITI(ITI<=0) = 1/T_theta_max;
            Theta = cumsum(ITI);

            % histogram(ITI,linspace(0,0.25,20),'FaceColor','k','EdgeColor','none','Normalization','probability')
            % axis square, xlabel('theta period (sec)'), ylabel('fraction'), set(gca,'fontsize',14), set(gcf,'color','w')
            % xlim([0 0.25])
            % keyboard
        end

        %theta modulation
            if spatialDim == 1
                theta = [Theta(1:end-1) Theta(1:end-1) + 0.6*diff(Theta)];
            else
                theta = [Theta(1:end-1) Theta(1:end-1) + 0.5*diff(Theta)];
            end
            theta_ff = [Theta(1:end-1) - 0.1*diff(Theta) Theta(1:end-1) + 0.1*diff(Theta)];
            
                %fix negative theta windows
                ind = find(theta(:,1)<0); if ~isempty(ind), if ind>1, theta(ind,1) = theta(ind-1,2); else, theta(ind,1) = 0; end, end
                ind = find(theta_ff(:,1)<0); if ~isempty(ind), if ind>1, theta_ff(ind,1) = theta_ff(ind-1,2); else, theta_ff(ind,1) = 0; end, end
            
            theta_ind = floor(theta/dt);
            theta_ff_ind = floor(theta_ff/dt);
            g_theta = zeros(T/dt,1);
            g_theta_FF = zeros(T/dt,1);
            for i = 1:size(theta,1)
                g_theta(theta_ind(i,1):theta_ind(i,2)) = 1;
                g_theta_FF(theta_ff_ind(i,1):theta_ff_ind(i,2)) = 1;
            end
            g_theta = g_theta(1:T/dt);
            g_theta_FF = g_theta_FF(1:T/dt);

            %used to modulate facilitation
            g_rest = zeros(T/dt,1);
            g_rest(v<0.01) = 1;


    %plot
        % subplot(211)
        % plot(dt:dt:1,g_theta(1:1/dt),'k','linewidth',2)
        % hold on, plot(dt:dt:1,g_theta_FF(1:1/dt),'r','linewidth',2), hold off
        % legend('recurrent','place'), ylabel('modulation'), xlabel('time (sec)'), box on, set(gca,'fontsize',14), axis tight, pbaspect([4 1 1]), set(gcf,'color','w')
        % 
        % subplot(223)
        % cplot(x(:,1),x(:,2),g_theta_FF,'linewidth',2), colormap((parula))
        % axis square, set(gca,'fontsize',14), xlabel('x (m)'), ylabel('y (m)')
        % cb = colorbar(); ylabel(cb,'Place input modulation','Rotation',270)
        % 
        % subplot(224)
        % cplot(x(:,1),x(:,2),g_theta,'linewidth',2), colormap((parula))
        % axis square, set(gca,'fontsize',14), xlabel('x (m)'), ylabel('y (m)')
        % cb = colorbar(); ylabel(cb,'Recurrent input modulation','Rotation',270)
        % 
        % keyboard

%synaptic weight matrics
    if spatialDim==2
        w = mvnpdf([X(:) Y(:)],(n/2+1)/n*[1 1],sig_rec*eye(2)); w = reshape(w,n,n)/max(w); w = w_rec*w;
        % w = 0.5*mvnpdf([X(:) Y(:)],(n/2+1)/n*[1 1],sig_rec*eye(2)/2) - 0.95*mvnpdf([X(:) Y(:)],(n/2+1)/n*[1 1],sig_rec*eye(2)); w = reshape(w,n,n)/max(w); w = w_rec*w;
        % plot(w(:,end/2)), title(sum(w(:,end/2)))
     
        % I1 = 1/2*(1+sin(1*2*pi*X*3 + 0*2*pi*Y*3)) ; 
        % I2 = 1/2*(1+sin(cos(-120)*-2*pi*X*3 + sin(-120)*2*pi*Y*3)) ; 
        % I3 = 1/2*(1+sin(cos(120)*-2*pi*X*3 + sin(120)*2*pi*Y*3)) ; 
        % imagesc(I1+I2+I3)

        if use_periodic == 1
            w_ft=fft2(fftshift(w));
        else
            w_ft=fft2((w),big,big);
        end
    else
        w = normpdf(X,(n/2+1)/n,sig_rec); w = w/max(w); w = w_rec*w;
        if use_periodic == 1
            w_ft=fft(fftshift(w));
        else
            w_ft=fft((w),big);
        end
    end


%initialize vectors
    activityCOM = nan(T/dt,spatialDim);
    activitySpread = nan(T/dt,1);
    if spatialDim==2
        r = zeros(n,n); a = zeros(n,n);
    else
        r = zeros(n,1); a = zeros(n,1); f = zeros(n,1);
    end
    spk_count = zeros(N,1);
    if use_plot == 1, spk_mat = zeros(N,T/dt); end

    spk_mat = nan(n,T/dt);
    r_mat = nan(n,T/dt);
    a_mat = nan(n,T/dt);
    f_mat = nan(n,T/dt);

%simulation
    for t=1:T/dt
        % if mod(t,10*T_plot/dt)==0, t/(T/dt), end
        
        %combined input conductance
            if spatialDim==2
                if use_periodic == 1
                    g_rec = real(ifft2(fft2(r).*w_ft));
                else
                    g_rec = real(ifft2(fft2(r,big,big).*w_ft)); g_rec = g_rec(n/2+1:big-n/2,n/2+1:big-n/2);
                end
                g_FF = reshape(mvnpdf([X(:) Y(:)],x(t,:),sig_FF*eye(2)),n,n); g_FF = w_FF*g_FF/max(g_FF(:));
            else
                if use_periodic == 1
                    g_rec = real(ifft(fft(r).*w_ft));
                else
                    g_rec = real(ifft(fft(r,big).*w_ft)); g_rec = g_rec(n/2+1:big-n/2);
                end
                g_FF = normpdf(X,x(t,:),sig_FF); g_FF = w_FF*g_FF/max(g_FF(:));
            end
            if use_feedforward==0, g_FF = 0; end
            g_feedbackInh = -w_inh*sum(g_rec(:));
            g_unifInput = beta_0;
            if spatialDim == 1
                g_adapt = -g_transition(t)*w_adapt*a;
            else
                g_adapt = -w_adapt*a;
            end

            if use_theta==1
                G = g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_theta_FF(t)*g_FF;
                if use_facilitation == 1 && spatialDim == 1
                    g_facil = w_facil*f;
                    G = g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_theta_FF(t)*g_FF + g_rest(t)*g_theta(t)*A.*g_facil;
                    % G = g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_theta_FF(t)*g_FF + g_theta(t)*A.*g_facil;
                    % G = g_transition(t)*(g_theta(t)*A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh)) + g_theta_FF(t)*g_FF + g_rest(t)*g_theta(t)*A.*g_facil;
                end
            else
                G = A.*(g_rec + g_adapt + g_unifInput + g_feedbackInh) + g_FF;
            end
    
        %pass conductance variables through nonlinearity to generate output rates, F
            F = 0.*(G<0) + G.*(G>=0);   %rectified-linear transfer function
            
        %track bump properties
            activitySpread(t) = compute_imageSpread(F,2);
      
            %circular center of mass
                bin_centers = linspace(-pi,pi,n);
                center_of_mass = nan(1,spatialDim);
                for i = 1:spatialDim
                    if spatialDim==2
                        hist_data = nanmean(r,i);
                    else
                        hist_data = r;
                    end
                    hist_data = hist_data / sum(hist_data);   
                    weighted_mean_sin = sum(dot(hist_data,sin(bin_centers)));
                    weighted_mean_cos = sum(dot(hist_data,cos(bin_centers)));            
                    center_of_mass(i) = (atan2(weighted_mean_sin, weighted_mean_cos)+pi)*n/(2*pi);
                end
                activityCOM(t,:) = center_of_mass;

        %update neural activities 
            if use_spiking == 1
                %spikes generated from inhomogeneous Poisson process
                % spk = poissrnd(F*dt);
    
                %subdivide interval m times and take mth spike (for generating spikes with CV = 1/sqrt(m) )
                spk_sub = poissrnd(repmat(F(:),1,m)*dt);
                spk_count = spk_count+sum(spk_sub,2);
                spk = floor(spk_count/m);
                spk_count = rem(spk_count,m);
                if spatialDim==2, spk = reshape(spk,n,n); end
                if use_plot == 1, spk_mat(:,t) = spk(:); end
                
                %update firing rates
                r = r + spk - r*dt/tau_s;
                
                %update adaptation dynamics
                a = a + spk - a*dt/tau_adapt;

                if use_facilitation == 1 && spatialDim == 1
                    f = f + spk*(1-g_rest(t)) - f*dt/tau_facil; 
                    % f = f + spk - f*dt/tau_facil; 
                end
                
            else
                %update firing rates
                r = r + F*dt - r*dt/tau_s;
                
                %update adaptation dynamics
                a = a + F*dt - a*dt/tau_adapt;
            end

            if spatialDim == 1
                spk_mat(:,t) = spk;
                r_mat(:,t) = r;
                a_mat(:,t) = a;
                f_mat(:,t) = f;
            end

            % if use_plot==1 &&  mod(t,T_plot/dt)==0% && t>1/dt  
            %     neuron = 20;
            % 
            %     set(gcf,'color','w')
            %     subplot(6,1,[1 2]),   
            %         activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread) | [compute_sequenceJumps(activityCOM_sub); 0]>1) = nan;
            %         plot(dt:dt:t*dt,activityCOM_sub(1:t),'k-','linewidth',2)
            %         hold on, plot(dt:dt:t*dt,n*x(1:t)-0.5,'r.','linewidth',4), hold off
            %         hold on, hline(neuron,'g:'), hold off
            %         ylabel('neuron'), 
            %         set(gca,'fontsize',14,'xticklabel',[])
            % 
            %     subplot(613), 
            %         tt = dt:dt:t*dt;
            %         vline_efficient(tt(spk_mat(neuron,1:t)>0),[0 1],[])
            %         ylabel('spikes')
            %         ylabel('\sigma^{spk}(t)')
            %         set(gca,'xticklabel',[],'fontsize',14)
            % 
            %     subplot(614), 
            %         plot(dt:dt:t*dt,r_mat(neuron,1:t),'k','linewidth',2)
            %         ylabel('s(t)')
            %         set(gca,'xticklabel',[],'fontsize',14)
            % 
            %     subplot(615), 
            %         plot(dt:dt:t*dt,a_mat(neuron,1:t),'k','linewidth',2)
            %         ylabel('a(t)')
            %         set(gca,'xticklabel',[],'fontsize',14)
            % 
            %     subplot(616), 
            %         plot(dt:dt:t*dt,f_mat(neuron,1:t),'k','linewidth',2)
            %         ylabel('f(t)'), xlabel('time (sec)')
            %         set(gca,'fontsize',14)
            % 
            %     set(findall(gcf,'Type','axes'),'xlim',[0 t*dt])
            %     % set(findall(gcf,'Type','axes'),'xlim',[44 104])
            %     drawnow
            % end
            % continue
        
        %Plot
        if use_plot==1 &&  mod(t,T_plot/dt)==0% && t>1/dt    
            if spatialDim == 1
                subplot(331), plot(a)
                subplot(331), plot(r,'k','linewidth',2), ylim([0 3]), vline(n*x(t),'r'), xlabel('neuron'), ylabel('firing rate')
                subplot(334), plot(a,'k','linewidth',2), vline(n*x(t),'r'), xlabel('neuron'), ylabel('adaptation'), %ylim([0 30])
                    if w_facil ~=0, hold on, plot(f,'b','linewidth',2), hold off, end
                % subplot(337), plot(dt:dt:T,v,'k','linewidth',2), hold on, plot(t*dt,v(t),'ro'), hold off, ylabel('speed'), xlabel('time') 
                subplot(1,3,[2 3])
                    activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread) | [compute_sequenceJumps(activityCOM_sub); 0]>1) = nan;
                    plot(dt*(1:t),activityCOM_sub(1:t),'k-','linewidth',2)
                    hold on, plot(t*dt,n*x(t)-0.5,'ro','linewidth',2), hold off
                    hold on, plot((1:t)*dt,n*x(1:t)-0.5,'r.','markersize',4), hold off
                    ylim([0 n]), xlabel('time (sec)'), ylabel('neuron')
                set(findobj(gcf,'type','axes'),'FontSize',12), set(gcf,'color','w')
            else
                % subplot(331), imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('firing rate'), clim([0 3])
                % subplot(334), imagesc(a), set(gca,'ydir','normal'), axis off, axis square, title('adaptation (-)'), clim([0 30])
                % subplot(337), plot(dt:dt:T,v), hold on, plot(t*dt,v(t),'ro'), hold off, ylabel('speed'), xlabel('time') 
                subplot(231), imagesc(r), set(gca,'ydir','normal'), axis off, axis square, title('Firing rate'), clim([0 3])
                subplot(234), imagesc(a), set(gca,'ydir','normal'), axis off, axis square, title('Adaptation'), clim([0 30])
                subplot(1,3,[2 3]), 
                    activityCOM_sub = activityCOM; activityCOM_sub(isnan(activitySpread(:,1)) | [compute_sequenceJumps(activityCOM_sub); 0]>1,:) = nan;
                    plot((activityCOM_sub(1:t,1)),(activityCOM_sub(1:t,2)),'.','color',0.8*ones(1,3),'markersize',8)
                    % hold on, cplot((activityCOM_sub(t-1/dt+1:t,1)),(activityCOM_sub(t-1/dt+1:t,2)),1:1/dt,'.','markersize',12), hold off
                    hold on, cplot((activityCOM_sub(1:t,1)),(activityCOM_sub(1:t,2)),1:t,'.','markersize',12), hold off
                    hold on, plot((n*x(t,1)-0.5),(n*x(t,2)-0.5),'ro','linewidth',2), hold off
                    hold on, plot((n*x(1:t,1)-0.5),(n*x(1:t,2)-0.5),'k:','linewidth',2), hold off
                    axis([0 n 0 n]), set(gca,'xtick',[],'ytick',[]), 
                    axis square, title('Center-of-mass trajectory')
                    % camroll(-90)

                % set(findobj(gcf,'type','axes'),'FontSize',12), set(gcf,'color','w')
                % xl = xlim; yl = ylim; text(xl(1)+0.01*diff(xl),yl(1),strcat('elapsed time:',{' '},num2str(t*dt,'%0.2f'),' sec'),'HorizontalAlignment','left','VerticalAlignment','top','fontsize',14);
                % writeVideo(V,getframe(gcf));
            end
            drawnow
        end
            
    end
    
    if use_plot == 1
        t = dt*(1:T/dt)';
        spk_cell = cell(N,1);
        for i = 1:N
            spk_cell{i} = t(spk_mat(i,:)==1);
            spk_cell{i} = spk_cell{i} + 0.1*(2*rand(size(spk_cell{i}))-1)/2;
        end
        plot_raster_times(spk_cell)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %rat position
            t = dt*(1:T/dt)';
            if spatialDim == 1, angle = [sign(diff(x)); 0];
            else, angle = [atan2(diff(x(:,2)),diff(x(:,1))); nan];
            end
            rat = [t x v angle];

        %Theta cycles
            [pks,locs] = findpeaks(abs(diff(g_theta)));
            theta_bounds = [locs(1:2:end-1)+1 locs(2:2:end)+1];
            theta_timeBounds = t(theta_bounds);
            theta_time = nanmean(theta_timeBounds,2);
            numThetaCycles = size(theta_bounds,1);

            [pks,locs] = findpeaks(abs(diff(g_theta_FF)));
            theta_bounds_FF = [locs(1:2:end-1)+1 locs(2:2:end)+1]; 

        %rat properties
            if spatialDim == 1, id_angleColumns = []; else, id_angleColumns = 5; end
            theta_rat = compute_dataInterpolation(rat,theta_time,id_angleColumns);

            [pks,locs] = findpeaks(abs(diff(g_theta_FF)));
            theta_bounds_FF = [locs(1:2:end-1)+1 locs(2:2:end)+1]; 
            theta_timeBounds_FF = t(theta_bounds_FF);
            theta_rat_start = compute_dataInterpolation(rat,theta_timeBounds_FF(:,1),id_angleColumns);

        %Rat running periods
            v_thr = 0.03;
            d = 0*(v>v_thr) + 1*(v<v_thr);
            [~,locs] = findpeaks(abs(diff(d)));
            rat_runningTime_start = t(locs(1:2:end));
            rat_runningTime_end = t(locs(2:2:end));
            rat_running_timeBounds = [rat_runningTime_start(1:min(length(rat_runningTime_start),length(rat_runningTime_end))) rat_runningTime_end(1:min(length(rat_runningTime_start),length(rat_runningTime_end)))];
            numRunningPeriods = size(rat_running_timeBounds,1);

            %time since rat running
            theta_timeSinceRunning = theta_time - repmat(rat_running_timeBounds(:,1)',numThetaCycles,1);
            theta_timeSinceRunning(theta_timeSinceRunning<0) = nan;
            theta_timeSinceRunning = min(theta_timeSinceRunning,[],2);
        
        %Rat stopping periods
            v_thr = 0.01;
            d = 0*(v>v_thr) + 1*(v<v_thr);
            [~,locs] = findpeaks(abs(diff(d)));
            rat_stoppingTime_start = t(locs(2:2:end-1));
            rat_stoppingTime_end = t(locs(3:2:end));
            rat_stopping_timeBounds = [rat_stoppingTime_start rat_stoppingTime_end];
            numStoppingPeriods = size(rat_stopping_timeBounds,1);

            %time since rat stopping
            theta_timeSinceStopping = theta_time - repmat(rat_stopping_timeBounds(:,1)',numThetaCycles,1);
            theta_timeSinceStopping(theta_timeSinceStopping<0) = nan;
            theta_timeSinceStopping = min(theta_timeSinceStopping,[],2);

            % plot(rat(:,1),rat(:,2))
            % vline(rat_running_timeBounds(:,1),'m--')
            % vline(rat_running_timeBounds(:,2),'b--')            
            % vline(rat_stopping_timeBounds(:,1),'r')
            % vline(rat_stopping_timeBounds(:,2),'g')
            % keyboard

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if use_plot == 0
            data_params(kk).T = T;
            data_params(kk).dt = dt;
            data_params(kk).spatialDim = spatialDim;
            data_params(kk).n = n;
            data_params(kk).x = x;
            data_params(kk).theta_bounds = theta_bounds;
            data_params(kk).theta_bounds_FF = theta_bounds_FF;
            data_params(kk).rat_stopping_timeBounds = rat_stopping_timeBounds;
            data_params(kk).rat_running_timeBounds = rat_running_timeBounds;
            data_params(kk).activityCOM = activityCOM;
            data_params(kk).T_transition = T_transition;
            data_params(kk).T_run = T_run;
            data_params(kk).T_stop = T_stop;
            data_params(kk).v_max = v_max;
            data_params(kk).T_theta_max = T_theta_max;
            data_params(kk).w_adapt = w_adapt;
            data_params(kk).tau_adapt = tau_adapt;
            if spatialDim == 1
                data_params(kk).w_facil = w_facil; 
                data_params(kk).tau_facil = tau_facil; 
            end

        end
end


%plot
for kk = 1:2

    if isempty(data_params(kk).T), continue, end

    T = data_params(kk).T;
    dt = data_params(kk).dt;
    spatialDim = data_params(kk).spatialDim;
    n = data_params(kk).n;
    x = data_params(kk).x;
    theta_bounds = data_params(kk).theta_bounds;
    theta_bounds_FF = data_params(kk).theta_bounds_FF;
    rat_stopping_timeBounds = data_params(kk).rat_stopping_timeBounds;
    rat_running_timeBounds = data_params(kk).rat_running_timeBounds;
    activityCOM = data_params(kk).activityCOM;
    T_transition = data_params(kk).T_transition;
    T_run = data_params(kk).T_run;
    T_stop = data_params(kk).T_stop;
    v_max = data_params(kk).v_max;
    T_theta_max = data_params(kk).T_theta_max;
    w_adapt = data_params(kk).w_adapt;
    tau_adapt = data_params(kk).tau_adapt;
    w_facil = data_params(kk).w_facil;
    tau_facil = data_params(kk).tau_facil;

    %rat position
        t = dt*(1:T/dt)';
        angle = sign(diff(x)); angle = [angle(1,:); angle];
        speed = diff(x(:,1))/mode(diff(t)); speed = [speed(1,:); speed];
        rat = [t x speed angle];

    %Theta cycles
        theta_timeBounds = t(theta_bounds);
        theta_time = theta_timeBounds(:,1);%nanmean(theta_timeBounds,2);
        numThetaCycles = size(theta_bounds,1);

    %rat properties
        theta_rat = compute_dataInterpolation(rat,theta_time,[]);

        theta_timeBounds_FF = t(theta_bounds_FF);
        theta_rat_start = compute_dataInterpolation(rat,theta_timeBounds_FF(:,1),[]);

    %Rat running and stopping periods
        numRunningPeriods = size(rat_running_timeBounds,1);
        numStoppingPeriods = size(rat_stopping_timeBounds,1);

        %time since rat running
        theta_timeSinceRunning = theta_time - repmat(rat_running_timeBounds(:,1)',numThetaCycles,1);
        theta_timeSinceRunning(theta_timeSinceRunning<0) = nan;
        theta_timeSinceRunning = min(theta_timeSinceRunning,[],2);

        %time since rat stopping
        theta_timeSinceStopping = theta_time - repmat(rat_stopping_timeBounds(:,1)',numThetaCycles,1);
        theta_timeSinceStopping(theta_timeSinceStopping<0) = nan;
        theta_timeSinceStopping = min(theta_timeSinceStopping,[],2);

    %sequence direction
        theta_seqDir = nan(numThetaCycles,1);
        theta_seqDist = nan(numThetaCycles,1);
        for i = 2:numThetaCycles
            theta_sub = unwrap(activityCOM(theta_bounds(i,1):theta_bounds(i,2),:));
            
            delta = diff(theta_sub,[],1);
            delta(prod(delta,2)==0,:) = [];
            theta_seqDir(i,1) = sign(nanmean(delta));
            theta_seqDist(i) = abs(theta_sub(end) - theta_sub(1));

            % hold on, plot(t(theta_bounds(i,1):theta_bounds(i,2)),theta_sub-(n*rat(theta_bounds(i,1),2)-0.5))
            % % title(theta_seqDist(i))
            % title(theta_seqDir(i,1))
            % drawnow, %keyboard
        end

    theta_seqDist_thr = 3;
    %theta_seqDist_thr = 3; % in manuscript

    %forward vs reverse rate
        [edges_time,centers_time] = load_timeBins2([0 12],0.5,0.5);
        edges_distance = linspace(0.15,0.46,20); centers_distance = edges_distance(1:end-1) + mean(diff(edges_distance))/2;
        data_distance = nan(numStoppingPeriods,length(edges_distance)-1);
        data_forward = nan(numStoppingPeriods,length(edges_time)-1);
        data_reverse = nan(numStoppingPeriods,length(edges_time)-1);
        data_difference = nan(numStoppingPeriods,length(edges_time)-1);
        data_sum = nan(numStoppingPeriods,length(edges_time)-1);
        data_eventRate = nan(numStoppingPeriods,length(edges_time)-1);
        data_length = nan(numStoppingPeriods,length(edges_time)-1);
        for i = 1:numStoppingPeriods
            data_sub = compute_dataTemporalConcatenation([theta_time theta_timeSinceStopping theta_seqDir theta_rat(:,4) theta_seqDist],rat_stopping_timeBounds(i,:));
            data_sub(data_sub(:,end)<theta_seqDist_thr,:) = [];
            if size(data_sub,1)<2, continue, end

            %distance
                data_distance(i,:) = histcounts(data_sub(:,5),edges_distance,'normalization','pdf');

            for j = 1:length(centers_time)
                %forward vs. reverse rate
                ind = find(data_sub(:,2)>edges_time(j,1) & data_sub(:,2)<=edges_time(j,2));
                if length(ind)>0
                    data_forward(i,j) = sum(data_sub(ind,3)==1)/mean(diff(edges_time,[],2));
                    data_reverse(i,j) = sum(data_sub(ind,3)==-1)/mean(diff(edges_time,[],2));
                    data_difference(i,j) = data_forward(i,j) - data_reverse(i,j);
                    data_sum(i,j) = data_forward(i,j) + data_reverse(i,j);
                    data_length(i,j) = nanmean(data_sub(ind,5));
                    data_eventRate(i,j) = length(ind);
                else
                    data_forward(i,j) = 0;
                    data_reverse(i,j) = 0;
                    data_eventRate(i,j) = 0;
                    data_difference(i,j) = 0;
                    data_sum(i,j) = 0;
                end
            end

            data_forward(i,:) = conv(data_forward(i,:),setUp_gaussFilt([1 10],1),'same');
            data_reverse(i,:) = conv(data_reverse(i,:),setUp_gaussFilt([1 10],1),'same');
            data_difference(i,:) = conv(data_difference(i,:),setUp_gaussFilt([1 10],1),'same');
            data_sum(i,:) = conv(data_sum(i,:),setUp_gaussFilt([1 10],1),'same');
        end

        % p_forVsRev = nan(1,length(centers_time));
        % for j = 1:length(centers_time)
        %     [~,p_forVsRev(j)] = ttest(data_forward(:,j),data_reverse(:,j));
        % end

    %boostraps and shuffles
        [~,theta_time_nan] = compute_dataTemporalConcatenation([theta_time theta_time],rat_stopping_timeBounds);
        ind_thetaCycles = find(theta_seqDist>theta_seqDist_thr & ~isnan(theta_time_nan(:,2)));
        replayEvents = [theta_time(ind_thetaCycles) theta_seqDir(ind_thetaCycles)];
        [edges_time_pairwise,centers_time_pairwise] = load_timeBins2([0 10],0.1,0.4);

        %mean
        delT = squareform(pdist(replayEvents(:,1))); delT = delT(:);
        delX = squareform(pdist(replayEvents(:,2))); delX = delX(:)/2;
        data_opposite = nan(1,length(centers_time_pairwise));

        %bootstraps
        % numBoostraps = 500;
        numBootstraps = 5;
        data_bootstraps = nan(numBoostraps,length(centers_time_pairwise));
        for j = 1:length(centers_time_pairwise)
            ind = find(delT>edges_time_pairwise(j,1) & delT<=edges_time_pairwise(j,2));
            data_opposite(j) = nanmean(delX(ind));

            for i = 1:numBoostraps
                data_bootstraps(i,j) = nanmean(datasample(delX(ind),length(ind)));
            end
        end
        bootstraps_margins = quantile(data_bootstraps,[0.0250 0.975]);

        % %shuffles
        % numShuffles = 500;
        % data_shuffles = nan(numShuffles,length(centers_time_pairwise));
        % for i = 1:numShuffles
        %     i
        %     replayEvents_shuffle = replayEvents;
        %     replayEvents_shuffle(:,1) = replayEvents(randperm(size(replayEvents,1))',1);
        %     delT = squareform(pdist(replayEvents_shuffle(:,1))); delT = delT(:);
        %     delX = squareform(pdist(replayEvents_shuffle(:,2))); delX = delX(:)/2;
        % 
        %     delX(delT>max(edges_time_pairwise(:))) = [];
        %     delT(delT>max(edges_time_pairwise(:))) = [];
        %     for j = 1:length(centers_time_pairwise)
        %         ind = find(delT>edges_time_pairwise(j,1) & delT<=edges_time_pairwise(j,2));
        %         data_shuffles(i,j) = nanmean(delX(ind));
        %     end
        % end
        % shuffles_margins = quantile(data_shuffles,[0.0250 0.975]);

    %plot
        % figure
        colors = [.4660 0.6740 0.1880;0.4940 0.1840 0.5560];
        subplot(2,3,1+3*(kk-1))
        shadedErrorBar(centers_time,nanmean(data_forward),nanstd(data_forward)./sqrt(size(data_forward,1)),'lineprops',{'color',colors(1,:),'linewidth',2})
        hold on, shadedErrorBar(centers_time,nanmean(data_reverse),nanstd(data_reverse)./sqrt(size(data_reverse,1)),'lineprops',{'color',colors(2,:),'linewidth',2}), hold off
        xlabel('Time since arrival (s)'), ylabel('Events/sec')
        axis square, set(gca,'fontsize',14), set(gcf,'color','w'), box on, axis tight
        % hold on, plot(centers_time(p_forVsRev<0.05),0.7,'k.','markersize',10), hold off, legend('forward','reverse')
        xlim([0 10])
        % ylim([0 0.8])

        subplot(2,3,2+3*(kk-1))
        shadedErrorBar(centers_time,nanmean(data_difference),nanstd(data_difference)./sqrt(size(data_difference,1)),'lineprops',{'k','linewidth',2})
        hold on, shadedErrorBar(centers_time,nanmean(data_sum),nanstd(data_sum)./sqrt(size(data_sum,1)),'lineprops',{'r','linewidth',2}), hold off
        xlabel('Time since arrival (s)'), ylabel('Events/sec'),
        axis square, set(gca,'fontsize',14), set(gcf,'color','w'), box on, axis tight
        hline(0,'k')
        % hold on, plot(centers_time(p_forVsRev<0.05),0.4,'k.','markersize',10), hold off
        xlim([0 10])
        ylim([-0.5 1.2])

        subplot(2,3,3+3*(kk-1))
        shadedErrorBar(centers_time_pairwise,100*data_opposite,100*abs(bootstraps_margins-data_opposite),'lineprops',{'k','linewidth',2})
        % hold on, plot(centers_time_pairwise,100*shuffles_margins,'k--'), hold off
        xlabel('time between events (sec)'), ylabel('% opposite'), set(gca,'fontsize',14)
        axis square, box on, axis tight
        axis([0 10 0 100])
        % hold on, plot(centers_time_pairwise(data_mean>shuffles_margins(2,:) | data_mean<shuffles_margins(1,:)),90,'k.','markersize',10), hold off
        drawnow

end



