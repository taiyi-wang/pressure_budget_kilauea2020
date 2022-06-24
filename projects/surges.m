
function varargout = surges(fn,varargin)
    [varargout{1:nargout}] = feval(fn,varargin{:});
    % Code to model surges in LERZ during 2018 Kilauea caldera collapse.
    % Forward model includes a single reservoir in the ERZ.  
    %  Matrix form where the pressure difference p(t=0) - p(t=T) = 
    %  [dP; 0]; that is the pressure in the ERZ reservoir returns to the
    %  initial value at the end of the cycle. While the pressure in the HMM
    %  reservoir changes by dP at the start of a cycle due to collapse.
    %  The model is described by three decay times.
    %  Inversion solves for the best fitting values of these decay times by

    %  such that:
    %   1)  P1  decays exponentially with a decay time of 0.5 days
    %   2)  Peak pressure in P2 is delayed by 2-3 hours
    %   3)  Max(P2)/min(P2) ~ 3

    % To run with nominal weights in the objective: 
    % [tau, fval] =  surges('invert');
    %
    % To run with other weights, for example 
        %   sigp = 0.001; % sigma on fitting summit pressure
        %   sigt = 1.5;   % sigma on time of peak flux
        %   sigq = 0.2;   % sigma on flux ratio
        %   sigmas = [sigp,sigt,sigq];
        %
        %   [tau, fval] =  surges('invert', sigmas);    

    % Paul Segall, 2021-2022
end

function [tau, fval] = invert(varargin)
    tau0 = [1/24,1.5/24,0.25]; % initial guess at decay times
   
    opts = optimoptions(@lsqcurvefit,'FunctionTolerance',1e-3,'Display','off');

    % function describes exponential decay of summit pressure with decay
    % time of 0.5 days. This is used in the objective function
    fun = @(x,xdata)x(1)*exp(-xdata/0.5) + x(2);

    % Run optimization
    if nargin == 0  % use nominal weights
    [tau,fval] = fmincon(@(x) surge_obj(x, opts, fun),tau0,[],[]);

    % print the values in the objective function
    surges('surge_obj',tau, opts, fun,1);

    elseif nargin == 1
    sigmas = varargin{1};
    [tau,fval] = fmincon(@(x) surge_obj(x, opts, fun, 0, sigmas),tau0,[],[]);

    % print the values in the objective function
    surge_obj(tau, opts, fun, 1, sigmas);
    end
end
 

function obj = surge_obj(varargin)
    %function obj = surge_obj(tau, opts, fun)
    % Objective function for the optimization

    % Input parameters
    tau = varargin{1};
    opts = varargin{2};
    fun = varargin{3};

    if nargin > 3
        plt = varargin{4}; % = 1 to plot solution
    else
        plt = 0;
    end

    if nargin > 4
        sigmas = varargin{5};
        sigp = sigmas(1);
        sigt = sigmas(2);
        sigq = sigmas(3); 
    else
        % Nominal Values of sigmas in objective
        sigp = 0.05;
        sigt = 0.5;
        sigq = 0.2;
    end

    
    % Hardwired parameters
    T = 1.4;  % Cycle time in days
    dP = 3;  % pressure change in MPa
    
    %unpack decay constants
    tau1 = tau(1); tau2 = tau(2); tau3 = tau(3);
    t = linspace(0,T,1000);
    
    % create non-dimensional matrix
    A = [-1/tau3  1/tau3;
        1/tau2    -1/tau1];
    
    [V,D] = eig(A);
    
    rgH = 5; % rho*g*H; where H is depth of reservoir (MPa)
    B = rgH*[1/tau3; -1/tau2];
    
    % particular solution
    pp = -inv(A)*B;
    
    % solve for constants
    C = (V - [V(:,1)*exp(D(1,1)*T) V(:,2)*exp(D(2,2)*T)])\[dP; 0];
    
    % full pressure solution
    p = C(1)*V(:,1)*exp(D(1,1)*t) + C(2)*V(:,2)*exp(D(2,2)*t) + pp;
    
    % Find time and value of maximum pressure in LERZ reservoir
    [pmax, I] = max(p(2,:));
    % display(['pressure ratio: ',num2str(pmax/p(2,1))])
    % display(['time of peak flux: ',num2str(t(I)*24)])
    
    if(plt)
    % Subtract hydrostatic term from p1 only because p2 already is
    % measured relative to hydrostatic
    figure; plot(t*24, p - [rgH;0], 'LineWidth',2); hold on
    set(gca,'FontSize',14)
    end


 % find best fit to exponential decay in summit reservoir  with fun where     
 % fun = @(x,xdata)x(1)*exp(-xdata/0.5) + x(2);
 
    x0 = [3,1];
    [x,~,~,~] = lsqcurvefit(fun,x0,t,p(1,:),[],[],opts);
    ppred = fun(x,t);  % predicted data for exponential decay


    if(plt)
    plot(t*24,ppred-rgH, 'k--', 'LineWidth',2)
        display(['pressure ratio: ',num2str(pmax/p(2,1))])
        display(['time of peak flux (hrs): ',num2str(t(I)*24)])
        display(['Minimum HMM driving pressure (MPa):  ', num2str(p(1,end)-rgH)])
        legend('p_{HMM}', 'p_{ERZ}', 'Exp Fit','AutoUpdate','off')
        xlabel('Time, hrs', 'FontSize',16)
        ylabel('Excess Pressure, MPa', 'FontSize',16)
    plot(t(I)*24, p(2,I), 'ro','LineWidth',2); hold on
   end

    %  difference between observed and predicted summit pressure
    r = p(1,:)-ppred; 
    N = length(r); 
    k = 1:N/20:N;
  
  % components of the objective function Quadratic form
  obj1 = norm(r(k)/sigp,2)/length(k);   % Fit to exponential decay
  obj2 = ( (t(I)*24 - 2.5)/sigt )^2;    % Fit to time of peak flux  
  obj3 = ( (pmax/p(2,1)-3)/sigq )^2;    % Fit to flux ratio
  
  obj = obj1 + obj2 + obj3;
 if(plt)
  disp('  ');
  disp('Components of Objective Function (Exponential fit; Time of peak flux; Flux Ratio');
  obj1, obj2, obj3  
  title([ '\sigma_p, \sigma_t, \sigma_q =   ', num2str(sigp), ', ',  num2str(sigt), ' , ', num2str(sigq)])
 end

end













