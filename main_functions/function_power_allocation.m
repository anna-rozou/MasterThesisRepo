%................Power Allocation & Sum Rate Calculation...................

function [p1, p2, q, SumRate,rate1,rate2 break_const_2] = function_power_allocation_01_11_20(s_match,G_pa,pt,BW)
global M W

N = 2*M;
Bc = BW/M;

%initialization
%G_pa = zeros(M,2);     %users Gamma
%R_pa = zeros(M,1);     %sum rate of each channel
%H_pa = zeros(M,2);     %channel coefficients
p_function = zeros(M,2);
omega = zeros(M,1);
break_const_2 = 0;

%channel coefficients & Gamma calculation
%{
for m1=1:M
    G_pa(m1,1) = G(s_match(m1,1),m1);
    G_pa(m1,2) = G(s_match(m1,2),m1);
end
%}
%omega calculation
for m1=1:M
    o1 = W(m1,2)*G_pa(m1,2) - W(m1,1)*G_pa(m1,1);
    o2 = G_pa(m1,1)*G_pa(m1,2) * (W(m1,1)-W(m1,2));
    omega(m1) = o1/o2;
    if omega(m1)<0
        fprintf('\no1 = : ');
        fprintf('%g ',o1);
        fprintf('         o2 = : ');
        fprintf('%g ',o2);
        fprintf('\n');
    end
        
end
if pt<2*sum(omega)
    fprintf('\npt(b1) = : ');
    fprintf('%g ', pt);
    fprintf('\nsum = : ');
    fprintf('%g ', 2*sum(omega));
    fprintf('\n pt(b1)<2*sum(omega)  -> ERROR ');
    
    break_const_2 = 1;
    
    fprintf('\nG_pa = : ');
    fprintf('%g ',G_pa(1:M,1));
    fprintf('\n         ');
    fprintf('%g ',G_pa(1:M,2));
    
    fprintf('\npt = : ');
    fprintf('%g ',pt);
    
    fprintf('\nomega = : ');
    fprintf('%g ',omega);
    fprintf('\n');
end
if break_const_2 == 0
    %...........fmincon.............
    A = ones(1,M);
    b=pt;
    lb = 2*omega';
    ub = pt*ones(1,M);
    %x0 = zeros(1,M);
    x0 = pt*ones(1,M)/M;
    Aeq = [];
    beq = [];
    options= optimoptions('fmincon','Algorithm','interior-point');
    %,'PlotFcns','optimplotx');
    
    o_f = objective_function_21_9_20(x0,G_pa,omega);
    %{
    fprintf('\nlb = : ');
    fprintf('%g ',lb);
    fprintf('\n');
     fprintf('\nub = : ');
    fprintf('%g ',ub);
    fprintf('\n');
    fprintf('\nobjective_function = : ');
    fprintf('%g ',o_f);
    fprintf('\n');
    %}
    [qval, fval]= fmincon(@(x) (objective_function_21_9_20(x,G_pa,omega)),x0,A,b,Aeq,beq,lb,ub,[],options);
    fval = -fval;
    rate1=zeros(M,1);
    rate2=zeros(M,1);
    %power coefficient and rate calculation
    for m1=1:M
        p_function(m1,1) = omega(m1);
        p_function(m1,2) = qval(m1)-omega(m1);
        %sum rate calculation
        rate1(m1) = W(m1,1)*Bc*log(1+omega(m1)*G_pa(m1,1));
        rate2(m1) = W(m1,2)*Bc*log((qval(m1)*G_pa(m1,2)+1) / (omega(m1)*G_pa(m1,2)+1));
        %R_pa(m1) = rate1+rate2;
    end
else
    fval=0;
    qval=0;
    p_function(1:M,1:2) = 0;
end
%{
%power constraints calculation
for m=1:M
    help_var = W(m,2)*G(m,2)-W(m,1)*G(m,1);
    if help_var>0
        G(m,2) = (W(m,1)*G(m,1) / W(m,2)) - 1;
    end
    %optimal solution of OP2
    omega(m) = (W(m,2)*G(m,2) - W(m,1)*G(m,1)) / ((G(m,1).*G(m,2)) * (W(m,1)-W(m,2)));
    q_function(m) =(W(m,2)*Bc)/lamda - (1/G(m,2));
    if q_function(m)< 2* omega(m)
        q_function(m) = 2* omega(m);
    end
    %optimal solution of OP3
    p_function(m,1) = omega(m);
    p_function(m,2) = q_function(m)-omega(m);
end
%}

%return values
%SumRate= sum(R_pa);
SumRate= fval;
p1 = p_function(:,1);
p2 = p_function(:,2);
q = qval;
end