%...................objective function ..................................

function f_sum = objective_function_21_9_20(x,G_pa,omega)
global M W
BW = 5*10^6;                 %System bandwidth
Bc = BW/M;

%fprintf('\n q_opt = ');
%fprintf('%g ', x);

h_total = zeros(M,1);
for m=1:M
    %h1=0; h2=0;
    h1 = W(m,1)*Bc*log(1+omega(m)*G_pa(m,1));
    h2 = W(m,2)*Bc*log( (x(m)*G_pa(m,2)+1) / (omega(m)*G_pa(m,2)+1) );
    %fprintf('\n h1 = ');
    %fprintf('%g ', h1);
    %fprintf('            h2 = ');
    %fprintf('%g ', h2);
    h_total(m) = -h1-h2;
    %fprintf('\n htotal = ');
    %fprintf('%g ', h_total);
end
f_sum = sum(h_total);

%fprintf('\n f_sum = ');
%fprintf('%g ', -f_sum);
end