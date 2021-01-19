%result = squeeze(num2cell(reshape(perms(1:6),[],2,3),2));
%result_2 = squeeze(num2cell(reshape(perms(1:4),[],2,2),2));

clear
clc
global M W G 

M = 4;                       %Number of channels
N = 2*M;                     %Number of users
eta = -1;                    %Path loss exponent

pt =[20];      %15;
repeat_number=1;

dmin_u = 30;              %Minimum distances between users
dmin_b = 40;              %Minimum distances of users from base station (BS)
radious = 300;            %Rdious of the cell

BW = 5*10^6;                 %System bandwidth
Bc = BW/M;
No = -174 ;                  %Noise power spectral density(dBm)
no = 10^((No-30)/10);
sigma_m = (BW*no)/M; 

%USERS Weights
W(1:M,1) = 0.9;     %user 1 weight
W(1:M,2) =  1.1;    %user 2 weight

%.....................smooth result..........................
g = zeros(N,M);     
for n=1:N
    for m=1:M
        g(n,m) = raylrnd(1);
    end
end
dis = function_distance_calc_30_10_20(dmin_u, dmin_b, radious);
H = zeros(N,M);
G = zeros(N,M);
for n=1:N
    for m=1:M
        H(n,m) = g(n,m)*(dis(n) ^eta);
        G(n,m)=(abs(H(n,m)^2)/sigma_m);
    end
end
%.............................................................
%}
[matrix1, gamma, help_matrix] = best_solution_helper_18_11_20();

p1_final=zeros(M,length(pt),length(matrix1));
p2_final=zeros(M,length(pt),length(matrix1));
q_final=zeros(M,length(pt),length(matrix1));
SumRate=zeros(length(pt),length(matrix1));

for b1=1:length(pt)
    %Initialization -> q(m)
    q_match=zeros(M,1);
    for m=1:M
        q_match(m) = pt(b1)/(M);
    end
    %g_calculation,dis_calculation,H,G
    %{
        g = zeros(N,M);     %channel coefficients that follows gaussian distribution
        for n=1:N
            for m=1:M
                g(n,m) = raylrnd(1);
            end
        end
        %dis_calculation
        dis = function_distance_calc_30_10_20(dmin_u, dmin_b, radious);

        H = zeros(N,M);
        G = zeros(N,M);
        for n=1:N
            for m=1:M
                H(n,m) = g(n,m)*(dis(n) ^eta);
                G(n,m)=(abs(H(n,m)^2)/sigma_m);
            end
        end
    %}
 
    %channel assignment -> best pairing
    s_match_total = zeros(M,2);
    gamma_match = zeros(M,2);
    for b2=1:length(matrix1)
        %bla bla bla
        fprintf('___________________________________________________');
        fprintf('');
        fprintf('\npt = ');
        fprintf('%g ', pt(b1));
        fprintf('       b2 = ');
        fprintf('%g ', b2);
        fprintf('\n');
        
        s_match_total=matrix1(:,:,b2);
        gamma_match = gamma(:,:,b2);
        [p1_final(:,b1,b2),p2_final(:,b1,b2), q_final(:,b1,b2), SumRate(b1,b2), break_const_2] = function_power_allocation_01_11_20(s_match_total(:,:),gamma_match(:,:),pt(b1));
    end
end

%.........delete invalid points..........
help_matrix_2 = zeros(M,2,length(matrix1));
SumRate_2 = zeros(length(pt),length(matrix1));
b3=0;
for b2=1:length(matrix1)
    check = 0;
    for b1=1:length(pt)
        if SumRate(b1,b2) == 0
            check = 1;
            break
        end
    end
    if check ==0
        b3=b3+1;
        help_matrix_2(:,:,b3) = matrix1(:,:,b2);
        SumRate_2(:,b3) = SumRate(:,b2);
    end
end
for b4=b3+1:length(matrix1)
    SumRate_2(:,b3) = [];
    help_matrix_2(:,:,b3) = [];
end

%..............find the best combination..............
[maximum_SumRate,idx] = max(SumRate_2,[],2);
best_combination = help_matrix_2(:,:,idx(1));

fprintf('\n Ladies & gendlemets the best combination for 6 users is ... \n');
fprintf('\n *************************');
%fprintf('%g ', best_combination);
best_combination
fprintf('\n *************************');
fprintf('\n');



%plot_all
Xaxis = pt;
Yaxis = SumRate_2/10^6;
plot(Xaxis,Yaxis,':')
title('Solutions of all possible combinations')
%xlabel('Combinations')
xlabel('Power of BS with 6 Users (Watt)')
ylabel('Total Sum Rate of System (Mbps)')

%plot the best
Xaxis = pt;
Yaxis = SumRate_2(:,idx(1))/10^6;
plot(Xaxis,Yaxis,':')
title('Best Solution')
%xlabel('Combinations')
xlabel('Power of BS with 6 Users (Watt)')
ylabel('Total Sum Rate of System (Mbps)')


