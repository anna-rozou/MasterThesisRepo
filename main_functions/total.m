clear
clc
global M W G 

M = 3;                       %Number of channels
N = 2*M;                     %Number of users
eta = -1;                    %Path loss exponent

pt =[15];      %15;
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

p1_final=zeros(M,repeat_number,length(pt));
p2_final=zeros(M,repeat_number,length(pt));
q_final=zeros(M,repeat_number,length(pt));
SumRate_Repeat = zeros(repeat_number,length(pt));
SumRate=zeros(length(pt),1);

%smooth result
g = zeros(N,M);     %channel coefficients that follows gaussian distribution
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
%}

b=1;
while b<=length(pt)
%for b=1:length(pt)
    b2=1;
    while b2<=repeat_number
    %for b2=1:repeat_number
        
        %Initialization -> q(m)
        q_match=zeros(M,1);
        for m=1:M
            q_match(m) = pt(b)/(M);
        end
        %         %g_calculation
        %         g = zeros(N,M);     %channel coefficients that follows gaussian distribution
        %         for n=1:N
        %             for m=1:M
        %                 g(n,m) = raylrnd(1);
        %             end
        %         end
        %dis_calculation
        %        dis = function_distance_calc_30_10_20(dmin_u, dmin_b, radious);
        
%         H = zeros(N,M);
%         G = zeros(N,M);
%         for n=1:N
%             for m=1:M
%                 H(n,m) = g(n,m)*(dis(n) ^eta);
%                 G(n,m)=(abs(H(n,m)^2)/sigma_m);
%             end
%         end
        
        p1_total=zeros(M,2);
        p2_total=zeros(M,2);
        q_total=zeros(M,2);
        rate = zeros(2);
        s_match_total = zeros(M,2);
        gamma_match= zeros(M,2);
        for b3=1:2
            %bla bla bla
            fprintf('___________________________________________________');
            fprintf('');
            fprintf('\npt = ');
            fprintf('%g ', pt(b));
            fprintf('       b2 = ');
            fprintf('%g ', b2);
            fprintf('       b3 = ');
            fprintf('%g ', b3);
            fprintf('\n');
            
            %channel_assignment for specific pt
            [s_match_total(:,:),gamma_match(:,:),break_const_1] = function_channel_assignment_19_10_20(q_match,BW);
            %if channel assignement-> succesful
            if break_const_1==0
                %power allocation for specific pt
                [p1_total(:,b3),p2_total(:,b3), q_total(:,b3), rate(b3), break_const_2] = function_power_allocation_01_11_20(s_match_total(:,:),gamma_match(:,:),pt(b));
                if break_const_2==0
                    % bla bla bla
                    fprintf('\n s_match_total = ');
                    fprintf('%g ',  s_match_total(:,:));
                    fprintf('\n p1_total = ');
                    fprintf('%g ',  p1_total(:,b3));
                    fprintf('\n p2_total = ');
                    fprintf('%g ',  p2_total(:,b3));
                    fprintf('\n q_total = ');
                    fprintf('%g ',  q_total(:,b3));
                    fprintf('\n');
                    
                    q_match(:) = q_total(:,b3);
                    %optional s_match and q infos -> like total_1_11_20.m line 117-123
                    %optional: while case infos -> like total_1_11_20.m line 136-142
                else
                    fprintf('\nPower allocation is impossible\n');
                    %{
                    b = b-1;
                    fprintf('\nI reduce b so as to try again\n');
                    fprintf('\n now pt = ');
                    fprintf('%g ',pt(b));
                    %}
                    break;   
                end
            else
                fprintf('\nChannel asignement is impossible\n');
                %{
                b = b-1;
                fprintf('\nI reduce b so as to try again\n');
                fprintf('\n now pt = ');
                fprintf('%g ',pt(b));
                %}
                break;
            end
        end
        if break_const_1==0 && break_const_2==0
            SumRate_Repeat(b2,b) = rate(b3);
            p1_final(:,b2,b)=p1_total(:,b3);
            p2_final(:,b2,b)=p2_total(:,b3);
            q_final(:,b2,b)=q_total(:,b3);
            b2 = b2+1;
        else
            fprintf('\n!!! next iteration !!!\n');
            %continue;
        end
    end
    if break_const_1==0 && break_const_2==0
        SumRate(b) = mean(SumRate_Repeat(:,b));
        b = b+1;
    else               %else repeat-try again for the same pt 
       fprintf('\n!!! one more try !!!\n');
       %continue;
    end
    
end



%......plot......
Xaxis = pt;
Yaxis = SumRate/10^6;

plot(Xaxis,Yaxis,':')
title('Jointly optimize channel assignmnet & power allocation')
xlabel('Power of BS with 4 Users (Watt)')
ylabel('Total Sum Rate of System (Mbps)')
