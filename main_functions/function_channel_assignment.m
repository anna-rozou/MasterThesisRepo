%..........................Channel Assignment...................

function [S_match,gamma_match, break_const] = function_channel_assignment_19_10_20(q,BW)
global M W G

N =2*M;
Bc = BW/M;

%users_PL = zeros(N,M);
[~,index] = sort(G,2,'descend');
users_PL = index;

%initialize S_match & S_unmatch
S_match = zeros(M,2);
gamma_match = zeros(M,2);
S_unmatch = 1:N;
s_match_counter = zeros(M,1);   %shows which sub channel to use: 1 or 2
while_counter = 0;

%combinations for 3 users
combos = [1 2 ; 1 3 ; 2 3];
break_const = 0;
while any(S_unmatch,'all')
    while_counter = while_counter +1; 
    for n=1:N
        %bla bla bla
        %{
        fprintf('.......................................................')
        X5 = sprintf('\ncounter=%d ',while_counter);
        disp(X5)
        X8 = sprintf('user=%d \n',n);
        disp(X8)
        %}
        prefered_channel = users_PL(n,1);

        %ismember(n,S_match)=1 if n exist on S_match
        if ismember(n,S_match)~=1              
            if s_match_counter(prefered_channel)==0
                s_match_counter(prefered_channel) = s_match_counter(prefered_channel) +1;
                S_match(prefered_channel,s_match_counter(prefered_channel)) = n;
                S_unmatch(n) = 0;

                gamma_match(prefered_channel,s_match_counter(prefered_channel))=G(n,prefered_channel);
                
            elseif s_match_counter(prefered_channel)==1 && S_match(prefered_channel,1)~=n
                %...........................................g1>g2................
                s_match_counter(prefered_channel) = s_match_counter(prefered_channel) +1;
                
                G1 = G(S_match(prefered_channel,1),prefered_channel); %old user's gamma of this channel
                G2 = G(n,prefered_channel);                           %the new user's gamma  
                
                if G2 > G1 * W(prefered_channel,2)/W(prefered_channel,1)
                    temp = S_match(prefered_channel,1);
                    S_match(prefered_channel,1)  = n;
                    S_match(prefered_channel,s_match_counter(prefered_channel)) = temp;
                    %fprintf('\n1-> G2>G1*W2/W1');
                    gamma_match(prefered_channel,1)=G2;
                    gamma_match(prefered_channel,2)=G1;
                    S_unmatch(n) = 0;
                elseif G1 > G2 * W(prefered_channel,2)/W(prefered_channel,1)
                    S_match(prefered_channel, s_match_counter(prefered_channel)) = n;
                    gamma_match(prefered_channel,2)=G2;
                    S_unmatch(n) = 0;
                    %fprintf('\n1-> G1>G2*W2/W1');
                else
                    if G2>G1
                        old_user = S_match(prefered_channel,1);
                        S_match(prefered_channel,1)  = n;
                        gamma_match(prefered_channel,1)=G2;
                        S_unmatch(n) = 0;
                        S_unmatch(old_user) = old_user;
                        %fprintf('\n1-> G2>G1');
                    end
                        s_match_counter(prefered_channel) = s_match_counter(prefered_channel) -1;
                        %fprintf('\n1-> Did not follow the constraint');
                end
                
            
            elseif s_match_counter(prefered_channel) == 2
                wrong_user=0;
                
                user1 = S_match(prefered_channel,1);
                user2 = S_match(prefered_channel,2);
                this_users = [user1, user2, n];
                
                R = zeros(3,3);
                G_ca = zeros(3,1);     %users Gamma,3 rows:3 users 2 coloums for different Weights
                p = zeros(3,2);        %power constraints
                omega= zeros(3,1);
                %BLA BLA BLA
                %{
                fprintf('')
                fprintf('\nThe this_users is: [');
                fprintf('%g ', this_users);
                fprintf(']\n');
                %}
                
                G_ca(1) = G(this_users(1),prefered_channel);
                G_ca(2) = G(this_users(2),prefered_channel);
                G_ca(3) = G(this_users(3),prefered_channel);
                
                const = zeros(3,1);
                for k3=1:3
                    if G_ca(combos(k3,2)) > G_ca(combos(k3,1))* W(prefered_channel,2)/W(prefered_channel,1)
                         temp1 = combos(k3,1);
                         combos(k3,1) = combos(k3,2);
                         combos(k3,2) = temp1;
                         o1 = W(prefered_channel,2)*G_ca(combos(k3,2)) - W(prefered_channel,1)*G_ca(combos(k3,1));
                         o2 = G_ca(combos(k3,1))*G_ca(combos(k3,2)) * (W(prefered_channel,1)-W(prefered_channel,2));
                         omega(k3)= o1/o2;
                         
                         if q(prefered_channel) < 2*omega(k3)
                             q(prefered_channel) = 2*omega(k3);
                         end
                         p(k3,1) = omega(k3);
                         p(k3,2) = q(prefered_channel)-omega(k3);
                         %fprintf('\n2-> G2>G1*W2/W1');
                    elseif G_ca(combos(k3,1)) > G_ca(combos(k3,2))* W(prefered_channel,2)/W(prefered_channel,1)
                        o1 = W(prefered_channel,2)*G_ca(combos(k3,2)) - W(prefered_channel,1)*G_ca(combos(k3,1));
                        o2 = G_ca(combos(k3,1))*G_ca(combos(k3,2)) * (W(prefered_channel,1)-W(prefered_channel,2));
                        omega(k3)= o1/o2;
                        if q(prefered_channel) < 2*omega(k3)
                            q(prefered_channel) = 2*omega(k3);
                        end
                        p(k3,1) = omega(k3);
                        p(k3,2) = q(prefered_channel)-omega(k3);
                         %fprintf('\n2-> G1>G2*W2/W1');
                     else
                         %........delete this combo..........
                         const(k3) = 1;
                         %fprintf('\n2-> Did not follow the constraint');
                     end 
                end
                omega
                p
                q                
                %sum rate calculation
                for k3=1:3

                    
                    if const(k3) == 0
                        rate1 = W(prefered_channel,1)*Bc*log(1+omega(k3)*G_ca(combos(k3,1)));
                        rate2 = W(prefered_channel,2)*Bc*log((q(prefered_channel)*G_ca(combos(k3,2))+1) / (omega(k3)*G_ca(combos(k3,2))+1));
                        R(combos(k3,1),combos(k3,2)) = rate1+rate2;
                    else
                        R(combos(k3,1),combos(k3,2)) = -1;
                    end
                end
                
                maximum_rate = max(R(:));
                [max1 ,max2] = find(R==maximum_rate);
                
                max_user1 = this_users(max1);
                max_user2 = this_users(max2);
                
                for k3=1:3
                    if (this_users(k3)~=max_user1) && (this_users(k3)~=max_user2)
                        wrong_user = this_users(k3);
                    end
                end
                
                S_unmatch(wrong_user) = wrong_user;
                %remove channel from wrong_user PL
                for m2=1:M
                    if users_PL(wrong_user,m2)==prefered_channel
                        for m3=m2:M-1
                            users_PL(wrong_user,m3) = users_PL(wrong_user,m3+1);
                        end
                        users_PL(wrong_user,M)=0;
                    end
                end
                if G_ca(max1)>G_ca(max2)
                    S_match(prefered_channel,1) = max_user1;
                    S_match(prefered_channel,2) = max_user2;
                    gamma_match(prefered_channel,1) = G_ca(max1);
                    gamma_match(prefered_channel,2) = G_ca(max2);
                else
                    S_match(prefered_channel,1) = max_user2;
                    S_match(prefered_channel,2) = max_user1;
                    gamma_match(prefered_channel,1) = G_ca(max2);
                    gamma_match(prefered_channel,2) = G_ca(max1);
                end
                S_unmatch(max_user1) = 0;
                S_unmatch(max_user2) = 0;
                
            end
        else
            %fprintf('This user has already been allocated')
        end
        %bla bla bla
        %{
        %display things
        S_match
        users_PL
        fprintf('')
        fprintf('\nThe S_unmatch is: [');
        fprintf('%g ', S_unmatch);
        fprintf(']\n');
        %}
    end
    if while_counter >10
        S_match = zeros(M,2);
        gamma_match = zeros(M,2);
        break_const = 1;
        fprintf('STOP');
        break;
    end
end

%................CHECK...................
if break_const == 0
    check1 = 0;
    for m=1:M
        if gamma_match(m,1) < gamma_match(m,2) * W(prefered_channel,2)/W(prefered_channel,1)
            check1 = check1 +1;
        end
    end
    if check1 == 0
        fprintf('\nCORRECT ALLOCATION 1 !!!')
    else
        fprintf('\nWRONG ALLOCATION 1 !!!')
    end
    
    check2 = 0;
    for m=1:M
        if gamma_match(m,1) < gamma_match(m,2) 
            check2 = check2 +1;
        end
    end
    
    if check2 == 0
        fprintf('\nCORRECT ALLOCATION 2 !!!')
    else
        fprintf('\nWRONG ALLOCATION 2 !!!')
    end
end

S_match
gamma_match
break_const

end
