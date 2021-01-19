%.....................Distance Calculation.................................

function dis = function_distance_calc_30_10_20(dmin_u, dmin_b, radious)
global M
N=2*M;
dmin_u=dmin_u^2;
dmin_b=dmin_b^2;
X = zeros(N, 1);
Y = zeros(N, 1);
R = zeros(N, 1);
x1 = 0;            %center of the circle
y1 = 0;            %center of the circle
nValid=0;
iLoop= 1; 

while nValid < N && iLoop < 1e6
    a=2*pi*rand;
    r=sqrt(rand);
    newX=(radious*r)*cos(a)+x1;
    newY=(radious*r)*sin(a)+y1;
    
    if all(((X(1:nValid) - newX).^2 + (Y(1:nValid) - newY).^2) > dmin_u)
        if (newX^2+newY^2)>dmin_b
            % Success: The new point does not touch existing points:
            nValid    = nValid + 1;  % Append this point
            X(nValid) = newX;
            Y(nValid) = newY;
            R(nValid) = sqrt((newX^2+newY^2));
        end
    end
  iLoop = iLoop + 1;
end

if nValid < N
  error('Cannot find wanted number of points in %d iterations.', iLoop)
end
    
dis = R;
%scatter(X,Y,'filled')

end