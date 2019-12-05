function W = mygraph(coord)
numOfnode = size(coord,1);
param.type = 'knn';
param.k = 6;
G = gsp_nn_graph(coord, param);
G.coords = coord;
for i = 1:numOfnode
    for j = 1:numOfnode
        if G.W(i,j)~=0 && G.W(i,j) >= 0.3865
            G.W(i,j) = 1;
        elseif G.W(i,j)~=0 && G.W(i,j) < 0.3865
            G.W(i,j) = 1/(sqrt(8)/2);
        end
    end
end

W = zeros(numOfnode,numOfnode);
for i = 1:numOfnode
    for j = 1:numOfnode
        W(i,j) = G.W(i,j);
        if abs(i-j) == 2 || abs(i-j) == 15 || abs(i-j) == 16
            W(i,j) = 0;
        end
    end
end
W(1,11) = 0; W(11,1) = 0;
W(8,14) = 0; W(14,8) = 0;
W(42,57) = 0; W(57,42) = 0;
W(47,64) = 0; W(64,47) = 0;