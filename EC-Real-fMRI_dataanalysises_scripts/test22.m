% test 1 so com posições
figure()
s = [1 1 1 1 1 2 2 7 7 9 3 3 1 4 10 8 4 5 6 8];
t = [2 3 4 5 7 6 7 5 9 6 6 10 10 10 11 11 8 8 11 9];
weights = [1 1 1 1 3 3 2 4 1 6 2 8 8 9 3 2 10 12 15 16];
G = graph(s,t,weights)

x = [0 0.5 -0.5 -0.5 0.5 0 1.5 0 2 -1.5 -2];
y = [0 0.5 0.5 -0.5 -0.5 2 0 -2 0 0 0];
z = [5 3 3 3 3 0 1 0 0 1 0];
plot(G,'XData',x,'YData',y,'ZData',z,'EdgeLabel',G.Edges.Weight)

% edge line by eight
G.Edges.Weight = randi([10 250],130,1);
G.Edges.LWidths = 7*G.Edges.Weight/max(G.Edges.Weight);
p.LineWidth = G.Edges.LWidths;


%test 2 com setinhas
figure()
G = digraph(1,2:5);
G = addedge(G,2,6:15);
G = addedge(G,15,16:20)
plot(G,'Layout','force')

figure ()
s = [1 1 1 2 2 3 3 4 5 5 6 7];
t = [2 4 8 3 7 4 6 5 6 8 7 8];
H = digraph(s,t);
plot(H,'Layout','force')

figure()
s = [1 2 1 2 2 3 4 5 6 7];
t = [2 1 8 3 7 6 5 6 7 8];
weights = [10 20 10 1 10 1 10 1 12 12];
names = {'A' 'B' 'C' 'D' 'E' 'F' 'G' 'H'};
G = digraph(s,t,weights,names)

x = [0 0.5 -0.5 -0.5 0.5 0 1.5 0 2 -1.5];
y = [0 0.5 0.5 -0.5 -0.5 2 0 -2 0 0];
plot(G,'Layout','force','EdgeLabel',G.Edges.Weight)
% plot(G,'XData',x,'YData',y,'Layout','force','EdgeLabel',G.Edges.Weight)