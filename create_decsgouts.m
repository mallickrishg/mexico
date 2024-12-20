function g = create_decsgouts(filename)

% data = readtable(filename,"FileType","text","NumHeaderLines",1,"ReadVariableNames",false);
% 
% npts = length(data{:,1});
% Rdata = [2,npts,data{:,1}',data{:,2}']';
% 
% g = decsg(Rdata);

data = readtable(filename);

x1 = data{:,1};
y1 = data{:,2};
x2 = data{:,3};
y2 = data{:,4};
l1 = data{:,5};
l2 = data{:,6};
npts = length(x1);

g = [ones(1,npts).*2;...
    x1';x2';y1';y2';...
    l1';...
    l2'];

end