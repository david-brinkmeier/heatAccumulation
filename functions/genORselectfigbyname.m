function currentfig = genORselectfigbyname(varargin)
% genORselectfigbyname(figname,identifier)
% two inputs with identifier if plot should persist 

% select figure by specifying its name
% or if it doesn't exist, generate it upon first call

% example: generate figure with a name
% figure('name','banana')
% select this figure for plotting
% selectfigbyname('banana')

if length(varargin) == 2
   % uniqueID = char(java.util.UUID.randomUUID.toString);
   figname = strcat(varargin{1},32,'-',32,varargin{2});
else
    figname = varargin{1};
end

currentfig = findobj('type','figure','name',figname);
if ~isempty(currentfig)
figure(currentfig)
else
    warning on
    disp(strcat('generating figure',32,'"',figname,'"'));
    currentfig = figure('name',figname,'color','w');
end

end

