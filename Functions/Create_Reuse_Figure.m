function [curFig] = Create_Reuse_Figure(NewFig,FigName,figpos,varargin)
%curFig = Create_Reuse_Figure(1_0_[],'Sample Figure',[100 100 300 400]);
% Create a new figure or use the existing one

clf_flag = 1;
if any(strcmpi(varargin,'NoCLF'))
    clf_flag = 0;   
end
    

curFigall = findobj('type','figure','name',FigName);

if ~isempty(NewFig) && NewFig == 1
	curFig = figure('Name',FigName,'Position',figpos);
	%curFig.Position = figpos;
else
	if ~isempty(curFigall)
		curFig = curFigall(1); 
        if clf_flag, clf(curFig); end
        figure(curFig);
		if curFig.Position(4) ~= figpos(4)
			curFig.Position(3:4) = figpos(3:4);
		end
	else
		curFig = figure('Name',FigName,'Position',figpos);
	end
end
end

