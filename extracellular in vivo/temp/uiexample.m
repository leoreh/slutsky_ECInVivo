function uiexample
f = figure('Visible','on');
ax = axes('Units','pixels');
plot(rand(10000, 1))
axpos = get(ax,'Position');
xylim = [ax.XLim; ax.YLim]';
XYlimitpos = [axpos(1)-3 axpos(2)-25 20 20;...
    axpos(1)+ax.Position(3)-10 axpos(2)-25 20 20;...
    axpos(1)-25 axpos(2)-8 20 20;...
    axpos(1)-25 axpos(2)+axpos(4)-8 20 20];
for ind = 1:4
    xyedit(ind) = uicontrol('style','edit','Position',XYlimitpos(ind,:),...
        'string',xylim(ind),'callback',{@updateplot,ax,ind});
end
end
function updateplot(hObject,event,axhandle,index)
XYLim = [axhandle.XLim;axhandle.YLim]';
XYLim(index)=str2num(hObject.String);
axhandle.XLim = XYLim(:,1)';
axhandle.YLim= XYLim(:,2)';
end