function [h, hl]=logscatterlines(ax,xli)
% FUntion to add 1:1, 2:1, 5:1 and 10:1 lines to scatter plots
%   Useful when comparing methods/ data that should fall on a 1:1 line
% Required input:
%   ax: axes handles to add lines
% Outputs:
%   h: Line handles to change their appearance
%   hl: legend handle
% Optional input:
%   xli: axes limits, both will be the same. By default it will use the
%   largest range of the xli or the yli to add the lines

if nargin<2
    tmpx=get(ax,'xlim');
    tmpy=get(ax,'ylim');
    xli(1)=min([tmpx(1) tmpy(1)]);
    xli(2)=max([tmpx(1) tmpy(2)]);
end

axes(ax)
hold on;

h(1)=loglog(xli,xli,'k-');
h(2)=loglog(xli,xli*2,'k--');
h(3)=loglog(xli,xli/2,'k--');

h(4)=loglog(xli,xli*5,'k-.');
h(5)=loglog(xli,xli/5,'k-.');

h(6)=loglog(xli,xli*10,'k:');
h(7)=loglog(xli,xli/10,'k:');

set(ax,'xlim',xli,'ylim',xli,'xscale','log','yscale','log');
hl=legend(h([1;2;4;6]),'1:1','2:1','5:1','10:1');
set(hl,'box','off','fontsize',8,'location','southeast');
box on;