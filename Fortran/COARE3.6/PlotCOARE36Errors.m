clear all; close all; 

% Set this to true to output plots to png files.
output_plots = false;

fid = fopen('test_36_output_matlab_072821.txt','r');
fgetl(fid);
a = reshape(fscanf(fid,'%g'), [21 116]);

fid = fopen('test_36_output_f90coolskin.txt','r');
fgetl(fid)
b = fscanf(fid,'%g'); disp(size(b))
b = reshape(b, [21 116]);
b(b==-999999) = NaN;
b(7,:) = 1e3*b(7,:); % qsr --> g/kg

abserr = abs((b-a));
relerr = abs((b-a)./a);

vv = {'usr','tau','hsb','hlb','hlwebb','tsr','qsr','zot','zoq','Cd','Ch','Ce','L','zet','dter','dqer','tkt','RF','Cdn10','Chn10','Cen10'};
Nvar = length(vv);
gg = {[1 2],[3 4],[6 7],[8 9],[10 11 12],[13 14],[15 17],[19 20 21]}

for n = 1:length(gg)
  varind = gg{n}
  for ii = 1:4
    expind = [1:29]+29*(ii-1);
    close all
    figure; clf
    subplot(211); plot(expind,abserr(varind,expind),'-o'); 
    legend(vv{varind}); ylabel('|error|'); 
    title('coare36flux\_coolskin vs coare36vn_zrf.m')
    subplot(212); plot(expind,relerr(varind,expind),'-o'); 
    ylabel('|rel error|'); 
    xlabel('experiment #')
    set(gca,'XMinorTick','on')

    fname = sprintf('Figures/COARE36CoolSkinCheckAgainstMatlabv36_%.2d_Exp%.3d-%.3d_%s.png', ...
                    n,expind(1),expind(end),datestr(today,['yyyy-mm-' ...
                        'dd']));
    if output_plots; eval(['print -dpng -r400 ' fname]); end
  end
end