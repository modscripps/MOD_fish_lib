function ax=QuickEpsiProfilePlotMHA_accel(Profile,params)
    %function ax=QuickEpsiProfilePlotMHA_accel(Profile,params)
    %MHA 8/2024.
    %As QuickEpsiProfilePlotMHA but with accel. Do a 6-panel plot of eps, chi, w, T and S for basic profile health.
    %basic plotting params can be specified by passing in a structure of the following form:
    %    params.epslims=[1e-10 1e-7];
    %    params.chilims=[1e-10 1e-6];
    %    params.wlims=[.5 .8];
    %    params.alims=[-.05 .05];
    
    if nargin < 2
        params.epslims=[1e-10 1e-7];
        params.chilims=[1e-10 1e-6];
        params.wlims=[.5 .8];
        params.alims=[-.05 .05];
        
    end
    
    Profile.Meta_Data.plot_properties = epsiSetup_set_plot_properties;
    cols = Profile.Meta_Data.plot_properties.Colors;
    
    %params.ylims
    clf
    ax=MySubplot(.1,0.02,0.02,.1,.15,0,6,1);
    
    % Epsilon
    axes(ax(1))
    plt.eps = semilogx(Profile.epsilon_co(:,1),Profile.z,...
        Profile.epsilon_co(:,2),Profile.z);axis ij
    plt.eps(1).Color = cols.s1;
    plt.eps(2).Color = cols.s2;
    xlim(params.epslims)
    grid
    legend('s1','s2','location','southeast')
    xlabel('\epsilon (W/kg)')
    ylabel('depth (m)')
    
    % Chi
    axes(ax(2))
    plt.chi = semilogx(Profile.chi,Profile.z);axis ij
    plt.chi(1).Color = cols.t1;
    plt.chi(2).Color = cols.t2;
    xlim(params.chilims)
    grid
    legend('t1','t2','location','southeast')
    ax(2).YTickLabel='';
    xlabel('\chi (K^2/s)')
    
    % w
    axes(ax(3))
    plot(Profile.w,Profile.z,'color',cols.dPdt);axis ij
    xlim([0 1])
    xlim(params.wlims)
    grid
    % title([Profile.deployment ', #' num2str(Profile.profNum) ', ' datestr(nanmin(Profile.dnum))],'Interpreter', 'none')
    title(['#' num2str(Profile.profNum) ', ' datestr(nanmin(Profile.dnum))],'Interpreter', 'none')
    %title([Profile.deployment ', #' num2str(Profile.profNum)],'Interpreter', 'none')
    ax(3).YTickLabel='';
    xlabel('w (m/s)')
    
    % Acceleration
    axes(ax(4))
    %Put accel onto the Profile grid
    a1=interp1(Profile.epsi.dnum,Profile.epsi.a1_g,Profile.dnum);
    a2=interp1(Profile.epsi.dnum,Profile.epsi.a2_g,Profile.dnum);
    a3=interp1(Profile.epsi.dnum,Profile.epsi.a3_g,Profile.dnum);
    a1m=nanmean(a1);
    a2m=nanmean(a2);
    a3m=nanmean(a3);
    
    plt.acc = plot(a1-a1m,Profile.z,a2-a2m,Profile.z,a3-a3m,Profile.z);axis ij
    plt.acc(1).Color = cols.a1;
    plt.acc(2).Color = cols.a2;
    plt.acc(3).Color = cols.a3;
    lg=legend(['a1-' num2str(a1m,1)],['a2-' num2str(a2m)],['a3-' num2str(a3m)],'location','southeast');
    pos=get(lg,'position');pos(2)=pos(2)+.25;set(lg,'position',pos);
    
    ax(4).YTickLabel='';
    xlim(params.alims)
    grid
    xlabel('acceleration (g)')
    
    % T
    axes(ax(5))
    plot(Profile.ctd.T,Profile.ctd.z,'color',cols.T);axis ij
    %xlim([0 1])
    grid
    ax(5).YTickLabel='';
    xlabel('T (^oC)')
    
    % S
    axes(ax(6))
    plot(Profile.ctd.S,Profile.ctd.z,'color',cols.S);axis ij
    %xlim([0 1])
    grid
    ax(6).YTickLabel='';
    xlabel('S (psu)')
    linkaxes(ax,'y');
    
    
    
    shg
    setpp(15,6)
    