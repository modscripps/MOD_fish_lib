i = 150;


%%
figure('units','inches','position',[0         0  21.1667   11.5556]);
ax(1)=subplot('Position',[.10 .08 .20 .65]);
ax(2)=subplot('Position',[.32 .08 .20 .65]);
ax(3)=subplot('Position',[0.5600    0.434    0.3500    0.3000]);
ax(4)=subplot('Position',[0.5600    0.08    0.3500    0.3000]);
ax(5)=subplot('Position',[0.1000    0.8700    0.8100    0.1000]);
ax(6)=subplot('Position',[0.1000    0.7600    0.8100    0.1000]);

idx_pump=double(Profile.sumPa.a3>5e-7);

a=1;
% Epsilon 1 and 2
semilogx(ax(a),Profile.epsilon_co(:,1),Profile.pr,'k.','linewidth',1)
hold(ax(a),'on')
semilogx(ax(a),Profile.epsilon_co(:,2),Profile.pr,'.','color',[0.5 0.5 0.5],'linewidth',1)
% Highlight fom < 1.15 for 1 and 2
id_good=Profile.epsi_fom(:,1)<1.15;
semilogx(ax(a),Profile.epsilon_co(~id_good,1),Profile.pr(~id_good),'r.','linewidth',1)
id_good=Profile.epsi_fom(:,2)<1.15;
semilogx(ax(a),Profile.epsilon_co(~id_good,2),Profile.pr(~id_good),'.','linewidth',1)
hold(ax(a),'off')
grid(ax(a),'on')
xlabel(ax(a),'\epsilon [W kg^{-1}]')
ylabel(ax(a),'Pr [dbar]')
axis(ax(a),'ij')

a=2;
% Chi 1 and 2
semilogx(ax(a),Profile.chi(:,1),Profile.pr,'k.','linewidth',1)
hold(ax(a),'on')
semilogx(ax(a),Profile.chi(:,2),Profile.pr,'.','color',[0.5 0.5 0.5],'linewidth',1)
% Highlight fom < 10 for 1 and 2
id_good=Profile.chi_fom(:,1)<10; %hacking for now - all chi are good
semilogx(ax(a),Profile.chi(~id_good,1),Profile.pr(~id_good),'r.','linewidth',1)
id_good=Profile.chi_fom(:,2)<10; %hacking for now - all chi are good
semilogx(ax(a),Profile.chi(~id_good,2),Profile.pr(~id_good),'r.','linewidth',1)
hold(ax(a),'off')
grid(ax(a),'on')
xlabel(ax(a),'\chi [˚C^2 s^{-1}]')
ylabel(ax(a),'Pr [dbar]')
axis(ax(a),'ij')

a=3;
grid(ax(a),'on')
xlabel(ax(a),'cpm')
ylabel(ax(a),'\Phi_{shear}')

a=4;
grid(ax(a),'on')
xlabel(ax(a),'cpm')
ylabel(ax(a),'\Phi_{TG}')

a=5;
grid(ax(a),'on')
ax(a).XTickLabel='';

a=6;
grid(ax(a),'on')


    local_pr         = Profile.pr(i);
    local_range      = Profile.ind_range_epsi(i,1):Profile.ind_range_epsi(i,2);
    local_epsi_1     = Profile.epsilon_co(i,1);
    local_epsi_2     = Profile.epsilon_co(i,2);
    local_chi_1      = Profile.chi(i,1);
    local_chi_2      = Profile.chi(i,2);
    local_epsi_fom_1 = Profile.epsi_fom(i,1);
    local_epsi_fom_2 = Profile.epsi_fom(i,2);
    local_chi_fom_1  = Profile.chi_fom(i,1);
    local_chi_fom_2  = Profile.chi_fom(i,2);
    local_kvis       = Profile.kvis(i);
    local_ktemp      = Profile.ktemp(i);
    local_kc_epsi_1  = Profile.sh_kc(i,1);
    local_kc_epsi_2  = Profile.sh_kc(i,2);
    local_kc_chi_1   = Profile.tg_kc(i,1);
    local_kc_chi_2   = Profile.tg_kc(i,2);
    local_kmin       = 3;
    local_k          = Profile.k(i,:);
    local_Pt_1       = Profile.Pt_Tg_k.t1(i,:);
    local_Pt_2       = Profile.Pt_Tg_k.t2(i,:);
    local_Ps_1       = Profile.Ps_shear_co_k.s1(i,:);
    local_Ps_2       = Profile.Ps_shear_co_k.s2(i,:);

    local_epsi_good_1 = local_epsi_fom_1<1.15 & idx_pump(i)==0;
    local_epsi_good_2 = local_epsi_fom_2<1.15 & idx_pump(i)==0;
    local_chi_good_1  = 1:length(local_chi_fom_1); %hacking for now - all chi are good
    local_chi_good_2  = 1:length(local_chi_fom_2); %hacking for now - all chi are good

    kin_epsi_1 = local_k>=local_kmin & local_k<=local_kc_epsi_1;
    kin_epsi_2 = local_k>=local_kmin & local_k<=local_kc_epsi_2;
    kin_chi_1 = local_k>=local_kmin & local_k<=local_kc_chi_1;
    kin_chi_2 = local_k>=local_kmin & local_k<=local_kc_chi_2;

    [kbatch_1,Pbatch_1]=batchelor(local_epsi_1,local_chi_1,local_kvis,local_ktemp);
    nanmask_1=isnan(Pbatch_1);
    [kbatch_2,Pbatch_2]=batchelor(local_epsi_2,local_chi_2,local_kvis,local_ktemp);
    nanmask_2=isnan(Pbatch_2);

    for N=1:2

        switch N
            case 1
                local_epsi = local_epsi_1;
                local_chi = local_chi_1;
                local_epsi_fom = local_epsi_fom_1;
                local_chi_fom = local_chi_fom_1;
                local_kc_epsi = local_kc_epsi_1;
                local_kc_chi = local_kc_chi_1;
                local_Pt = local_Pt_1;
                local_Ps = local_Ps_1;
                local_epsi_good = local_epsi_good_1;
                local_chi_good = local_chi_good_1;
                kin_epsi = kin_epsi_1;
                kin_chi = kin_chi_1;
                local_epsi = local_epsi_1;
                local_chi = local_chi_1;
                kbatch = kbatch_1;
                Pbatch = Pbatch_1;
                nanmask = nanmask_1;

                color_data = 'k';
                color_epsi='g';
                if ~local_epsi_good
                    color_epsi='r';
                end
                color_chi='g';
                if ~local_chi_good
                    color_chi='r';
                end
            case 2
                local_epsi = local_epsi_2;
                local_chi = local_chi_2;
                local_epsi_fom = local_epsi_fom_2;
                local_chi_fom = local_chi_fom_2;
                local_kc_epsi = local_kc_epsi_2;
                local_kc_chi = local_kc_chi_2;
                local_Pt = local_Pt_2;
                local_Ps = local_Ps_2;
                local_epsi_good = local_epsi_good_2;
                local_chi_good = local_chi_good_2;
                kin_epsi = kin_epsi_2;
                kin_chi = kin_chi_2;
                local_epsi = local_epsi_2;
                local_chi = local_chi_2;
                kbatch = kbatch_2;
                Pbatch = Pbatch_2;
                nanmask = nanmask_2;

                color_data = [0.6 0.6 0.6];
                color_epsi=[0.4 1 0.4];
                if ~local_epsi_good
                    color_epsi=[0.4 1 0.4];
                end
                color_chi=[0.4 1 0.4];
                if ~local_chi_good
                    color_chi=[0.4 1 0.4];
                end
        end
        if sum(~nanmask)>2
            local_time     = Profile.epsi.dnum(local_range)-Profile.epsi.dnum(1);
            switch N
                case 1
                    local_s      = Profile.epsi.s1_volt(local_range);
                    local_t       = Profile.epsi.t1_volt(local_range);
                case 2
                    local_s      = Profile.epsi.s2_volt(local_range);
                    local_t       = Profile.epsi.t2_volt(local_range);
            end

            Pbatch=interp1(kbatch(~nanmask),Pbatch(~nanmask),local_k);
            [Pnasm,~] = nasmyth(local_epsi,local_kvis,local_k);

            a=1;
            hold(ax(a),'on')
            sc(a,N)=scatter(ax(a),local_epsi,local_pr,100,color_epsi,'filled','MarkerEdgeColor',color_data);
            hold(ax(a),'off')
            a=2;
            hold(ax(a),'on')
            sc(a,N)=scatter(ax(a),local_chi,local_pr,100,color_chi,'filled','MarkerEdgeColor',color_data);
            hold(ax(a),'off')
            a=3;
            hold(ax(a),'on')
            ll{1}=loglog(ax(a),local_k,local_Ps,'color',color_data,'LineWidth',1);
            ll{2}=loglog(ax(a),local_k(kin_epsi),local_Ps(kin_epsi),'color',color_epsi,'LineWidth',2);
            ll{3}=loglog(ax(a),local_k,Pnasm,'--','color',color_data,'LineWidth',2);
            hold(ax(a),'off')
            a=4;
            hold(ax(a),'on')
            ll{4}=loglog(ax(a),local_k,local_Pt,'color',color_data,'LineWidth',2);
            ll{5}=loglog(ax(a),local_k(kin_chi),local_Pt(kin_chi),'color',color_chi,'LineWidth',2);
            ll{6}=loglog(ax(a),local_k,Pbatch,'--','color',color_data,'LineWidth',2);
            hold(ax(a),'off')

            switch N
                case 1
                    text(ax(3),2,2e-7,...
                        ['\' sprintf('epsilon=%1.2e [W kg^{-1}]\n FOM= %1.2f ',...
                        local_epsi,local_epsi_fom)])
                    text(ax(4),2,2e-5,...
                        ['\' sprintf('chi=%1.2e [˚C^2 s^{-1}]\n FOM= %1.2f ',...
                        local_chi,local_chi_fom)])
                case 2
                    text(ax(3),2,3e-8,...
                        ['\' sprintf('epsilon=%1.2e [W kg^{-1}]\n FOM= %1.2f ',...
                        local_epsi,local_epsi_fom)])
                    text(ax(4),2,3e-6,...
                        ['\' sprintf('chi=%1.2e [˚C^2 s^{-1}]\n FOM= %1.2f ',...
                        local_chi,local_chi_fom)])
            end
            

            set(ax(3:4),'XScale','log','YScale','log')
            set(ax(3),'XLim',[1 1e3],'YLim',[1e-8 1e0])
            set(ax(4),'XLim',[1 1e3],'YLim',[1e-6 1e-1])
            


            a=5;
            ax(a).NextPlot = 'add';
            plot(ax(a),local_time,local_s,'color',color_data);
            grid(ax(a),'on')
            title(ax(5),sprintf("Profile %i, scan %i, Pr %3.2f",profile_id,i,local_pr))


            a=6;
            ax(a).NextPlot = 'add';
            plot(ax(a),local_time,detrend(local_t,'constant'),'color',color_data);
            grid(ax(a),'on')
           
            % set(ax(5),'YLim',[-1e-8 1e-3])
            
            % set(ax(6),'XLim',[1 1e3],'YLim',[1e-6 1e-1])

        end
    end %end for N=1:2

    [ax(:).FontSize] = deal(16);
    ax(5).YLabel.String = 'shear';
    ax(6).YLabel.String = 'fpo7';i = 150;


%%
figure('units','inches','position',[0         0  21.1667   11.5556]);
ax(1)=subplot('Position',[.10 .08 .20 .65]);
ax(2)=subplot('Position',[.32 .08 .20 .65]);
ax(3)=subplot('Position',[0.5600    0.434    0.3500    0.3000]);
ax(4)=subplot('Position',[0.5600    0.08    0.3500    0.3000]);
ax(5)=subplot('Position',[0.1000    0.8700    0.8100    0.1000]);
ax(6)=subplot('Position',[0.1000    0.7600    0.8100    0.1000]);

idx_pump=double(Profile.sumPa.a3>5e-7);

a=1;
% Epsilon 1 and 2
semilogx(ax(a),Profile.epsilon_co(:,1),Profile.pr,'k.','linewidth',1)
hold(ax(a),'on')
semilogx(ax(a),Profile.epsilon_co(:,2),Profile.pr,'.','color',[0.5 0.5 0.5],'linewidth',1)
% Highlight fom < 1.15 for 1 and 2
id_good=Profile.epsi_fom(:,1)<1.15;
semilogx(ax(a),Profile.epsilon_co(~id_good,1),Profile.pr(~id_good),'r.','linewidth',1)
id_good=Profile.epsi_fom(:,2)<1.15;
semilogx(ax(a),Profile.epsilon_co(~id_good,2),Profile.pr(~id_good),'.','linewidth',1)
hold(ax(a),'off')
grid(ax(a),'on')
xlabel(ax(a),'\epsilon [W kg^{-1}]')
ylabel(ax(a),'Pr [dbar]')
axis(ax(a),'ij')

a=2;
% Chi 1 and 2
semilogx(ax(a),Profile.chi(:,1),Profile.pr,'k.','linewidth',1)
hold(ax(a),'on')
semilogx(ax(a),Profile.chi(:,2),Profile.pr,'.','color',[0.5 0.5 0.5],'linewidth',1)
% Highlight fom < 10 for 1 and 2
id_good=Profile.chi_fom(:,1)<10; %hacking for now - all chi are good
semilogx(ax(a),Profile.chi(~id_good,1),Profile.pr(~id_good),'r.','linewidth',1)
id_good=Profile.chi_fom(:,2)<10; %hacking for now - all chi are good
semilogx(ax(a),Profile.chi(~id_good,2),Profile.pr(~id_good),'r.','linewidth',1)
hold(ax(a),'off')
grid(ax(a),'on')
xlabel(ax(a),'\chi [˚C^2 s^{-1}]')
ylabel(ax(a),'Pr [dbar]')
axis(ax(a),'ij')

a=3;
grid(ax(a),'on')
xlabel(ax(a),'cpm')
ylabel(ax(a),'\Phi_{shear}')

a=4;
grid(ax(a),'on')
xlabel(ax(a),'cpm')
ylabel(ax(a),'\Phi_{TG}')

a=5;
grid(ax(a),'on')
ax(a).XTickLabel='';

a=6;
grid(ax(a),'on')


    local_pr         = Profile.pr(i);
    local_range      = Profile.ind_range_epsi(i,1):Profile.ind_range_epsi(i,2);
    local_epsi_1     = Profile.epsilon_co(i,1);
    local_epsi_2     = Profile.epsilon_co(i,2);
    local_chi_1      = Profile.chi(i,1);
    local_chi_2      = Profile.chi(i,2);
    local_epsi_fom_1 = Profile.epsi_fom(i,1);
    local_epsi_fom_2 = Profile.epsi_fom(i,2);
    local_chi_fom_1  = Profile.chi_fom(i,1);
    local_chi_fom_2  = Profile.chi_fom(i,2);
    local_kvis       = Profile.kvis(i);
    local_ktemp      = Profile.ktemp(i);
    local_kc_epsi_1  = Profile.sh_kc(i,1);
    local_kc_epsi_2  = Profile.sh_kc(i,2);
    local_kc_chi_1   = Profile.tg_kc(i,1);
    local_kc_chi_2   = Profile.tg_kc(i,2);
    local_kmin       = 3;
    local_k          = Profile.k(i,:);
    local_Pt_1       = Profile.Pt_Tg_k.t1(i,:);
    local_Pt_2       = Profile.Pt_Tg_k.t2(i,:);
    local_Ps_1       = Profile.Ps_shear_co_k.s1(i,:);
    local_Ps_2       = Profile.Ps_shear_co_k.s2(i,:);

    local_epsi_good_1 = local_epsi_fom_1<1.15 & idx_pump(i)==0;
    local_epsi_good_2 = local_epsi_fom_2<1.15 & idx_pump(i)==0;
    local_chi_good_1  = 1:length(local_chi_fom_1); %hacking for now - all chi are good
    local_chi_good_2  = 1:length(local_chi_fom_2); %hacking for now - all chi are good

    kin_epsi_1 = local_k>=local_kmin & local_k<=local_kc_epsi_1;
    kin_epsi_2 = local_k>=local_kmin & local_k<=local_kc_epsi_2;
    kin_chi_1 = local_k>=local_kmin & local_k<=local_kc_chi_1;
    kin_chi_2 = local_k>=local_kmin & local_k<=local_kc_chi_2;

    [kbatch_1,Pbatch_1]=batchelor(local_epsi_1,local_chi_1,local_kvis,local_ktemp);
    nanmask_1=isnan(Pbatch_1);
    [kbatch_2,Pbatch_2]=batchelor(local_epsi_2,local_chi_2,local_kvis,local_ktemp);
    nanmask_2=isnan(Pbatch_2);

    for N=1:2

        switch N
            case 1
                local_epsi = local_epsi_1;
                local_chi = local_chi_1;
                local_epsi_fom = local_epsi_fom_1;
                local_chi_fom = local_chi_fom_1;
                local_kc_epsi = local_kc_epsi_1;
                local_kc_chi = local_kc_chi_1;
                local_Pt = local_Pt_1;
                local_Ps = local_Ps_1;
                local_epsi_good = local_epsi_good_1;
                local_chi_good = local_chi_good_1;
                kin_epsi = kin_epsi_1;
                kin_chi = kin_chi_1;
                local_epsi = local_epsi_1;
                local_chi = local_chi_1;
                kbatch = kbatch_1;
                Pbatch = Pbatch_1;
                nanmask = nanmask_1;

                color_data = 'k';
                color_epsi='g';
                if ~local_epsi_good
                    color_epsi='r';
                end
                color_chi='g';
                if ~local_chi_good
                    color_chi='r';
                end
            case 2
                local_epsi = local_epsi_2;
                local_chi = local_chi_2;
                local_epsi_fom = local_epsi_fom_2;
                local_chi_fom = local_chi_fom_2;
                local_kc_epsi = local_kc_epsi_2;
                local_kc_chi = local_kc_chi_2;
                local_Pt = local_Pt_2;
                local_Ps = local_Ps_2;
                local_epsi_good = local_epsi_good_2;
                local_chi_good = local_chi_good_2;
                kin_epsi = kin_epsi_2;
                kin_chi = kin_chi_2;
                local_epsi = local_epsi_2;
                local_chi = local_chi_2;
                kbatch = kbatch_2;
                Pbatch = Pbatch_2;
                nanmask = nanmask_2;

                color_data = [0.6 0.6 0.6];
                color_epsi=[0.4 1 0.4];
                if ~local_epsi_good
                    color_epsi=[0.4 1 0.4];
                end
                color_chi=[0.4 1 0.4];
                if ~local_chi_good
                    color_chi=[0.4 1 0.4];
                end
        end
        if sum(~nanmask)>2
            local_time     = Profile.epsi.dnum(local_range)-Profile.epsi.dnum(1);
            switch N
                case 1
                    local_s      = Profile.epsi.s1_volt(local_range);
                    local_t       = Profile.epsi.t1_volt(local_range);
                case 2
                    local_s      = Profile.epsi.s2_volt(local_range);
                    local_t       = Profile.epsi.t2_volt(local_range);
            end

            Pbatch=interp1(kbatch(~nanmask),Pbatch(~nanmask),local_k);
            [Pnasm,~] = nasmyth(local_epsi,local_kvis,local_k);

            a=1;
            hold(ax(a),'on')
            sc(a,N)=scatter(ax(a),local_epsi,local_pr,100,color_epsi,'filled','MarkerEdgeColor',color_data);
            hold(ax(a),'off')
            a=2;
            hold(ax(a),'on')
            sc(a,N)=scatter(ax(a),local_chi,local_pr,100,color_chi,'filled','MarkerEdgeColor',color_data);
            hold(ax(a),'off')
            a=3;
            hold(ax(a),'on')
            ll{1}=loglog(ax(a),local_k,local_Ps,'color',color_data,'LineWidth',1);
            ll{2}=loglog(ax(a),local_k(kin_epsi),local_Ps(kin_epsi),'color',color_epsi,'LineWidth',2);
            ll{3}=loglog(ax(a),local_k,Pnasm,'--','color',color_data,'LineWidth',2);
            hold(ax(a),'off')
            a=4;
            hold(ax(a),'on')
            ll{4}=loglog(ax(a),local_k,local_Pt,'color',color_data,'LineWidth',2);
            ll{5}=loglog(ax(a),local_k(kin_chi),local_Pt(kin_chi),'color',color_chi,'LineWidth',2);
            ll{6}=loglog(ax(a),local_k,Pbatch,'--','color',color_data,'LineWidth',2);
            hold(ax(a),'off')

            switch N
                case 1
                    text(ax(3),2,2e-7,...
                        ['\' sprintf('epsilon=%1.2e [W kg^{-1}]\n FOM= %1.2f ',...
                        local_epsi,local_epsi_fom)])
                    text(ax(4),2,2e-5,...
                        ['\' sprintf('chi=%1.2e [˚C^2 s^{-1}]\n FOM= %1.2f ',...
                        local_chi,local_chi_fom)])
                case 2
                    text(ax(3),2,3e-8,...
                        ['\' sprintf('epsilon=%1.2e [W kg^{-1}]\n FOM= %1.2f ',...
                        local_epsi,local_epsi_fom)])
                    text(ax(4),2,3e-6,...
                        ['\' sprintf('chi=%1.2e [˚C^2 s^{-1}]\n FOM= %1.2f ',...
                        local_chi,local_chi_fom)])
            end
            

            set(ax(3:4),'XScale','log','YScale','log')
            set(ax(3),'XLim',[1 1e3],'YLim',[1e-8 1e0])
            set(ax(4),'XLim',[1 1e3],'YLim',[1e-6 1e-1])
            


            a=5;
            ax(a).NextPlot = 'add';
            plot(ax(a),local_time,local_s,'color',color_data);
            grid(ax(a),'on')
            title(ax(5),sprintf("Profile %i, scan %i, Pr %3.2f",profile_id,i,local_pr))


            a=6;
            ax(a).NextPlot = 'add';
            plot(ax(a),local_time,detrend(local_t,'constant'),'color',color_data);
            grid(ax(a),'on')
           
            % set(ax(5),'YLim',[-1e-8 1e-3])
            
            % set(ax(6),'XLim',[1 1e3],'YLim',[1e-6 1e-1])

        end
    end %end for N=1:2

    [ax(:).FontSize] = deal(16);
    ax(5).YLabel.String = 'shear';
    ax(6).YLabel.String = 'fpo7';