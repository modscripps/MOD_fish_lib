function PressureTimeseries = get_profile_from_rbr(PressureTimeseries)
% similar to ctd except the I only get upcast
% requires at least
% PressureTimeseries.meta
% PressureTimeseries.P
% PressureTimeseries.dnum
% PressureTimeseries.meta.smfilt=30;
% PressureTimeseries.meta.cutspeed=.1;

% smooth P times series the 30*60 in smoothdata should be a parameter
smdepth=smoothdata(PressureTimeseries.P,'movmean',PressureTimeseries.meta.smfilt); % 30*60 s
% compute speed
w=diff(smdepth)./diff(PressureTimeseries.dnum)/86400;
% get rid of spikes
 w=filloutliers(w,'linear');
% smooth speed just in case we still have high frequency features
smw=smoothdata(w,'movmedian',ceil(100*PressureTimeseries.meta.smfilt/3));
% define upcast  and dowcast as section of speed above and below a threshold
% here .1 m/s The Threshold should be a parameter
downcast(smw>PressureTimeseries.meta.cutspeed)=1; % down
upcast(smw<-PressureTimeseries.meta.cutspeed)=1; % up
% find index where upcast==1 i.e. select only the upcast
idx_upcast=find(upcast==1);
% difference between the upcast index so we can see the jump in ondex
% between cast
didx_upcast=diff(idx_upcast);
% identify and nan sections where indexes are NOT contiguous (i.e. didx_upcast>1)
didx_upcast(didx_upcast>1)=nan;
% I only keep the sections where indexes are contiguous (continuous cast)
idx_upcast=idx_upcast(~isnan(didx_upcast));

% % same procedure for downcast
idx_downcast=find(downcast==1);
didx_downcast=diff(idx_downcast);
didx_downcast(didx_downcast>1)=nan;
idx_downcast=idx_downcast(~isnan(didx_downcast));

% Now I identifiy all the downcast and give a ID number. 
section_number=PressureTimeseries.P.*0;
down_nb=0;
disp("start downcast")
for j=2:length(idx_downcast)
    if diff(idx_downcast([j-1 j]))>1 % this happens when we change cast
        down_nb=down_nb+1;
    else
        section_number(idx_downcast(j))=down_nb; %give all the element of the downcast the same ID number
    end
    if mod(length(idx_downcast)-j,100)==0
    fprintf(".")
    end
    if mod(length(idx_downcast)-j,1000)==0
    fprintf("%i\r\n",length(idx_downcast)-j)
    end
end
disp("done downcast")



disp("start upcast")
up_nb=0;
for j=2:length(idx_upcast)
    if diff(idx_upcast([j-1 j]))>1
        up_nb=up_nb-1;
    else
        section_number(idx_upcast(j))=up_nb;
    end
    if mod(length(idx_upcast)-j,100)==0
    fprintf(".")
    end
    if mod(length(idx_upcast)-j,1000)==0
        fprintf("%i\r\n",length(idx_upcast)-j)
    end

end
disp("done upcast")

%% cleaning /removing short profiles
%downcast
id_section=1;
for j=1:max(section_number)
    P0=PressureTimeseries.P(section_number==id_section);
    if  isempty(P0)
        section_number(section_number==id_section)=0;
        section_number(section_number>id_section)=section_number(section_number>id_section)-1;
        id_section=id_section-1;
    else
        if  ((P0(1)<5) || (diff(P0([1 end]))<20)) % if the profile is short <20m or too shallow 5m
            section_number(section_number==id_section)=0;
            section_number(section_number>id_section)=section_number(section_number>id_section)-1;
            id_section=id_section-1;
        end
    end
    id_section=id_section+1;
    if mod(j,100)==0
        fprintf('% i Profile clean, still %i \r\n',j,max(section_number)-j)
    end
end

%upcast
id_section=-1;
for j=-1:-1:min(section_number)
    P0=PressureTimeseries.P(section_number==id_section);
    if  isempty(P0)
        section_number(section_number==id_section)=0;
        section_number(section_number<id_section)=section_number(section_number<id_section)+1;
        id_section=id_section+1;

    else
        if  ((P0(1)<5) || (diff(P0([end 1]))<20)) % if the profile is short <20m or too shallow 5m
            section_number(section_number==id_section)=0;
            section_number(section_number<id_section)=section_number(section_number<id_section)+1;
            id_section=id_section+1;
        end
    end
    id_section=id_section-1;
    if mod(j,100)==0
        fprintf('% i Profiles cleaned, still %i Profiles. \r\n',-j,abs(min(section_number))+j)
    end

end
PressureTimeseries.section_number=section_number;
nb_upcast=0;
for i=-1:-1:min(PressureTimeseries.section_number)
    nb_upcast=nb_upcast+1;
    idx_upcast=find(PressureTimeseries.section_number==i);
    PressureTimeseries.startprof(nb_upcast)=idx_upcast(1);
    PressureTimeseries.endprof(nb_upcast)=idx_upcast(end);
end


end