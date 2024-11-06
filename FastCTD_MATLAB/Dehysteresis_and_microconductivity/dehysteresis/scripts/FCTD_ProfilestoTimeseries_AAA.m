function FCTD_corrected_timeseries = FCTD_ProfilestoTimeseries_AAA(FCTD,datapath)
% This function is for converting FCTD profile data back into timeseries
% after having done the processing code
%
% Alex Andriatis
% 2021-07-21

FCTD_corrected_timeseries=[];

directions={'up','down'};
variables=fieldnames(FCTD.up);
nprof=length(FCTD.up.time);

T=0;
for k=1:nprof
    T = T+length(FCTD.up.time{k});
end

for n=1:length(variables)
    variable=variables{n};
    FCTD_corrected_timeseries.(variable)=NaN(T,size(FCTD.up.(variable){1},2));
end

for i=1:length(directions)
    direction=directions{i};
    for n=1:length(variables)
        variable=variables{n};
        nrec = 1;
        for k=1:nprof
            tmp=FCTD.(direction).(variable){k};
            nk = size(tmp,1);
            nvar = size(tmp,2);
            FCTD_corrected_timeseries.(variable)(nrec:nrec+nk-1,1:nvar)=FCTD.(direction).(variable){k};
            nrec = nrec+nk;
        end
    end
end


headvars=fieldnames(FCTD.header);

for n=1:length(headvars)
    headvar=headvars{n};
    FCTD_corrected_timeseries.(headvar)=NaN(T,1);
    for i=1:length(directions)
        direction=directions{i};
        nrec = 1;
        for k=1:nprof
            tmp = FCTD.header.(headvar)(k)*ones(size(FCTD.(direction).time{k}));
            nk = size(tmp,1);
            nvar = size(tmp,2);            
            FCTD_corrected_timeseries.(headvar)(nrec:nrec+nk-1,1:nvar)=tmp;
            nrec = nrec+nk;
        end
    end
end


FCTD_corrected_timeseries.profile=NaN(T,1);
for i=1:length(directions)
    direction = directions{i};
    nrec = 1;
    for k=1:nprof
        tmp = k*ones(size(FCTD.(direction).time{k}));
        nk = size(tmp,1);
        nvar = size(tmp,2);             
        FCTD_corrected_timeseries.profile(nrec:nrec+nk-1,1:nvar)=tmp;
        nrec = nrec+nk;
    end
end


[~,I]=sort(FCTD_corrected_timeseries.time);
variables=fieldnames(FCTD_corrected_timeseries);
for n=1:length(variables)
    variable=variables{n};
    FCTD_corrected_timeseries.(variable)=FCTD_corrected_timeseries.(variable)(I,:);
end

save(datapath,'-struct','FCTD_corrected_timeseries','-v7.3');
disp(['Saved converted FCTD profiles as timeseries in ' datapath]);
