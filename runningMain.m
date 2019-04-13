function M = runningMain()
M = 1;
%N = [3 5 7 10 15 20 50 100 150 256 500];
%T = [0.25 0.5 0.75 1.00 1.50 2.00 4.00 5.00];
%DeltaT = [0.25 0.5 0.75 1 2 5 8 10 15 20 25 50];

N = 10;
T= [0.67 2.08]; 
DeltaT = [0.10 0.25 0.50 1.00];

realDT = DeltaT.*(1/(2.17e-12*1e15))

path = pwd;
disp(path)
for dd=1:length(DeltaT)
    dt = DeltaT(dd);
    dtpass = realDT(dd);
    newpath = strcat(path,'/data/',num2str(dt));
%    disp(newpath);
    mkdir(newpath);
    for nn=1:length(N)
        newpath1 = strcat(newpath,'/', num2str(N(nn)));
        mkdir(newpath1);
        for tt=1:length(T)
            newpath2 = strcat(newpath1,'/', num2str(T(tt)));
            mkdir(newpath2);
            densityVector = densities(N(nn),2^(1/6));
            for ee=1:length(densityVector)
                newpath3 = strcat(newpath2,'/',num2str(densityVector(ee)));
                mkdir(newpath3);
                fid = fopen(strcat(newpath3,'/energy_temp_vs_time.txt'), 'w');
                fid2 = fopen(strcat(newpath3,'/velocity_distro.txt'), 'w');
                TotalSteps = floor(0.5 + (5.0e-11/(dt*1.0e-15)));
                [allx,L] = Main(N(nn),densityVector(ee),T(tt),dtpass,TotalSteps,fid,fid2);
%                return;
                fclose(fid);fclose(fid2);
                videopathname = strcat(newpath3, '/points.avi');
                mov = Visualization(allx,L,2^(1/6),TotalSteps+1,videopathname);
            end
        end
    end
end
