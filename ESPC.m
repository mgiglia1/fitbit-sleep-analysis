function ESPC(pid)
%function a = ESPC(pid)

    rawGPSPath='/Users/mollygiglia/Thesis/';
    %% Read RAW GPS Location Points i.e. Lp
    File = [rawGPSPath pid '.csv'];
    %disp(File);
    fid1 = fopen(File);
    % <ts,lat,long,ACCURACY> 
    A = textscan(fid1,'%f %f %f %*d','delimiter',',','headerlines',1);
    fclose(fid1);
    [ts,lat,long]=A{1,1:3};
    %ts = ts./1000; % sec
    lat_long = [lat,long];
    N = length(ts);
    

    %disp(lat_long);
   %disp(N);
    
    Tmax = 10*60; % sec
    Tmin = 5*60; % Sec 
    Tmin_opportunistic = Tmin; % Sec 
%     Dmax = .1; % km 
    Dmax = 250; % meter 
    
    %% First cluster Lp into sp
    sp_counter = 0;
    sp = zeros(N,4+1);  % stay points <lat,long,st_ts,end_ts>   true_sp(1)/missing(-1)
    sp_data_missing_Tmax_test = []; %zeros(N,2);
    sp_cluster_Dmax_Ttotal_size = zeros(N,1+3+1); % radius (distSPAN timeSPAN sample_count,i.e., sp_cluster_size)   true_sp(1)/missing(-1)
    i = 1;
    while i<N
%         i
        j = i+1;
        while j<N
%            j
           t = (ts(j)-ts(j-1));
           if(t>Tmax) 
                %% This satisfies spatial bound (Dmax = 250) to form clusters, but we still need to check for temporal bound (Tmin_new)
                t = (ts(j-1)-ts(i)); 
                if(t>Tmin_opportunistic)
                    centroid_lat = median(lat_long(i:j-1,1));
                    centroid_long = median(lat_long(i:j-1,2));
    %                 sp = [sp; [centroid_lat centroid_long ts(i) ts(j-1)]];
                    sp_counter = sp_counter + 1;
                    sp(sp_counter,:) = [centroid_lat centroid_long ts(i) ts(j-1) -1];
                    [d_max,d_ind] = max( max_haversine(lat_long(i:j-1,:),'max') ); % max distance SPAN between any 2 sps of a sp cluster
                    r = max( haversine(lat_long(i:j-1,1),lat_long(i:j-1,2),centroid_lat,centroid_long) ); % radius of the centroid 
                    t = (ts(j-1)-ts(i)); 
                    sp_cluster_Dmax_Ttotal_size(sp_counter,:) = [r d_max t (j-i) -1]; 
                end
                i = j; 
                break; 
           end;
           d = haversine(lat_long(i,1),lat_long(i,2),lat_long(j,1), lat_long(j,2)); % km
           if(d>Dmax)
                t = (ts(j-1)-ts(i)); 
                if(t>Tmin)
                    centroid_lat = median(lat_long(i:j-1,1));
                    centroid_long = median(lat_long(i:j-1,2));
%                     sp = [sp; [centroid_lat centroid_long ts(i) ts(j-1)]];
                    sp_counter = sp_counter + 1;
                    sp(sp_counter,:) = [centroid_lat centroid_long ts(i) ts(j-1) 1];
                    [d_max,d_ind] = max( max_haversine(lat_long(i:j-1,:),'max') ); % max distance SPAN between any 2 sps of a sp cluster
                    r = max( haversine(lat_long(i:j-1,1),lat_long(i:j-1,2),centroid_lat,centroid_long) ); % radius of the centroid 
                    t = (ts(j-1)-ts(i)); 
                    sp_cluster_Dmax_Ttotal_size(sp_counter,:) = [r d_max t (j-i) 1]; 
%                     pause
                end
                i = j;
                break;
           end
           j = j+1;
        end
        if(j==N)
            t = (ts(j-1)-ts(i)); 
            if(t>Tmin)
                centroid_lat = median(lat_long(i:j-1,1));
                centroid_long = median(lat_long(i:j-1,2));
    %             sp = [sp; [centroid_lat centroid_long ts(i) ts(j-1)]];
                sp_counter = sp_counter + 1;
                sp(sp_counter,:) = [centroid_lat centroid_long ts(i) ts(j-1) 2];
                [d_max,d_ind] = max( max_haversine(lat_long(i:j-1,:),'max') ); % max distance SPAN between any 2 sps of a sp cluster
                r = max( haversine(lat_long(i:j-1,1),lat_long(i:j-1,2),centroid_lat,centroid_long) ); % radius of the centroid 
                t = (ts(j-1)-ts(i)); 
                sp_cluster_Dmax_Ttotal_size(sp_counter,:) = [r d_max t (j-i) 2]; 
            end
            i = j;
        end
    end
    %% Re-size sp with ONLY sp clusters and discard non-filled entries
    ind = find(sp(:,1)>0);
%     length(ind)
%     size(sp)
%     size(sp_cluster_Dmax_Tmax_size)
    sp = sp(ind,:);
    sp_cluster_Dmax_Ttotal_size = sp_cluster_Dmax_Ttotal_size(ind,:);

    %% Write Stay Point (sp) clusters into files
    File_Name_Initial = 'sp';
    File_Dir = [rawGPSPath 'Full_cluster_info/sp_median_Tmin' num2str(Tmin/60) 'minutes_Dmax' num2str(Dmax) 'meters'];
    if(~exist(File_Dir,'dir')) 
%             mkdir([File_Dir '\']);
        mkdir([File_Dir '/']);
    end
    File = [File_Dir '/' File_Name_Initial '_pid_' num2str(pid) '.csv'];
    if(size(sp,1)>0)
%         sp
        File
%         pause
    
        fileID = fopen(File,'w+');
        fprintf(fileID,'lat,long,st_ts,end_ts,true_sp(1),radius(m),d_max(m),duration(sec),Lp_count\n');
        for i=1:size(sp,1)
            fprintf(fileID,'%f,%f,%0.0f,%0.0f,%d,%f,%f,%f,%d\n',sp(i,:),sp_cluster_Dmax_Ttotal_size(i,1:end-1));
        end
        fclose(fileID);
    end

end