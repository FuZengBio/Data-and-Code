clear;clc;
filePath=['D:\Paper_recalibration']; %%% Create a file path
cd(filePath);

%%% download data file : "Data_session.rar" and unzip to above file path
%%% download neuron timecourse list file : "Neuron_timecourse_list.rar" and unzip to above file path

load MSTd_PreData_timecourse.mat; load MSTd_PostData_timecourse.mat;
filePath=['D:\Paper_recalibration\Data_session\Timecourse\MSTd\'];

% load PIVC_PreData_timecourse.mat; load PIVC_PostData_timecourse.mat;
% filePath=['D:\Paper_recalibration\Data_session\Timecourse\PIVC\'];

% load VIP_PreData_timecourse.mat; load VIP_PostData_timecourse.mat;
% filePath=['D:\Paper_recalibration\Data_session\Timecourse\VIP\'];

for kk=1:length(Pre_fileName)
    clear MatFile;
    MatFile=[filePath Pre_fileName{kk}(1:end-4) 'ch' num2str(Pre_SpkChan(kk)) '_Pre_timecourse.mat'];
    load(MatFile);
    
    clear MatFile;
    MatFile=[filePath Post_fileName{kk}(1:end-4) 'ch' num2str(Post_SpkChan(kk)) '_Post_timecourse.mat'];
    load(MatFile);
    
    DeltaOffset=Post_timecourse.offset(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:length(Pre_timecourse.stim_type) %%% 1:vestibular, 2:visual
        for n=1:length(Pre_timecourse.time)%timecourse
            %%%%%%%%% Pre
            for i=1:length(Pre_timecourse.heading)
                clear x; x=Pre_timecourse.heading';
                FR(i)=Pre_timecourse.relative_timecourse_mean{k,i}(n);
            end
            clear X; X=[ones(length(x),1) x];
            y=FR';
            [b,bint,r,rint,stats] = regress(y,X);
            p_val_pre(k,n)=stats(3);
            reg_pre= X\y;
            intercept_pre(k,n)=reg_pre(1);
            slope_pre(k,n)=reg_pre(2);
            r_val_pre(k,n)=b(2);
            
            %%%%%%%%% Post
            for i=1:length(Pre_timecourse.heading)
                clear x; x=Post_timecourse.heading';
                FR(i)=Post_timecourse.relative_timecourse_mean{k,i}(n);
            end
            clear X; X=[ones(length(x),1) x];
            y=FR';
            [b,bint,r,rint,stats] = regress(y,X);
            p_val_post(k,n)=stats(3);
            reg_post= X\y;
            intercept_post(k,n)=reg_post(1);
            slope_post(k,n)=reg_post(2);
            r_val_post(k,n)=b(2);
            
            if p_val_pre(k,n)<0.05 & p_val_post(k,n)<0.05
                mean_slope(k,n) = mean([atan(slope_pre(k,n)),atan(slope_post(k,n))]);
                shift(k,n)=(intercept_pre(k,n)-intercept_post(k,n))/mean_slope(k,n);
            else
                shift(k,n)=NaN;
            end
            clear pre_spikes;clear post_spikes;
            for i=1:length(Pre_timecourse.heading) %%% heading
                pre_spikes(i,:)=Pre_timecourse.relative_timecourse{k,i,n};
                post_spikes(i,:)=Post_timecourse.relative_timecourse{k,i,n};
            end
            
            z_pre=(pre_spikes-mean(pre_spikes(:)))/std(pre_spikes(:));
            z_post=(post_spikes-mean(pre_spikes(:)))/std(pre_spikes(:));
            
            for i=1:length(Pre_timecourse.heading)
                if slope_pre > 0
                    Neuro_correct(k,i) =  1-rocN(0,z_pre(i,:),100);
                    Neuro_correct_post(k,i) =  1-rocN(0,z_post(i,:),100);
                else
                    Neuro_correct(k,i) =  rocN(0,z_pre(i,:),100);
                    Neuro_correct_post(k,i) =  rocN(0,z_post(i,:),100);
                end
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%fit the neurometric function
            addpath Z:\Users\ZengFu\For_Adam\psignifit
            for i=1:size(Neuro_correct,2)
                fit_data_pre(i, 1) = Pre_timecourse.heading(1,i);
                fit_data_pre(i, 2) = Neuro_correct(k,i);
                fit_data_pre(i, 3) = 100;
                
                fit_data_post(i, 1) = Post_timecourse.heading(1,i);
                fit_data_post(i, 2) = Neuro_correct_post(k,i);
                fit_data_post(i, 3) = 100;
            end
            
            wichman_pre = pfit(fit_data_pre,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
            Thresh_pre(k,n) = wichman_pre.params.est(2);
            if Thresh_pre(k,n)<0 | Thresh_pre(k,n)> 300
                Thresh_pre(k,n) = 300;
                wichman_pre.params.est(2) = 300;
            end
            Bias_pre(k,n) = wichman_pre.params.est(1);
            
            wichman_post = pfit(fit_data_post,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
            Thresh_post(k,n) = wichman_post.params.est(2);
            if Thresh_post(k,n)<0 | Thresh_post(k,n)> 300
                Thresh_post(k,n) = 300;
                wichman_post.params.est(2) = 300;
            end
            Bias_post(k,n) = wichman_post.params.est(1);
            
            Neu_Bias_shift(k,n)= Bias_post(k,n)-Bias_pre(k,n);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % if you wish to save the plots, uncomment the followig line:
    %%% Create a file path to save the figures
    %%% the OutputPath had better be consistent with the filePath for one brain area
    %%% if you wish to save data for further analysis, uncomment the followig line:
    
%     TimecourseData.Delta=DeltaOffset;
%     TimecourseData.r_pre=r_val_pre;
%     TimecourseData.p_pre=p_val_pre;
%     TimecourseData.r_post=r_val_post;
%     TimecourseData.p_post=p_val_post;    
%     TimecourseData.Bias_pre=Bias_pre;
%     TimecourseData.Thresh_pre=Thresh_pre;
%     TimecourseData.Bias_post=Bias_post;
%     TimecourseData.Thresh_post=Thresh_post;
%     TimecourseData.Neu_Bias_shift=Neu_Bias_shift;
%     
%     OutputPath=['D:\Paper_recalibration\Results\MSTd\Timecourse\'];
%     %OutputPath=['D:\Paper_recalibration\Results\PIVC\Timecourse\'];
%     %OutputPath=['D:\Paper_recalibration\Results\VIP\Timecourse\'];    
%     
%     SaveFileName=[OutputPath Pre_fileName{kk}(1:end-6) 'ch' num2str(Pre_SpkChan(kk)) '_TimecourseData'];
%     save(SaveFileName, 'TimecourseData');  
end



