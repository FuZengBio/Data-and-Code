clear;clc;
filePath=['D:\Paper_recalibration']; %%% Create a file path
cd(filePath);

%%% download data file : "Data_session.rar" and unzip to above file path
%%% download neuron list file : "Neuron_list.rar" and unzip to above file path

load MSTd_PreData.mat; load MSTd_PostData.mat;
filePath=['D:\Paper_recalibration\Data_session\MSTd\'];
  
% load PIVC_PreData.mat; load PIVC_PostData.mat; 
% filePath=['D:\Paper_recalibration\Data_session\PIVC\'];

% load VIP_PreData.mat; load VIP_PostData.mat; 
% filePath=['D:\Paper_recalibration\Data_session\VIP\'];

for kk=1:length(Pre_fileName)
    clear MatFile;
    MatFile=[filePath Pre_fileName{kk}(1:end-4) 'ch' num2str(Pre_SpkChan(kk)) '_HeadingPSTH.mat'];
    load(MatFile);
    
    clear MatFile;
    MatFile=[filePath Post_fileName{kk}(1:end-4) 'ch' num2str(Post_SpkChan(kk)) '_PostHeadingPSTH.mat'];
    load(MatFile);
    
    DeltaOffset=PostHeadingPSTH.Moog_Visual_offset;
    FigureIndex=1;figure(FigureIndex);
    set(FigureIndex,'Position', [300,200 800,600], 'Name', 'Pre & Post');
    orient landscape;text(-0.05,1.05,[Pre_fileName{kk}(1:end-4) 'ch' num2str(Pre_SpkChan(kk)) '  offset: ' num2str(DeltaOffset)]); axis off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:size(HeadingPSTH.firingrate,2) %%% 1:vestibular, 2:visual
        %%% Linear regression for firing rate
        clear x; x=HeadingPSTH.heading';
        clear y; y=HeadingPSTH.firingrate{1,k}';
        clear spont; spont=cell2mat(HeadingPSTH.SpikeBase);
        y=y-spont(:,k);
        clear X; X=[ones(length(x),1) x];
        [b,bint,r,rint,stats] = regress(y,X);
        clear p_val_pre;p_val_pre=stats(3);
        reg_pre= X\y;
        slope_pre=reg_pre(2);
        r_val_pre=b(2);
        
        clear x; x=PostHeadingPSTH.heading';
        clear y; y=PostHeadingPSTH.firingrate{1,k}';
        clear spont; spont=cell2mat(PostHeadingPSTH.SpikeBase);
        y=y-spont(:,k);
        clear X; X=[ones(length(x),1) x];
        [b,bint,r,rint,stats] = regress(y,X);
        clear p_val_post;p_val_post=stats(3);
        reg_post= X\y;
        slope_post=reg_post(2);
        r_val_post=b(2);
        
        r_pre(1,k)=r_val_pre;
        r_post(1,k)=r_val_post;
        p_pre(1,k)=p_val_pre;
        p_post(1,k)=p_val_post;
        
        %%%% baseline change
        Baseline_pre(k) = mean(cell2mat(HeadingPSTH.SpikeBase(1,k)));
        Baseline_post(k) = mean(cell2mat(PostHeadingPSTH.SpikeBase(1,k)));
        Baseline_Change(k)=mean(cell2mat(PostHeadingPSTH.SpikeBase(1,k))-cell2mat(HeadingPSTH.SpikeBase(1,k)))/mean(cell2mat(HeadingPSTH.SpikeBase(1,k)));
        
        %%% neurometric function        
        clear z_pre z_post pre_spikes post_spikes;
        for i=1:size(HeadingPSTH.heading,2)% heading
            pre_spikes(i,:)=HeadingPSTH.relativeFRtrial{k,i};
            post_spikes(i,:)=PostHeadingPSTH.relativeFRtrial{k,i};
        end
        z_pre=(pre_spikes-mean(pre_spikes(:)))/std(pre_spikes(:));
        z_post=(post_spikes-mean(pre_spikes(:)))/std(pre_spikes(:));
        
        for i = 1 : size(HeadingPSTH.heading,2)
            if slope_pre > 0
                Neuro_correct(k,i) =  1-rocN(0,z_pre(i,:),100);
                Neuro_correct_post(k,i) =  1-rocN(0,z_post(i,:),100);
            else
                Neuro_correct(k,i) =  rocN(0,z_pre(i,:),100);
                Neuro_correct_post(k,i) =  rocN(0,z_post(i,:),100);
            end
            
        end
        %%% fit the neurometric function
        addpath Z:\Users\ZengFu\For_Adam\psignifit %% you can also download and install the psignifit4 toolbox (from: https://github.com/wichmann-lab/psignifit) and make sure it is included in the Matlab path
        for i=1:size(Neuro_correct,2)
            fit_data_pre(i, 1) = HeadingPSTH.heading(1,i);
            fit_data_pre(i, 2) = Neuro_correct(k,i);
            fit_data_pre(i, 3) = 100;
            
            fit_data_post(i, 1) = PostHeadingPSTH.heading(1,i);
            fit_data_post(i, 2) = Neuro_correct_post(k,i);
            fit_data_post(i, 3) = 100;
        end
        wichman_pre = pfit(fit_data_pre,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
        Thresh_pre = wichman_pre.params.est(2);
        if Thresh_pre<0 | Thresh_pre> 300
            Thresh_pre = 300;
            wichman_pre.params.est(2) = 300;
        end
        Bias_pre = wichman_pre.params.est(1);
        pre_perf(k,:) = [Bias_pre,Thresh_pre];
        wichman_post = pfit(fit_data_post,'plot_opt','no plot','shape','cumulative gaussian','n_intervals',1,'FIX_LAMBDA',0.001,'sens',0,'compute_stats','false','verbose','false');
        Thresh_post = wichman_post.params.est(2);
        if Thresh_post<0 | Thresh_post> 300
            Thresh_post = 300;
            wichman_post.params.est(2) = 300;
        end
        Bias_post = wichman_post.params.est(1);
        post_perf(k,:) = [Bias_post,Thresh_post];
        
        Neu_Bias_shift(1,k) = post_perf(k,1)-pre_perf(k,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT FIGURE
        %%% behavior
        axes('position',[0.1+(k-1)*0.2 0.7, 0.16 0.2]);
        plot(HeadingPSTH.PsyOri{1,k}(:,1),HeadingPSTH.PsyOri{1,k}(:,2),'ok','Linewidth',2.0);
        hold on;plot(HeadingPSTH.PsyFit{1,k}(:,1),HeadingPSTH.PsyFit{1,k}(:,2),'k','Linewidth',2.0);
        if k==1
            hold on;plot(PostHeadingPSTH.PsyOri{1,k}(:,1),PostHeadingPSTH.PsyOri{1,k}(:,2),'sb','Linewidth',2.0);
            hold on;plot(PostHeadingPSTH.PsyFit{1,k}(:,1),PostHeadingPSTH.PsyFit{1,k}(:,2),'b','Linewidth',2.0);
            title('psy: vest');
            ylabel('Rightward Choices');
        elseif k==2
            hold on;plot(PostHeadingPSTH.PsyOri{1,k}(:,1),PostHeadingPSTH.PsyOri{1,k}(:,2),'sr','Linewidth',2.0);
            hold on;plot(PostHeadingPSTH.PsyFit{1,k}(:,1),PostHeadingPSTH.PsyFit{1,k}(:,2),'r','Linewidth',2.0);
            title('psy: vis');
        end
        xlim([min(HeadingPSTH.PsyOri{1,k}(:,1)),max(HeadingPSTH.PsyOri{1,k}(:,1))]);
        set(gca,'XTick',min(HeadingPSTH.heading):0.5*max(HeadingPSTH.heading):max(HeadingPSTH.heading));set(gca,'XTickLabel',{num2str(min(HeadingPSTH.heading)),num2str(0.5*min(HeadingPSTH.heading)),'0',num2str(0.5*max(HeadingPSTH.heading)),num2str(max(HeadingPSTH.heading))});
        ylim([0,1]);
        
        %%% FR tuning
        axes('position',[0.1+(k-1)*0.2 0.4, 0.16 0.2]);
        plot(HeadingPSTH.heading,HeadingPSTH.relativefiringrate{1,k},'ko','Linewidth',2.0);
        hold on;
        errorbar(HeadingPSTH.heading,HeadingPSTH.relativefiringrate{1,k},HeadingPSTH.relativefiringrate_err{1,k},'k','Linewidth',2.0);
        if k==1
            hold on;plot(PostHeadingPSTH.heading,PostHeadingPSTH.relativefiringrate{1,k},'sb','Linewidth',2.0);
            hold on;errorbar(PostHeadingPSTH.heading,PostHeadingPSTH.relativefiringrate{1,k},PostHeadingPSTH.relativefiringrate_err{1,k},'b','Linewidth',2.0);
            title('tuning curve: vest');
            ylabel('Firing rate');
        elseif k==2
            hold on;plot(PostHeadingPSTH.heading,PostHeadingPSTH.relativefiringrate{1,k},'sr','Linewidth',2.0);
            hold on;errorbar(PostHeadingPSTH.heading,PostHeadingPSTH.relativefiringrate{1,k},PostHeadingPSTH.relativefiringrate_err{1,k},'r','Linewidth',2.0);
            title('tuning curve: vis');
        end
        xlim([min(HeadingPSTH.heading),max(HeadingPSTH.heading)]);
        set(gca,'XTick',min(HeadingPSTH.heading):0.5*max(HeadingPSTH.heading):max(HeadingPSTH.heading));set(gca,'XTickLabel',{num2str(min(HeadingPSTH.heading)),num2str(0.5*min(HeadingPSTH.heading)),'0',num2str(0.5*max(HeadingPSTH.heading)),num2str(max(HeadingPSTH.heading))});
        
        %%% neurometric
        axes('position',[0.1+(k-1)*0.2 0.1, 0.16 0.2]);
        xi = min(HeadingPSTH.PsyFit{1,1}(:,1)) : 0.1 : max(HeadingPSTH.PsyFit{1,1}(:,1));
        plot(HeadingPSTH.PsyOri{1,1}(:,1),Neuro_correct(k,:),'ok','Linewidth',2.0);
        hold on; plot(xi, cum_gaussfit(pre_perf(k,:), xi),'--k','Linewidth',2.0)
        Neu_pre_fit(k,:)=[cum_gaussfit(pre_perf(k,:), xi)];
        if k==1
            hold on;plot(PostHeadingPSTH.PsyOri{1,1}(:,1),Neuro_correct_post(k,:),'bs','Linewidth',2.0);
            hold on;plot(xi, cum_gaussfit(post_perf(k,:), xi),'--b','Linewidth',2.0)
            
        elseif k==2
            hold on;plot(PostHeadingPSTH.PsyOri{1,1}(:,1),Neuro_correct_post(k,:),'rs','Linewidth',2.0);
            hold on;plot(xi, cum_gaussfit(post_perf(k,:), xi),'--r','Linewidth',2.0)
        end
        Neu_post_fit(k,:)=[cum_gaussfit(post_perf(k,:), xi)];
        if k==1
            title('neurometric: vest');
            ylabel('Rightward Choices');
        elseif k==2
            title('neurometric: vis');
        end
        xlabel('Heading Angle');
        xlim([min(HeadingPSTH.PsyOri{1,k}(:,1)),max(HeadingPSTH.PsyOri{1,k}(:,1))]);
        set(gca,'XTick',min(HeadingPSTH.heading):0.5*max(HeadingPSTH.heading):max(HeadingPSTH.heading));set(gca,'XTickLabel',{num2str(min(HeadingPSTH.heading)),num2str(0.5*min(HeadingPSTH.heading)),'0',num2str(0.5*max(HeadingPSTH.heading)),num2str(max(HeadingPSTH.heading))});
        ylim([0,1]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % output some text of basic parameters in the figure
    out=[];
    out = ' Vest | Psy: u | thresh | Neu: u  | thresh ';
    out = strvcat(out, sprintf('--------------------------------------------'));
    OutputValue=[HeadingPSTH.Bias_psy(1) HeadingPSTH.Thresh_psy(1) pre_perf(1,:)];
    out=strvcat(out, sprintf(' %s  | %6.2f | %6.2f | %7.2f | %6.2f  ', 'Pre',OutputValue));
    OutputValue=[PostHeadingPSTH.Bias_psy(1) PostHeadingPSTH.Thresh_psy(1) post_perf(1,:)];
    out=strvcat(out, sprintf(' %s | %6.2f | %6.2f | %7.2f | %6.2f  ', 'Post',OutputValue));
    
    OutputValue=[PostHeadingPSTH.Bias_psy(1)-HeadingPSTH.Bias_psy(1) Neu_Bias_shift(1,1)];
    out=strvcat(out, sprintf(' %s| %6.2f          | %7.2f | %6.2f  ', 'Shift',OutputValue));
    figure(FigureIndex);axes('position',[0.60 0.82 0.2 0.1]);
    text(0,0.9,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');
    set(gca,'box','off','visible','off');
    
    out=[];
    out = ' Vis  | Psy: u | thresh | Neu: u  | thresh ';
    out = strvcat(out, sprintf('--------------------------------------------'));
    OutputValue=[HeadingPSTH.Bias_psy(2) HeadingPSTH.Thresh_psy(2) pre_perf(2,:)];
    out=strvcat(out, sprintf(' %s  | %6.2f | %6.2f | %7.2f | %6.2f  ', 'Pre',OutputValue));
    OutputValue=[PostHeadingPSTH.Bias_psy(2) PostHeadingPSTH.Thresh_psy(2) post_perf(2,:)];
    out=strvcat(out, sprintf(' %s | %6.2f | %6.2f | %7.2f | %6.2f  ', 'Post',OutputValue));
    
    OutputValue=[PostHeadingPSTH.Bias_psy(2)-HeadingPSTH.Bias_psy(2) Neu_Bias_shift(1,2)];
    out=strvcat(out, sprintf(' %s| %6.2f          | %7.2f | %6.2f  ', 'Shift',OutputValue));
    
    figure(FigureIndex);axes('position',[0.60 0.65 0.2 0.1]);
    text(0,0.9,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');
    set(gca,'box','off','visible','off');
    
    out=[];
    out = '                      Vest      vis    ';
    out = strvcat(out, sprintf('--------------------------------------------'));
    OutputValue=[  Baseline_Change  ];
    out=strvcat(out, sprintf(' %s  | %6.2f | %6.2f   ', 'Baseline Change',OutputValue));
    out = strvcat(out, sprintf('--------------------------------------------'));
    figure(FigureIndex);axes('position',[0.60 0.47 0.2 0.1]);
    text(0,0.9,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');
    set(gca,'box','off','visible','off');
    
    out=[];
    out = '       r(vest)  p(vest)     r(vis)  p(vis)   ';
    out = strvcat(out, sprintf('--------------------------------------------'));
    OutputValue=[r_pre(1,1) p_pre(1,1)  r_pre(1,2) p_pre(1,2)  ];
    out=strvcat(out, sprintf(' %s  | %6.2f | %6.2f | %6.2f | %6.2f  ', 'Pre',OutputValue));
    OutputValue=[r_post(1,1) p_post(1,1)  r_post(1,2) p_post(1,2)  ];
    out=strvcat(out, sprintf(' %s  | %6.2f | %6.2f | %6.2f | %6.2f  ', 'Post',OutputValue));
    
    figure(FigureIndex);axes('position',[0.60 0.34 0.2 0.1]);
    text(0,0.9,out,'fontsize',8,'fontname','courier','horizontalalignment','left','verticalalignment','top');
    set(gca,'box','off','visible','off');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % if you wish to save the plots, uncomment the followig line:
    %%% Create a file path to save the figures
    %%% the OutputPath had better be consistent with the filePath for one brain area
    %%% if you wish to save the plots, uncomment the followig line:
    
%         OutputPath=['D:\Paper_recalibration\Results\MSTd\Figure\'];
%     %     OutputPath=['D:\Paper_recalibration\Results\PIVC\Figure\'];
%     %     OutputPath=['D:\Paper_recalibration\Results\VIP\Figure\'];
%     
%     set(gcf, 'PaperOrientation', 'portrait');
%     saveas(gcf,[OutputPath Pre_fileName{kk}(1:end-6) 'ch' num2str(Pre_SpkChan(kk)) '.png'],'png');
%     close(FigureIndex);
    
end








