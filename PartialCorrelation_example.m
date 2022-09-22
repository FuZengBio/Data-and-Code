%%% download data file : "Fig7_PartialCorrelation_example.mat" 
%%% download code: "MR_Heading.m" and make sure it is included in the Matlab path
for k=1:size(Spike_rate_relative,1) %%% 1:vestibular, 2:visual   
    tmpp = cell2mat(heading_regr(k,:))';
    XX{k}(:,1) = tmpp./sqrt(sum(tmpp.^2)/length(tmpp));
    XX{k}(:,2) = cell2mat(choice_regre(k,:))';
    YY{k}=cell2mat(Spike_rate_relative(k,:))';
    
    [parresult(k,1),resp_hat_tmp,p_h] = MR_Heading(YY{k},XX{k},1);
    r_relative_heading(k)= parresult(k,1);
    p_relative_parcorr(k,1)=p_h;
    [parresult(k,2),resp_hat_tmp,p_c] = MR_Heading(YY{k},XX{k},2);
    r_relative_choice(k)= parresult(k,2);
    p_relative_parcorr(k,2)=p_c;
end