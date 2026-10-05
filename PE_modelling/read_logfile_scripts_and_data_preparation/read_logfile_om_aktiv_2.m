function [all_data,learn_data,resp_ti_data,stat_outcome]=read_logfile_om_aktiv_2(vp,filename)
%vp = vp-Kürzel als string
%filename = kompletter filename OHNE extension als string



path='logfiles\';
ext='.log';

fullfilename=strcat(path,filename,ext);
                
fid = fopen(fullfilename,'r');
zeile = fgets(fid);

trialnum_ges=0;
trialnum_learn=0;
trialnum_learn_session=0;
trialnum_test=0;
num_session=0;

pattern_presented=0;
test_pattern_presented=0;

ende_uebung=0; %wenn die Ende_Uebung Folie präsentiert wurde geht das Experiment erst richtig los.

learn_session=0;
test_session=0;

start_session=0;


while zeile >= 0
   %zeile = fgets(fid);
   if zeile < 0;break;end;
   [token,zeile] = strtok(zeile); %suche bis zum ersten Space
   								  %in token steht das erste Wort, in 
        token;                          %zeile steht der Rest, inkl. d. ersten Space   
    if (strcmp(token,vp)) %level1 %wenn das erste Wort gleich der VP ist
        token;                
        [token,zeile] = strtok(zeile); %2.Level: hier steht auf jeden Fall der Trial
            
            %trialnum=num2str(token(1,1));
                
            [token,zeile] = strtok(zeile); %3.Level: hier steht die Art des Ereignisses
            token;
            if (strcmp(token,'Picture')) & start_session==0 % hier sind jetzt 2 Loops nötig - die erste, damit das erste start_learn_session entdeckt wird
                [token,zeile] = strtok(zeile); %
                token; 
                if strcmp(token,'fixlearn_start_learn_session')
                    num_session=num_session+1;
                    start_session=1;
                    
                    if num_session==1
                        trialnum_learn=0;
                        trialnum_ges=0;
                    end
                    
                    trialnum_learn_session=0;
                    trialnum_test=0;
                    
                end
            end
               
            
            if (strcmp(token,'Picture')) & num_session>0 % diese zweite Loop soll nur starten, wenn Ende_Uebung schon präsentiert wurde
                [token,zeile] = strtok(zeile);
                 
                
                if length(token)>7 & strcmp(token(1,1:8),'fixlearn') & start_session==1
                   pattern_presented=0; 
                    
                   test_pattern_presented=0;
                   
                end
                
                if strcmp(token,'pause')
                      start_session=0;
                end    
                                   
                if length(token)>=4
                    if strcmp(token(1,1),'p') & strcmp(token(1,3),'p')  
                        if length(token)==4
                            test_pattern_presented=1;
                                                       
                            trialnum_ges=trialnum_ges+1;
                                      
                            trialnum_test=trialnum_test+1;
                            all_data{trialnum_ges,2}='test';
                            all_data{trialnum_ges,4}=trialnum_test;                   
                   
                            all_data{trialnum_ges,1}=vp;
                            all_data{trialnum_ges,3}=trialnum_ges;
                            all_data{trialnum_ges,6}=9; %default-Einstellung, wenn keine Response erfolgt
                            all_data{trialnum_ges,7}=9; 
                            all_data{trialnum_ges,8}=99999;
                            %folgender Wert bleibt immer gleich
                            all_data{trialnum_ges,9}=9; 
                   
                        else
                            pattern_presented=1;
                            
                            %if trialnum_learn_session<60
                                trialnum_ges=trialnum_ges+1;                                           
                                trialnum_learn=trialnum_learn+1;   
                                trialnum_learn_session=trialnum_learn_session+1;
                   
                                all_data{trialnum_ges,1}=vp;
                                all_data{trialnum_ges,2}=num_session;
                                all_data{trialnum_ges,3}=trialnum_ges;
                                all_data{trialnum_ges,4}=trialnum_learn_session;
                                all_data{trialnum_ges,6}=9; %default-Einstellung, wenn keine Response erfolgt
                                all_data{trialnum_ges,7}=9; %default-Einstellung, wenn keine Response erfolgt
                                all_data{trialnum_ges,8}=99999;
                                all_data{trialnum_ges,9}=9;
                                all_data{trialnum_ges,10}=9; %points
                            %end
                        end
                            pattern=token(1,2);                                                
                            trialnum_ges;
                            token;
                            all_data{trialnum_ges,5}=str2num(pattern);
                            
                        [token,zeile] = strtok(zeile);
                        pattern_time(trialnum_ges,1)=str2num(token)/10;
                                                                        
                    end
                end                     
                                                                                                        
                
                if strcmp(token,'right')
                        all_data{trialnum_ges,9}=11;
                        
                        if trialnum_learn_session>1;
                            all_data{trialnum_ges,10}=all_data{trialnum_ges-1,10}+20;
                       elseif trialnum_learn_session==1
                           all_data{trialnum_ges,10}=20;
                       end;   
                       
                elseif strcmp(token,'om_right')
                        all_data{trialnum_ges,9}=1;
                        
                        if trialnum_learn_session>1;
                            all_data{trialnum_ges,10}=all_data{trialnum_ges-1,10};
                       elseif trialnum_learn_session==1
                           all_data{trialnum_ges,10}=0;
                       end;          
                       
                elseif strcmp(token,'wrong')
                       all_data{trialnum_ges,9}=-11;
                       
                       if trialnum_learn_session>1;
                            all_data{trialnum_ges,10}=all_data{trialnum_ges-1,10}-20;
                       elseif trialnum_learn_session==1
                           all_data{trialnum_ges,10}=-20;
                       end;
                       
                elseif strcmp(token,'om_wrong')
                       all_data{trialnum_ges,9}=-1;
                       
                       if trialnum_learn_session>1                            
                           all_data{trialnum_ges,10}=all_data{trialnum_ges-1,10};
                       elseif trialnum_learn_session==1
                           all_data{trialnum_ges,10}=0;
                       end; 
                            
                elseif strcmp(token,'no_resp_no_rew')
                       all_data{trialnum_ges,6}=9;       
                       all_data{trialnum_ges,7}=9;
                       all_data{trialnum_ges,8}=9999;
                       all_data{trialnum_ges,9}=9;
                       
                       if trialnum_learn_session>1
                           all_data{trialnum_ges,10}=all_data{trialnum_ges-1,10};
                       elseif trialnum_learn_session==1
                           all_data{trialnum_ges,10}=0;
                       end;
                end
                
                if strcmp(token,'fixtest') |   strcmp(token,'pause_vor_test')
                    test_pattern_presented=0; 
                            
                end
                                            
                
            elseif (strcmp(token,'Response')) & start_session==1                
                    
                if pattern_presented==1                   
                    [token,zeile] = strtok(zeile);                    
                    resp=str2num(token); 
                    
                    trialnum_ges;
                    resp;
                    
                    if resp>110
                        all_data{trialnum_ges,7}=1;
                    else
                        all_data{trialnum_ges,7}=0;
                    end
                    
                    if resp==101 | resp==111
                        all_data{trialnum_ges,6}=1;
                    elseif resp==102 | resp==112
                        all_data{trialnum_ges,6}=2;
                    end
                    
                    [token,zeile] = strtok(zeile); 
                    resp_time(trialnum_ges,1)=str2num(token)/10;
                    all_data{trialnum_ges,8}=(resp_time(trialnum_ges,1)-pattern_time(trialnum_ges,1));
                     pattern_presented=0; %nächste Reaktion soll nicht mehr gezählt werden
                end
                
                if test_pattern_presented==1                  
                    [token,zeile] = strtok(zeile);                    
                    resp=str2num(token);
                    trialnum_ges;                   
                                       
                    if resp>110
                        all_data{trialnum_ges,7}=1;
                    else
                        all_data{trialnum_ges,7}=0;
                    end                                                                               
                    
                    if resp==101 | resp==111
                        all_data{trialnum_ges,6}=1;
                    elseif resp==102 | resp==112
                        all_data{trialnum_ges,6}=2;
                    end
                    
                        
                        
                    [token,zeile] = strtok(zeile); 
                    resp_time(trialnum_ges,1)=str2num(token)/10;
                    resp_time(trialnum_ges,1);
                    all_data{trialnum_ges,8}=(resp_time(trialnum_ges,1)-pattern_time(trialnum_ges,1));  
                    test_pattern_presented=0; %nächste Reaktion soll nicht mehr gezählt werden
                end                               
              
            end
                        
    end
    zeile=fgets(fid);   
end

fclose(fid);

size(all_data);
all_data;

numtrials_learn=trialnum_learn;
numtrials_test=trialnum_test;

stimuli=all_data(:,5);
response_buttons=all_data(:,6);
responses=all_data(:,7);
resp_times=all_data(:,8);
feedback=all_data(:,9);

num_misses=0;
for g=1:length(responses)
    if responses{g,1}==9
    num_misses=num_misses+1;
    end
end



for i=1:trialnum_ges
    
          
    stimuli{i,1};
    st_le(i,1)=stimuli{i,1};
    resp_bu_le(i,1)=response_buttons{i,1};
    resp_le(i,1)=responses{i,1};    
    resp_ti_le(i,1)=resp_times{i,1};
    fb_le(i,1)=feedback{i,1};
    
end

   
for k=1:numtrials_learn/60
    st_le_block=st_le((k-1)*60+1:k*60,1);
    resp_bu_le_block=resp_bu_le((k-1)*60+1:k*60,1);
    resp_ti_le_block=resp_ti_le((k-1)*60+1:k*60,1);
    
    resp_bu_le_block_a=resp_bu_le_block(st_le_block(:,1)==1); 
    acc_resp_le_block_a=resp_bu_le_block_a(resp_bu_le_block_a(:,1)==2);
    resp_ti_le_block_a=resp_ti_le_block(st_le_block(:,1)==1 & resp_bu_le_block(:,1)~=9); 
    
    resp_bu_le_block_b=resp_bu_le_block(st_le_block(:,1)==2); 
    acc_resp_le_block_b=resp_bu_le_block_b(resp_bu_le_block_b(:,1)==1);
    resp_ti_le_block_b=resp_ti_le_block(st_le_block(:,1)==2 & resp_bu_le_block(:,1)~=9);
    
    resp_bu_le_block_c=resp_bu_le_block(st_le_block(:,1)==3); 
    acc_resp_le_block_c=resp_bu_le_block_c(resp_bu_le_block_c(:,1)==1);
    resp_ti_le_block_c=resp_ti_le_block(st_le_block(:,1)==3 & resp_bu_le_block(:,1)~=9);
    
    resp_bu_le_block_d=resp_bu_le_block(st_le_block(:,1)==4); 
    acc_resp_le_block_d=resp_bu_le_block_d(resp_bu_le_block_d(:,1)==2);
    resp_ti_le_block_d=resp_ti_le_block(st_le_block(:,1)==4 & resp_bu_le_block(:,1)~=9);
    
    resp_bu_le_block_e=resp_bu_le_block(st_le_block(:,1)==5); 
    acc_resp_le_block_e=resp_bu_le_block_e(resp_bu_le_block_e(:,1)==1);
    resp_ti_le_block_e=resp_ti_le_block(st_le_block(:,1)==5 & resp_bu_le_block(:,1)~=9);
    
    resp_bu_le_block_f=resp_bu_le_block(st_le_block(:,1)==6); 
    acc_resp_le_block_f=resp_bu_le_block_f(resp_bu_le_block_f(:,1)==1);
    resp_ti_le_block_f=resp_ti_le_block(st_le_block(:,1)==6 & resp_bu_le_block(:,1)~=9);
    
    misses_block=resp_bu_le_block(resp_bu_le_block(:,1)==9);
    
    learn_data(k,1)=length(acc_resp_le_block_a);
    learn_data(k,2)=length(acc_resp_le_block_b);
    learn_data(k,3)=length(acc_resp_le_block_c);
    learn_data(k,4)=length(acc_resp_le_block_d);
    learn_data(k,5)=length(acc_resp_le_block_e);
    learn_data(k,6)=length(acc_resp_le_block_f);
    % combined conditions
    learn_data(k,7)=length(misses_block);
    
    resp_ti_data(k,1)=mean(resp_ti_le_block_a);
    resp_ti_data(k,2)=mean(resp_ti_le_block_b);
    resp_ti_data(k,3)=mean(resp_ti_le_block_c);
    resp_ti_data(k,4)=mean(resp_ti_le_block_d);
    resp_ti_data(k,5)=mean(resp_ti_le_block_e);
    resp_ti_data(k,6)=mean(resp_ti_le_block_f);
    
end


%check Contingency
resp_bu_le_a=resp_bu_le(st_le(:,1)==1);
resp_bu_le_b=resp_bu_le(st_le(:,1)==2);
resp_bu_le_c=resp_bu_le(st_le(:,1)==3);
resp_bu_le_d=resp_bu_le(st_le(:,1)==4);
resp_bu_le_e=resp_bu_le(st_le(:,1)==5);
resp_bu_le_f=resp_bu_le(st_le(:,1)==6);

fb_le_a=fb_le(st_le(:,1)==1);
fb_le_b=fb_le(st_le(:,1)==2);
fb_le_c=fb_le(st_le(:,1)==3);
fb_le_d=fb_le(st_le(:,1)==4);
fb_le_e=fb_le(st_le(:,1)==5);
fb_le_f=fb_le(st_le(:,1)==6);

%korrekte Reaktionen
acc_resp_bu_le_a=resp_bu_le_a(resp_bu_le_a(:,1)==2);
acc_resp_bu_le_b=resp_bu_le_b(resp_bu_le_b(:,1)==1);
acc_resp_bu_le_c=resp_bu_le_c(resp_bu_le_c(:,1)==1);
acc_resp_bu_le_d=resp_bu_le_d(resp_bu_le_d(:,1)==2);
acc_resp_bu_le_e=resp_bu_le_e(resp_bu_le_e(:,1)==2);
acc_resp_bu_le_f=resp_bu_le_f(resp_bu_le_f(:,1)==1);

%falsche Reaktionen
accf_resp_bu_le_a=resp_bu_le_a(resp_bu_le_a(:,1)==1);
accf_resp_bu_le_b=resp_bu_le_b(resp_bu_le_b(:,1)==2);
accf_resp_bu_le_c=resp_bu_le_c(resp_bu_le_c(:,1)==2);
accf_resp_bu_le_d=resp_bu_le_d(resp_bu_le_d(:,1)==1);
accf_resp_bu_le_e=resp_bu_le_e(resp_bu_le_e(:,1)==1);
accf_resp_bu_le_f=resp_bu_le_f(resp_bu_le_f(:,1)==2);

%pos FB nach korrekten Reaktionen
acc_pfb_resp_bu_le_a=resp_bu_le_a(resp_bu_le_a(:,1)==2 & fb_le_a(:,1)==11);
acc_pfb_resp_bu_le_b=resp_bu_le_b(resp_bu_le_b(:,1)==1 & fb_le_b(:,1)==1);
acc_pfb_resp_bu_le_c=resp_bu_le_c(resp_bu_le_c(:,1)==1 & fb_le_c(:,1)==11);
acc_pfb_resp_bu_le_d=resp_bu_le_d(resp_bu_le_d(:,1)==2 & fb_le_d(:,1)==1);
acc_pfb_resp_bu_le_e=resp_bu_le_e(resp_bu_le_e(:,1)==2 & fb_le_e(:,1)==11);
acc_pfb_resp_bu_le_f=resp_bu_le_f(resp_bu_le_f(:,1)==1 & fb_le_f(:,1)==1);

%pos FB nach inkorrekten Reaktionen
acc_pfbf_resp_bu_le_a=resp_bu_le_a(resp_bu_le_a(:,1)==1 & fb_le_a(:,1)==11);
acc_pfbf_resp_bu_le_b=resp_bu_le_b(resp_bu_le_b(:,1)==2 & fb_le_b(:,1)==1);
acc_pfbf_resp_bu_le_c=resp_bu_le_c(resp_bu_le_c(:,1)==2 & fb_le_c(:,1)==11);
acc_pfbf_resp_bu_le_d=resp_bu_le_d(resp_bu_le_d(:,1)==1 & fb_le_d(:,1)==1);
acc_pfbf_resp_bu_le_e=resp_bu_le_e(resp_bu_le_e(:,1)==1 & fb_le_e(:,1)==11);
acc_pfbf_resp_bu_le_f=resp_bu_le_f(resp_bu_le_f(:,1)==2 & fb_le_f(:,1)==1);


%neg FB nach korrekten Reaktionen
acc_nfb_resp_bu_le_a=resp_bu_le_a(resp_bu_le_a(:,1)==2 & fb_le_a(:,1)==-1);
acc_nfb_resp_bu_le_b=resp_bu_le_b(resp_bu_le_b(:,1)==1 & fb_le_b(:,1)==-11);
acc_nfb_resp_bu_le_c=resp_bu_le_c(resp_bu_le_c(:,1)==1 & fb_le_c(:,1)==-1);
acc_nfb_resp_bu_le_d=resp_bu_le_d(resp_bu_le_d(:,1)==2 & fb_le_d(:,1)==-11);
acc_nfb_resp_bu_le_e=resp_bu_le_e(resp_bu_le_e(:,1)==2 & fb_le_e(:,1)==-1);
acc_nfb_resp_bu_le_f=resp_bu_le_f(resp_bu_le_f(:,1)==1 & fb_le_f(:,1)==-11);

%neg FB nach inkorrekten Reaktionen
acc_nfbf_resp_bu_le_a=resp_bu_le_a(resp_bu_le_a(:,1)==1 & fb_le_a(:,1)==-1);
acc_nfbf_resp_bu_le_b=resp_bu_le_b(resp_bu_le_b(:,1)==2 & fb_le_b(:,1)==-11);
acc_nfbf_resp_bu_le_c=resp_bu_le_c(resp_bu_le_c(:,1)==2 & fb_le_c(:,1)==-1);
acc_nfbf_resp_bu_le_d=resp_bu_le_d(resp_bu_le_d(:,1)==1 & fb_le_d(:,1)==-11);
acc_nfbf_resp_bu_le_e=resp_bu_le_e(resp_bu_le_e(:,1)==1 & fb_le_e(:,1)==-1);
acc_nfbf_resp_bu_le_f=resp_bu_le_f(resp_bu_le_f(:,1)==2 & fb_le_f(:,1)==-11);

check=[length(fb_le_a) length(fb_le_b) length(fb_le_c) length(fb_le_d) length(fb_le_e) length(fb_le_f);
       length(acc_resp_bu_le_a) length(acc_resp_bu_le_b) length(acc_resp_bu_le_c) length(acc_resp_bu_le_d) length(acc_resp_bu_le_e) length(acc_resp_bu_le_f);
       length(acc_pfb_resp_bu_le_a) length(acc_pfb_resp_bu_le_b) length(acc_pfb_resp_bu_le_c) length(acc_pfb_resp_bu_le_d) length(acc_pfb_resp_bu_le_e) length(acc_pfb_resp_bu_le_f);
       length(acc_nfb_resp_bu_le_a) length(acc_nfb_resp_bu_le_b) length(acc_nfb_resp_bu_le_c) length(acc_nfb_resp_bu_le_d) length(acc_nfb_resp_bu_le_e) length(acc_nfb_resp_bu_le_f);
       length(accf_resp_bu_le_a) length(accf_resp_bu_le_b) length(accf_resp_bu_le_c) length(accf_resp_bu_le_d) length(accf_resp_bu_le_e) length(accf_resp_bu_le_f);
       length(acc_pfbf_resp_bu_le_a) length(acc_pfbf_resp_bu_le_b) length(acc_pfbf_resp_bu_le_c) length(acc_pfbf_resp_bu_le_d) length(acc_pfbf_resp_bu_le_e) length(acc_pfbf_resp_bu_le_f);
       length(acc_nfbf_resp_bu_le_a) length(acc_nfbf_resp_bu_le_b) length(acc_nfbf_resp_bu_le_c) length(acc_nfbf_resp_bu_le_d) length(acc_nfbf_resp_bu_le_e) length(acc_nfbf_resp_bu_le_f)];



    
num_misses;

sum_reward=sum(check(3,[1,3,5]))+sum(check(6,[1,3,5]));
sum_punishment=sum(check(4,[2,4,6]))+sum(check(7,[2,4,6]));

outcome=sum_reward*0.2+sum_punishment*(-0.2);
stat_outcome=[sum_reward,sum_punishment,outcome];
