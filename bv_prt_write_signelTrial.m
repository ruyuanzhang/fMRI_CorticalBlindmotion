%make prt file for single tr
%clear all;close all;clc;
%dir = '/motionDotsFMRIPrt';
%cd(dir);
%filename = ;
%load(filename)
%cd
Params = p;

Params.General.Run             = Params.run_nr;
Params.General.ExpName         = Params.ExpName;
Params.General.Subject         = Params.subject_nr;
Params.General.N_direction     = Params.n_direction;
Params.General.N_location      = Params.n_location;
Params.General.Trials_condition= Params.trials_condition;

%prt for individual trials
prtSingleTrial.name = [Params.General.ExpName '_Sub' int2str(Params.General.Subject) '_run' int2str(Params.General.Run) '_SingleTrial.prt'];
prtSingleTrial.NrOfConditions = Params.General.N_direction * Params.General.N_location * Params.General.Trials_condition;
prtSingleTrial.Condition{:}.ntpts=Params.General.Trials_condition;
prtSingleTrial.Condition{1}.name='DirLeft_FieldLeft_trial1'; %direction_field_trial1;
prtSingleTrial.Condition{2}.name='DirLeft_FieldLeft_trial2'; %direction_field_trial2;
prtSingleTrial.Condition{3}.name='DirLeft_FieldLeft_trial3'; %direction_field_trial3;
prtSingleTrial.Condition{4}.name='DirLeft_FieldLeft_trial4'; %direction_field_trial4;
prtSingleTrial.Condition{5}.name='DirLeft_FieldRight_trial1'; %direction_field_trial1;
prtSingleTrial.Condition{6}.name='DirLeft_FieldRight_trial2'; %direction_field_trial2;
prtSingleTrial.Condition{7}.name='DirLeft_FieldRight_trial3'; %direction_field_trial3;
prtSingleTrial.Condition{8}.name='DirLeft_FieldRight_trial4'; %direction_field_trial4;
prtSingleTrial.Condition{9}.name='DirRight_FieldLeft_trial1'; %direction_field_trial1;
prtSingleTrial.Condition{10}.name='DirRight_FieldLeft_trial2'; %direction_field_trial2;
prtSingleTrial.Condition{11}.name='DirRight_FieldLeft_trial3'; %direction_field_trial3;
prtSingleTrial.Condition{12}.name='DirRight_FieldLeft_trial4'; %direction_field_trial4;
prtSingleTrial.Condition{13}.name='DirRight_FieldRight_trial1'; %direction_field_trial1;
prtSingleTrial.Condition{14}.name='DirRight_FieldRight_trial2'; %direction_field_trial2;
prtSingleTrial.Condition{15}.name='DirRight_FieldRight_trial3'; %direction_field_trial3;
prtSingleTrial.Condition{16}.name='DirRight_FieldRight_trial4'; %direction_field_trial4;


% generate the colors (jet is a 64x3 color table blue-green-yellow-red):
colors = 255*jet; close(gcf);
prtSingleTrial.colors = round(linspace(1,length(colors),prtSingleTrial.NrOfConditions));
prtSingleTrial.colors = round(colors(prtSingleTrial.colors,:));

%compute prt for all individual trials
for i =1: prtSingleTrial.NrOfConditions
    prtSingleTrial.Condition{i}.ntpts = 1;
    prtSingleTrial.Condition{i}.color=prtSingleTrial.colors(i,:);
    
    
    condition_temp=fix((i-1)/4)+1;
    rem_temp = rem(i-1,4)+1;
    estart_temp = prt.data(prt.data(:,3)==condition_temp,4);
    eend_temp   = prt.data(prt.data(:,3)==condition_temp,5);
    
    prtSingleTrial.Condition{i}.estart =estart_temp(rem_temp) ;
    prtSingleTrial.Condition{i}.eend = eend_temp(rem_temp);
    
end
prtSingleTrial.data=prt.data;
bv_prt_write (prtSingleTrial);