if rat == 1 % Clover
    directory = '/home/caitlin/Data//Processed_Data/Clover/open_field/';
    dayFiles = {'20210419_session1','20210419_session2','20210420','20210421/Mclust_combined_sessions_1_2_3_4','20210422','20210423','20210424'};
% elseif rat == 2 % Bo
%     directory = '/home/caitlin/Data//Processed_Data/Bo/open_field/';
%     dayFiles = {'20210516','20210518','20210519','20210520'};
elseif rat == 4 % Bolt
    directory = '/home/caitlin/Data//Processed_Data/Bolt/open_field/';
    dayFiles = {'20220513_sessions_1-2/msort_clustered_full','20220515/msort_clustered_full','20220516/msort_clustered_full',...
        '20220520/msort_clustered_full','20220523/msort_clustered_full','20220524/msort_clustered_full',...
        '20220606_novel_openfield','20220608_novel_open_field_day2/msort_clustered_full'};
    % The following have been cleaned up/each reward consumption periods is
    % defined:
    dayFiles =  {'20220513_sessions_1-2/msort_clustered_full','20220515/msort_clustered_full','20220516/msort_clustered_full','20220520/msort_clustered_full','20220523/msort_clustered_full','20220524/msort_clustered_full'};
 
elseif rat == 6 % CM1
    directory = '/home/caitlin/Data//Processed_Data/CM1/open_field';
    dayFiles = {'20230306_01','20230306_02','20230307_01','20230307_02','20230308_01','20230309_01','20230310_01','20230312_01','20230313_01','20230315_01','20230317'};

elseif rat == 12 % Billy3
    directory = '/home/caitlin/Data//Processed_Data/Billy3/';
    dayFiles = {'20181205','20181206','20181207','20181208','20181212'};
  
%'20181207',
elseif rat == 13 %Curly2
     directory = '/home/caitlin/Data//Processed_Data/Curly2/';
     dayFiles = {'20191007','20191008','20191009','20191010','20191011','20191013','20191104'};

elseif rat == 14 %Goathe2
  directory = '/home/caitlin/Data//Processed_Data/Goethe2/';
  dayFiles = {'20180607','20180608','20180609','20180614'};

else
    directory = '';
    dayFiles = {};
end