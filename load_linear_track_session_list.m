    if rat == 1 %Clover
        directory = '/home/caitlin/Data//Processed_Data/Clover/linear_track';
        dayFiles = {'20210427/hand_clustered','20210428/hand_clustered','20210429/kilosort_clustered','20210430/kilosort_clustered','20210501/hand_clustered','20210507/kilosort_clustered','20210508/hand_clustered','20210509/kilosort_clustered','20210512/kilosort_clustered','20210517/hand_clustered','20210519/kilosort_clustered'};
    elseif rat == 2 %Bo
        directory = '/home/caitlin/Data//Processed_Data/Bo/linear_track';
        %best decoding
        %dayFiles = {'20210505','20210507','20210508','20210509/kilosort_clustered'} % first 4 sessions are for cell properties or LFP only
        dayFiles = {'20210511/kilosort_clustered','20210512/kilosort_clustered','20210513/kilosort_clustered','20210517/kilosort_clustered'};
    elseif rat == 3 % MEC1
        directory = '/home/caitlin/Data//Processed_Data/MEC1';
        %best decoding:
        %         dayFiles = {'20201030','20201103','20201104'}
        dayFiles = {'20200926/kilosort_clustered','20200927/hand_clustered2'};
    elseif rat ==4
        directory = '/home/caitlin/Data//Processed_Data/Bolt/linear_track';
        %dayFiles = {'20220714_161350/processed/msort_clustered','20220715_063413/processed/msort_clustered','20220718_070721/processed/msort_clustered','20220718_083127/processed/msort_clustered','20220722_073352/processed/msort_clustered_full','20220729_093340/processed'};
        dayFiles = {'20220714/msort_clustered','20220715/msort_clustered','20220718_1/msort_clustered','20220718_2/msort_clustered','20220719_1_opto/msort_clustered_full','20220720/msort_clustered','20220722/msort_clustered_full','20220729/msort_clustered'};
    elseif rat == 5
        directory = '/home/caitlin/Data//Processed_Data/Dash';
        dayFiles = {'20221025_113913', '20221025_170200'};
    elseif rat == 6
        directory = '/home/caitlin/Data//Processed_Data/CM1';
        dayFiles = ...
            {'20230217','20230219','20230219_01','20230220','20230221','20230221_2_bad','20230221_3','20230222_1','20230222_2','20230224_2','20230228','20230301_01','20230302_01','20230302_02','20230303_01','20230303_02','20230305_01','20230305_02'};
    elseif rat == 7 % Janni
        directory = '/home/caitlin/Data//Processed_Data/Janni';
        dayFiles = ...
            {'20100408_run1','20100410_run1','20100410_run2','20100410_run3','20100412_run1','20100413_run1','20100413_run2'};
    elseif rat == 8 % Harpy
        directory = '/home/caitlin/Data//Processed_Data/Harpy';
        dayFiles = ...
            {'20100114_run1','20100115_run1','20100115_run2'};
    elseif rat == 9 % Imp
        directory = '/home/caitlin/Data//Processed_Data/Imp';
        dayFiles = ...
            {'20100222_run1','20100222_run2','20100225_run1'};
    elseif rat == 10 % W-18
        directory = '/home/caitlin/Data//Processed_Data/W18';
        dayFiles = ...
            {'20140817_run1','20140817_run2','20140822_run1'};
    elseif rat == 11 % W-19
        directory = '/home/caitlin/Data//Processed_Data/W19';
        dayFiles = ...
            {'20141105_run1','20141106_run1','20141107_run1','20141107_run2'};
    else
                directory = {};
        dayFiles = ...
            {};
    end