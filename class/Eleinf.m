classdef Eleinf % functions to draw plots based on eleinf

    methods(Static)

        function draw_details_farfrag_2d_distribution(dir_details,all_eleinf_include_bos_tut_spe)%dir_farfrags,all_eleinf_include_bos_tut_spe)

            dbstop if error
           all_eleinf = all_eleinf_include_bos_tut_spe;

           for n = 1: length(all_eleinf)
               all_eleinf(n).fullid = sprintf('%s-%02u',all_eleinf(n).songname,all_eleinf(n).fragid);
           end

            dir_details = "E:\Shared_Stimuli_Pool\Detail-CON-B346-03";
            wavnames = Extract.filename(dir_details,'*.wav');
            fragnames = [wavnames{[find(~cellfun(@isempty,regexp([wavnames{:}].','Frag|frag')))]}].';
            shortnames = {};
            for n = 1:length(fragnames)
                [~,shortnames{n},~] = fileparts(fragnames(n)); % fraglist of the stimuli folder
            end

            detail_id = [];

            for n = 1: length(all_eleinf)
                if ~isempty(find(~cellfun(@isempty,regexp([shortnames{:}],all_eleinf(n).fullid))))
                    detail_id = [detail_id,n];
                end
            end

            detail_eleinf = all_eleinf(detail_id);


            figure
            hold on
            scatter([all_eleinf.coor_1].',[all_eleinf.coor_2].',[],'k','filled');

            scatter([detail_eleinf.coor_1].',[detail_eleinf.coor_2].',[],'b','filled');
            hold off

        end
        
        function draw_conspe_eleinf_vs_alleleinf_2d_distribution
            % draw all_eleinf 分布
            conspe_ids = find(~cellfun(@isempty, regexp([all_eleinf.songname].','CON|SPE')));
            ote_ids = find(~cellfun(@isempty, regexp([all_eleinf.songname].','OTE')));

            conspe_eleinf = all_eleinf(conspe_ids);
            ote_eleinf = all_eleinf(ote_ids);


            figure;
            hold on
            scatter([ote_eleinf.coor_1].',[ote_eleinf.coor_2].',[],'k','filled');
            scatter([conspe_eleinf.coor_1].',[conspe_eleinf.coor_2].',[],'g','filled');
            %hold off

            num_all_songs = length(unique([all_eleinf.songname].'));
            num_conspe = length(unique([conspe_eleinf.songname].'));

            set(gca,'xticklabel',{[]})
            xlabel('Dim-1')
            ylabel('Dim-2')

            title(sprintf('%u representatives from %u songs',num_conspe,num_all_songs));

        end
    
    
    end




end


