function [subject_averages] = subjectAverageCalculator(allSubjectData)
%subjectAverageCalculator: This function calculates an average for each
%subject

    % Input: a nx2 vector where the col. 1 is subj # and col. 2 is data

    subjIDs = unique(allSubjectData(:,1));
    subject_averages = zeros(size(subjIDs));

    for i = 1:length(subjIDs)
        subjData = allSubjectData(allSubjectData(:,1) == subjIDs(i), 2);
        subject_averages(i) = mean(subjData, 'omitnan');
    end

end