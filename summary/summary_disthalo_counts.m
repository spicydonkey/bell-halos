function nAtomsDet=summary_disthalo_counts(k)
%nAtomsDet=summary_disthalo_counts(k)
%

nAtomsDet=cellfun(@(x) size(x,1),k);        % number of atoms detected in halo

% histogram of population
figure;
for ii=1:size(k,2)
    hold on;
    histogram(nAtomsDet(:,ii));
    xlabel('num atoms detected');
    ylabel('N');
end

end