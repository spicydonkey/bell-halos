% Update figure3: quick and dirty script to edit figure properties from postprocessed data
% 20181005
% DKS
%
% * ONLY edit figure annotations!
%

% load figures to pretty up
h_rabi=openfig('C:\Users\HE BEC\Dropbox\PhD\thesis\projects\bell_epr_steering\fig3\rabi_halo_20180823.fig');
h_ramsey=openfig('C:\Users\HE BEC\Dropbox\PhD\thesis\projects\bell_epr_steering\fig3\ramsey_fringe_halo_20180823.fig');

%% figure changes
% 1. label update
ylabel_new='Population fraction';

figure(h_rabi);
ylabel(ylabel_new);

figure(h_ramsey);
ylabel(ylabel_new);

% 2. axis linewidth
ax_linewidth=1.1;

figure(h_rabi);
set(gca,'LineWidth',ax_linewidth);

figure(h_ramsey);
set(gca,'LineWidth',ax_linewidth);

% 3. renderer
renderer='painters';

set(h_rabi,'Renderer',renderer);
set(h_ramsey,'Renderer',renderer);
