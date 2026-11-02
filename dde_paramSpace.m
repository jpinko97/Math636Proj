%parameter space graph M636

B = logspace(-2, 1, 400); %log axes


% overdamped curve - real roots
lnR = 1+exp(-(B+1))./B;
r = exp(lnR);
aLow = B.*r;

maskLow = aLow < 1e3 & isfinite(aLow); %no divergence
bLow = B(maskLow);
aLow = aLow(maskLow);


% limit cycle curve - imaginary eigenvalues
Nw   = 8000;
ep = 1e-4;
w = linspace(pi/2+ep, pi-ep, Nw);

tw = tan(w);
cw = cos(w);

bUp = -w./tw;              
rUp = exp(1-1./cw);      
AUp = bUp.*rUp;     

% within lims
mask_top = (bUp > 1e-2) & (bUp < 10) & (AUp > 2) & (AUp < 1e3) & ...
    isfinite(bUp) & isfinite(AUp);
bUp = bUp(mask_top);
AUp = AUp(mask_top);

% plot
figure;
hold on;
loglog(bLow,  aLow,  'k-', 'LineWidth', 1.5);
loglog(bUp, AUp, 'k-', 'LineWidth', 1.5);
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([1e-2, 10]);
ylim([2, 1e3]);
xlabel('\delta T_D');
ylabel('P T_D');
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'Box', 'on');
set(gcf, 'Color', 'w');
grid on;
title('Qualitative behaviour of the model given parameters ( \delta T_D , P T_D )', 'Color','k');


% labels
text(0.015, 5, 'Stable overdamped:', 'FontSize', 10, 'Color', 'k');
text(0.015, 3.8,'no cycles', 'FontSize', 10, 'Color', 'k');
text(2.2, 5, 'Stable overdamped:', 'FontSize', 10, 'Color', 'k');
text(2.2, 3.8,'no cycles', 'FontSize', 10, 'Color', 'k');

text(0.14, 20, 'Stable', 'FontSize', 10, 'Color', 'k');
text(0.14, 15, 'underdamped:', 'FontSize', 10, 'Color', 'k');
text(0.14, 12, 'driven quasi-cycles', 'FontSize', 10, 'Color', 'k');

text(0.75, 200, 'Deterministically', 'FontSize', 10, 'Color', 'k');
text(0.75, 130, 'unstable:', 'FontSize', 10, 'Color', 'k');
text(0.75, 85, 'limit cycles', 'FontSize', 10, 'Color', 'k');

hold off;
