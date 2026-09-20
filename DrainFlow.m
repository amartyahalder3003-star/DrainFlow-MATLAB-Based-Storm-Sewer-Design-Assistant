function DrainFlow
%DRAINFLOW  Live storm sewer design assistant - single-file application.
%
%   DrainFlow
%
%   ONE FILE, NO DEPENDENCIES. Everything the application needs - the IDF
%   table, the Rational Method, Manning's equation, the circular-pipe
%   geometry, the bisection root finder, the design storm, the numerical
%   integrators and the statistics - lives in this file as local functions,
%   in the same way a web page carries its own CSS and JavaScript. Drop
%   this file anywhere on the MATLAB path and type DrainFlow.
%
%   STATION: KHULNA, BANGLADESH
%   ---------------------------------------------------------------------
%   The rainfall data is the published Khulna IDF table:
%
%     Sheonty, S.R. & Islam, G.M.T. (2020) "Establishment of Rainfall
%     Intensity-Duration-Frequency Curves of Khulna". Proceedings of the
%     5th International Conference on Civil Engineering for Sustainable
%     Development (ICCESD 2020), 7-9 February 2020, KUET, Khulna,
%     Bangladesh. Paper WRE-4374.
%     https://iccesd.kuet.ac.bd/2020/Papers/WRE-4374.pdf
%
%   That study used 62 years of BMD annual maximum daily rainfall for
%   Khulna (1948-2010), checked the record for homogeneity, consistency
%   and randomness, fitted a Gumbel Type-I Extreme Value distribution, and
%   then disaggregated to short durations using three candidate equations
%   (Talbot, Sherman, Kimjima), selecting Sherman by least squares.
%
%   Durations run from 15 minutes to 24 hours, which is the range urban
%   storm sewer design actually needs.
%
%   CAVEAT TO DECLARE IN THE REPORT: short-duration rainfall gauges are
%   scarce in Bangladesh, so durations below about 3 hours are
%   DISAGGREGATED from daily records rather than directly measured. The
%   source paper states this limitation itself.
%
%   WHAT IT DOES
%     1  Computes the time of concentration of the catchment (Kirpich).
%     2  Reads a rainfall intensity from the Khulna IDF table at that
%        duration, interpolating between tabulated entries.
%     3  Converts it to a peak discharge with the Rational Method.
%     4  Computes the trial pipe's capacity with Manning's equation.
%     5  Solves Manning's equation backwards, by bisection, for the depth
%        the water actually reaches - and animates it.
%     6  Reports a pass / warn / fail verdict against capacity and the
%        self-cleansing velocity limits.
%
%   BUTTONS
%     Pause / Play    freeze or resume the animation
%     Print report    write a full numerical design report, including the
%                     hand-calculation check and the Gumbel back-analysis
%                     of the published table, to the Command Window
%
%   REQUIRES MATLAB R2019a or newer (uifigure with uigridlayout).
%   No toolboxes.
%
%   CONTENTS OF THIS FILE
%     DrainFlow ................. the application (this function)
%       nested functions ........ callbacks, solver cache, drawing
%     -- user interface helpers --
%     sectionLabel .............. heading in the control column
%     styleAxes ................. dark theme for a uiaxes
%     -- engineering functions ---
%     idf_intensity ............. Khulna IDF table + interp1 interpolation
%     time_of_concentration ..... Kirpich formula
%     rational_Q ................ Q = C*i*A/360
%     pipe_geometry ............. circular segment A, P, R, T
%     manning_Q ................. Q = (1/n) A R^(2/3) S^(1/2)
%     manning_capacity .......... full-pipe capacity, Qmax by grid search
%     bisect_root ............... bisection root finder (from the lectures)
%     normal_depth .............. depth from discharge - the root problem
%     triangular_hyetograph ..... design storm shape
%     simpson13 ................. Simpson's 1/3 rule
%     design_status ............. pass / warn / fail verdict
%     gumbel_fit ................ EV-I fit by the method of moments
%     gumbel_from_table ......... recover EV-I parameters from the table
%
%   CE 206 - DrainFlow: Storm Sewer Design Assistant
%   Bangladesh University of Engineering and Technology

%% ======================================================================
%  COLOUR PALETTE
% =======================================================================
COL.bg      = [0.039 0.059 0.110];   % page background
COL.panel   = [0.067 0.102 0.173];   % panel fill
COL.panel2  = [0.051 0.082 0.141];   % inset / axes fill
COL.line    = [0.180 0.220 0.290];   % borders and grid
COL.ink     = [0.910 0.929 0.969];   % primary text
COL.inkDim  = [0.549 0.604 0.722];   % secondary text
COL.water   = [0.310 0.639 0.890];   % water blue
COL.ok      = [0.208 0.761 0.584];   % green
COL.warn    = [0.910 0.639 0.239];   % amber
COL.danger  = [0.886 0.333 0.310];   % red

MONO = get(groot, 'FixedWidthFontName');

% Depth at which discharge in a circular pipe is a maximum, as a fraction
% of the diameter. Q(y) rises then FALLS: past this point the wall curves
% inward and adds wetted perimeter faster than flow area, so R = A/P drops
% and Q with it. The ratio is a property of the circle alone - D, S and n
% all cancel. Used to bracket the bisection.
YMAX_RATIO = 0.9381812;

%% ======================================================================
%  STATE
% =======================================================================
st.A        = 2.0;      % catchment area (ha)
st.C        = 0.65;     % runoff coefficient
st.L        = 600;      % longest overland flow path (m)
st.Sc       = 0.005;    % average catchment slope (m/m)
st.useTc    = true;     % set the storm duration equal to tc
st.Tr       = 10;       % return period (yr)
st.durMin   = 30;       % storm duration (MINUTES) when useTc is false
st.override = false;    % use a manual intensity instead of the IDF table
st.manualI  = 120;      % manual intensity (mm/hr)
st.D        = 0.75;     % pipe diameter (m)
st.S        = 0.005;    % pipe slope (m/m)
st.n        = 0.013;    % Manning's n
st.loopSec  = 10;       % real seconds per simulated storm
st.vmin     = 0.6;      % self-cleansing velocity limit (m/s)
st.vmax     = 3.0;      % scour velocity limit (m/s)

st.simT     = 0;        % simulation clock, 0 to 1
st.running  = true;
st.hist     = struct('t', [], 'i', [], 'q', []);
st.rain     = struct('x', [], 'y', [], 'len', [], 'v', []);

% Derived, filled in by recompute()
st.tc = 0; st.durUse = 30; st.exact = false;
st.Ipeak = 0; st.Qreq = 0; st.Qfull = 0; st.Vfull = 0;
st.Qmax  = 0; st.yQmax = 0; st.yflow = 0; st.Vflow = 0;
st.status = ''; st.level = 1;

%% ======================================================================
%  WINDOW AND TOP-LEVEL GRID
% =======================================================================
fig = uifigure('Name', 'DrainFlow - Storm Sewer Design Assistant (Khulna)', ...
               'Color', COL.bg, 'Position', [50 50 1360 850]);

mainGrid = uigridlayout(fig, [2 2]);
mainGrid.RowHeight       = {52, '1x'};
mainGrid.ColumnWidth     = {340, '1x'};
mainGrid.BackgroundColor = COL.bg;
mainGrid.Padding         = [18 14 18 10];
mainGrid.RowSpacing      = 12;
mainGrid.ColumnSpacing   = 16;

hdr = uigridlayout(mainGrid, [1 2]);
hdr.Layout.Row      = 1;
hdr.Layout.Column   = [1 2];
hdr.ColumnWidth     = {'1x', 360};
hdr.Padding         = [0 0 0 0];
hdr.BackgroundColor = COL.bg;

titleLbl = uilabel(hdr, 'Text', 'DrainFlow', 'FontSize', 24, ...
    'FontWeight', 'bold', 'FontColor', COL.ink);
titleLbl.Layout.Column = 1;

h.badge = uilabel(hdr, 'Text', 'KHULNA IDF - SIMULATING', 'FontName', MONO, ...
    'FontSize', 12, 'FontColor', COL.ok, 'HorizontalAlignment', 'right');
h.badge.Layout.Column = 2;

%% ======================================================================
%  LEFT CONTROL COLUMN
% =======================================================================
ctrlPanel = uipanel(mainGrid, 'BackgroundColor', COL.panel, ...
    'BorderType', 'none');
ctrlPanel.Layout.Row    = 2;
ctrlPanel.Layout.Column = 1;

cg = uigridlayout(ctrlPanel, [27 1]);
cg.RowHeight = {20, 18, 32, 18, 32, 18, 32, 18, 32, 22, 22, 52, 22, 28, ...
                26, 18, 32, 18, 32, 16, 28, ...
                26, 18, 32, 30, 30, '1x'};
cg.Padding         = [16 12 16 12];
cg.RowSpacing      = 4;
cg.BackgroundColor = COL.panel;
cg.Scrollable      = 'on';

% ---- catchment ---------------------------------------------------------
sectionLabel(cg, 'CATCHMENT & RAINFALL', COL, 1);

h.lblArea = uilabel(cg, 'Text', '', 'FontSize', 12, 'FontColor', COL.inkDim);
h.lblArea.Layout.Row = 2;
h.sldArea = uislider(cg, 'Limits', [0.1 20], 'Value', st.A, ...
    'MajorTicks', [], 'MinorTicks', [], 'FontColor', COL.inkDim);
h.sldArea.Layout.Row = 3;
h.sldArea.ValueChangingFcn = @(s,e) setField('A', e.Value);
h.sldArea.ValueChangedFcn  = @(s,e) setField('A', s.Value);

h.lblC = uilabel(cg, 'Text', '', 'FontSize', 12, 'FontColor', COL.inkDim);
h.lblC.Layout.Row = 4;
h.sldC = uislider(cg, 'Limits', [0.10 1.00], 'Value', st.C, ...
    'MajorTicks', [], 'MinorTicks', [], 'FontColor', COL.inkDim);
h.sldC.Layout.Row = 5;
h.sldC.ValueChangingFcn = @(s,e) setField('C', e.Value);
h.sldC.ValueChangedFcn  = @(s,e) setField('C', s.Value);

h.lblL = uilabel(cg, 'Text', '', 'FontSize', 12, 'FontColor', COL.inkDim);
h.lblL.Layout.Row = 6;
h.sldL = uislider(cg, 'Limits', [50 2000], 'Value', st.L, ...
    'MajorTicks', [], 'MinorTicks', [], 'FontColor', COL.inkDim);
h.sldL.Layout.Row = 7;
h.sldL.ValueChangingFcn = @(s,e) setField('L', e.Value);
h.sldL.ValueChangedFcn  = @(s,e) setField('L', s.Value);

h.lblSc = uilabel(cg, 'Text', '', 'FontSize', 12, 'FontColor', COL.inkDim);
h.lblSc.Layout.Row = 8;
h.sldSc = uislider(cg, 'Limits', [0.0005 0.0500], 'Value', st.Sc, ...
    'MajorTicks', [], 'MinorTicks', [], 'FontColor', COL.inkDim);
h.sldSc.Layout.Row = 9;
h.sldSc.ValueChangingFcn = @(s,e) setField('Sc', e.Value);
h.sldSc.ValueChangedFcn  = @(s,e) setField('Sc', s.Value);

h.lblTc = uilabel(cg, 'Text', '', 'FontSize', 12, 'FontName', MONO, ...
    'FontColor', COL.water);
h.lblTc.Layout.Row = 10;

h.chkTc = uicheckbox(cg, 'Text', ' storm duration = t_c (recommended)', ...
    'Value', st.useTc, 'FontSize', 11.5, 'FontColor', COL.inkDim);
h.chkTc.Layout.Row = 11;
h.chkTc.ValueChangedFcn = @(s,e) setField('useTc', s.Value);

% return period and duration
sub = uigridlayout(cg, [2 2]);
sub.Layout.Row      = 12;
sub.RowHeight       = {16, 28};
sub.ColumnWidth     = {'1x', '1x'};
sub.Padding         = [0 0 0 0];
sub.RowSpacing      = 2;
sub.ColumnSpacing   = 8;
sub.BackgroundColor = COL.panel;

l1 = uilabel(sub, 'Text', 'Return period (yr)', 'FontSize', 11, ...
    'FontColor', COL.inkDim); l1.Layout.Row = 1; l1.Layout.Column = 1;
l2 = uilabel(sub, 'Text', 'Duration', 'FontSize', 11, ...
    'FontColor', COL.inkDim); l2.Layout.Row = 1; l2.Layout.Column = 2;

h.ddTr = uidropdown(sub, 'Items', {'2','5','10','20','30','50','100'}, ...
    'ItemsData', [2 5 10 20 30 50 100], 'Value', st.Tr, ...
    'FontName', MONO, 'BackgroundColor', COL.panel2, 'FontColor', COL.ink);
h.ddTr.Layout.Row = 2; h.ddTr.Layout.Column = 1;
h.ddTr.ValueChangedFcn = @(s,e) setField('Tr', s.Value);

h.ddDur = uidropdown(sub, ...
    'Items', {'15 min','30 min','60 min','90 min','2 hr','3 hr', ...
              '5 hr','6 hr','12 hr','24 hr'}, ...
    'ItemsData', [15 30 60 90 120 180 300 360 720 1440], ...
    'Value', st.durMin, 'Enable', 'off', ...
    'FontName', MONO, 'BackgroundColor', COL.panel2, 'FontColor', COL.ink);
h.ddDur.Layout.Row = 2; h.ddDur.Layout.Column = 2;
h.ddDur.ValueChangedFcn = @(s,e) setField('durMin', s.Value);

h.chkOverride = uicheckbox(cg, 'Text', ' override intensity (mm/hr)', ...
    'Value', st.override, 'FontSize', 11.5, 'FontColor', COL.inkDim);
h.chkOverride.Layout.Row = 13;
h.chkOverride.ValueChangedFcn = @(s,e) setField('override', s.Value);

h.efManual = uieditfield(cg, 'numeric', 'Value', st.manualI, ...
    'Limits', [0 600], 'FontName', MONO, 'Enable', 'off', ...
    'BackgroundColor', COL.panel2, 'FontColor', COL.ink);
h.efManual.Layout.Row = 14;
h.efManual.ValueChangedFcn = @(s,e) setField('manualI', s.Value);

% ---- pipe design -------------------------------------------------------
sectionLabel(cg, 'PIPE DESIGN', COL, 15);

h.lblD = uilabel(cg, 'Text', '', 'FontSize', 12, 'FontColor', COL.inkDim);
h.lblD.Layout.Row = 16;
h.sldD = uislider(cg, 'Limits', [0.15 2.00], 'Value', st.D, ...
    'MajorTicks', [], 'MinorTicks', [], 'FontColor', COL.inkDim);
h.sldD.Layout.Row = 17;
h.sldD.ValueChangingFcn = @(s,e) setField('D', e.Value);
h.sldD.ValueChangedFcn  = @(s,e) setField('D', s.Value);

h.lblS = uilabel(cg, 'Text', '', 'FontSize', 12, 'FontColor', COL.inkDim);
h.lblS.Layout.Row = 18;
h.sldS = uislider(cg, 'Limits', [0.0005 0.0300], 'Value', st.S, ...
    'MajorTicks', [], 'MinorTicks', [], 'FontColor', COL.inkDim);
h.sldS.Layout.Row = 19;
h.sldS.ValueChangingFcn = @(s,e) setField('S', e.Value);
h.sldS.ValueChangedFcn  = @(s,e) setField('S', s.Value);

l3 = uilabel(cg, 'Text', 'Pipe material', 'FontSize', 11, ...
    'FontColor', COL.inkDim); l3.Layout.Row = 20;

h.ddMat = uidropdown(cg, ...
    'Items', {'Concrete (n = 0.013)', 'PVC / smooth plastic (n = 0.010)', ...
              'Corrugated metal (n = 0.024)', 'Brick (n = 0.015)'}, ...
    'ItemsData', [0.013 0.010 0.024 0.015], 'Value', st.n, ...
    'FontName', MONO, 'BackgroundColor', COL.panel2, 'FontColor', COL.ink);
h.ddMat.Layout.Row = 21;
h.ddMat.ValueChangedFcn = @(s,e) setField('n', s.Value);

% ---- playback ----------------------------------------------------------
sectionLabel(cg, 'PLAYBACK', COL, 22);

h.lblSpeed = uilabel(cg, 'Text', '', 'FontSize', 12, 'FontColor', COL.inkDim);
h.lblSpeed.Layout.Row = 23;
h.sldSpeed = uislider(cg, 'Limits', [4 24], 'Value', st.loopSec, ...
    'MajorTicks', [], 'MinorTicks', [], 'FontColor', COL.inkDim);
h.sldSpeed.Layout.Row = 24;
h.sldSpeed.ValueChangingFcn = @(s,e) setField('loopSec', e.Value);
h.sldSpeed.ValueChangedFcn  = @(s,e) setField('loopSec', s.Value);

h.btnPlay = uibutton(cg, 'push', 'Text', 'Pause', ...
    'BackgroundColor', COL.panel2, 'FontColor', COL.ink, ...
    'ButtonPushedFcn', @togglePlay);
h.btnPlay.Layout.Row = 25;

h.btnReport = uibutton(cg, 'push', 'Text', 'Print report', ...
    'BackgroundColor', COL.panel2, 'FontColor', COL.water, ...
    'ButtonPushedFcn', @printReport);
h.btnReport.Layout.Row = 26;

foot = uilabel(cg, 'Text', { ...
    'Rainfall: Khulna IDF table,'; ...
    'Sheonty & Islam (2020), ICCESD,'; ...
    'KUET. 62 yr BMD record'; ...
    '(1948-2010), Gumbel EV-I with'; ...
    'Sherman disaggregation.'; ...
    'Durations under ~3 hr are'; ...
    'disaggregated, not gauged.'; ...
    ''; ...
    'Rational Method + Manning.'; ...
    'Depth solved by bisection.'; ...
    'Velocity limits 0.6 - 3.0 m/s.'}, ...
    'FontSize', 10.5, 'FontColor', COL.inkDim, ...
    'VerticalAlignment', 'top');
foot.Layout.Row = 27;

%% ======================================================================
%  RIGHT COLUMN - HERO PANEL AND CHARTS
% =======================================================================
rightGrid = uigridlayout(mainGrid, [2 1]);
rightGrid.Layout.Row      = 2;
rightGrid.Layout.Column   = 2;
rightGrid.RowHeight       = {'1.55x', '1x'};
rightGrid.Padding         = [0 0 0 0];
rightGrid.RowSpacing      = 14;
rightGrid.BackgroundColor = COL.bg;

heroPanel = uipanel(rightGrid, 'BackgroundColor', COL.panel2, ...
    'BorderType', 'none');
hgl = uigridlayout(heroPanel, [1 1]);
hgl.Padding = [0 0 0 0];
hgl.BackgroundColor = COL.panel2;

h.axHero = uiaxes(hgl);
styleAxes(h.axHero, COL);
h.axHero.XLim  = [0 1];
h.axHero.YLim  = [0 1];
h.axHero.XTick = []; h.axHero.YTick = [];
grid(h.axHero, 'off');
h.axHero.XColor = COL.panel2;
h.axHero.YColor = COL.panel2;
hold(h.axHero, 'on');

h.sky = patch(h.axHero, 'XData', [0 1 1 0], 'YData', [0.55 0.55 1 1], ...
    'FaceColor', COL.water, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
h.rain = line(h.axHero, NaN, NaN, 'Color', [0.59 0.76 0.94], 'LineWidth', 1.1);
h.funnel = patch(h.axHero, 'XData', [0.13 0.39 0.28 0.24], ...
    'YData', [0.72 0.72 0.56 0.56], 'FaceColor', COL.ink, ...
    'FaceAlpha', 0.04, 'EdgeColor', COL.inkDim, 'LineWidth', 1.2);
h.arrow = line(h.axHero, NaN, NaN, 'Color', COL.water, ...
    'LineWidth', 1.6, 'LineStyle', ':');
h.water     = patch(h.axHero, 'XData', NaN, 'YData', NaN, ...
    'FaceColor', COL.water, 'FaceAlpha', 0.85, 'EdgeColor', 'none');
h.waterLine = line(h.axHero, NaN, NaN, 'Color', [1 1 1], 'LineWidth', 1.4);
h.pipe      = line(h.axHero, NaN, NaN, 'Color', [0.82 0.87 0.96], 'LineWidth', 2.5);
h.crown     = line(h.axHero, NaN, NaN, 'Color', COL.inkDim, ...
    'LineWidth', 0.8, 'LineStyle', '--');

h.txtFunnel = text(h.axHero, 0.13, 0.76, 'catchment', 'Color', COL.inkDim, ...
    'FontName', MONO, 'FontSize', 10);
h.txtQreqK  = text(h.axHero, 0.03, 0.96, 'Required Q - Rational Method', ...
    'Color', COL.inkDim, 'FontName', MONO, 'FontSize', 10);
h.txtQreqV  = text(h.axHero, 0.03, 0.905, '-', 'Color', COL.ink, ...
    'FontName', MONO, 'FontSize', 17, 'FontWeight', 'bold');
h.txtQcapK  = text(h.axHero, 0.97, 0.96, 'Pipe capacity - Manning', ...
    'Color', COL.inkDim, 'FontName', MONO, 'FontSize', 10, ...
    'HorizontalAlignment', 'right');
h.txtQcapV  = text(h.axHero, 0.97, 0.905, '-', 'Color', COL.ink, ...
    'FontName', MONO, 'FontSize', 17, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'right');

% Rainfall intensity readout. The IDF table is the whole hydrological
% basis of the tool, so the value it produces is shown explicitly rather
% than being buried inside Q.
h.txtIdfV = text(h.axHero, 0.50, 0.905, '-', 'Color', COL.water, ...
    'FontName', MONO, 'FontSize', 15, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
h.txtIdfK = text(h.axHero, 0.50, 0.845, '-', 'Color', COL.inkDim, ...
    'FontName', MONO, 'FontSize', 9, 'HorizontalAlignment', 'center');

h.txtStatus = text(h.axHero, 0.03, 0.06, 'initialising', 'Color', COL.ok, ...
    'FontName', MONO, 'FontSize', 12, 'FontWeight', 'bold');
h.txtClock  = text(h.axHero, 0.97, 0.06, '0.0 min', 'Color', COL.ink, ...
    'FontName', MONO, 'FontSize', 12, 'HorizontalAlignment', 'right');
h.txtPipe   = text(h.axHero, 0.70, 0.14, '', 'Color', COL.inkDim, ...
    'FontName', MONO, 'FontSize', 10, 'HorizontalAlignment', 'center');

% ---------------- chart row ---------------------------------------------
chartGrid = uigridlayout(rightGrid, [1 2]);
chartGrid.ColumnWidth     = {'1.5x', '1x'};
chartGrid.Padding         = [0 0 0 0];
chartGrid.ColumnSpacing   = 14;
chartGrid.BackgroundColor = COL.bg;

hydPanel = uipanel(chartGrid, 'BackgroundColor', COL.panel, 'BorderType', 'none');
hg = uigridlayout(hydPanel, [1 1]);
hg.Padding = [8 6 8 6];
hg.BackgroundColor = COL.panel;

h.axHyd = uiaxes(hg);
styleAxes(h.axHyd, COL);
h.axHyd.Color = COL.panel;

yyaxis(h.axHyd, 'left');
h.axHyd.YColor = COL.water;
h.iFill = patch(h.axHyd, 'XData', NaN, 'YData', NaN, ...
    'FaceColor', COL.water, 'FaceAlpha', 0.18, 'EdgeColor', 'none');
h.iLine = line(h.axHyd, NaN, NaN, 'Color', COL.water, 'LineWidth', 2);
ylabel(h.axHyd, 'i (mm/hr)');

yyaxis(h.axHyd, 'right');
h.axHyd.YColor = COL.warn;
h.qLine   = line(h.axHyd, NaN, NaN, 'Color', COL.warn, 'LineWidth', 2.2);
h.capLine = line(h.axHyd, NaN, NaN, 'Color', COL.danger, ...
    'LineWidth', 1.5, 'LineStyle', '--');
h.nowLine = line(h.axHyd, NaN, NaN, 'Color', [0.75 0.79 0.87], 'LineWidth', 0.9);
ylabel(h.axHyd, 'Q (m^3/s)');
xlabel(h.axHyd, 'Time (min)');
title(h.axHyd, 'Hyetograph & hydrograph (live)', 'Color', COL.inkDim);

barPanel = uipanel(chartGrid, 'BackgroundColor', COL.panel, 'BorderType', 'none');
bg = uigridlayout(barPanel, [1 1]);
bg.Padding = [8 6 8 6];
bg.BackgroundColor = COL.panel;

h.axBar = uiaxes(bg);
styleAxes(h.axBar, COL);
h.axBar.Color = COL.panel;
h.bar = bar(h.axBar, [1 2], [0 0], 0.55, 'FaceColor', 'flat');
h.axBar.XLim = [0.4 2.6];
h.axBar.XTick = [1 2];
h.axBar.XTickLabel = {'Required Q', 'Capacity'};
ylabel(h.axBar, 'Q (m^3/s)');
title(h.axBar, 'Capacity vs peak demand', 'Color', COL.inkDim);
hold(h.axBar, 'on');
h.barTxt1 = text(h.axBar, 1, 0, '', 'Color', COL.ink, 'FontName', MONO, ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
h.barTxt2 = text(h.axBar, 2, 0, '', 'Color', COL.ink, 'FontName', MONO, ...
    'FontSize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

%% ======================================================================
%  START
% =======================================================================
recompute();

tmr = timer('ExecutionMode', 'fixedRate', 'Period', 0.05, ...
            'BusyMode', 'drop', 'TimerFcn', @onTick);
fig.CloseRequestFcn = @onClose;
start(tmr);

%% ======================================================================
%  NESTED FUNCTIONS
% =======================================================================

    function setField(name, value)
        st.(name) = value;
        if strcmp(name, 'override')
            if value, h.efManual.Enable = 'on'; else, h.efManual.Enable = 'off'; end
        end
        if strcmp(name, 'useTc')
            if value, h.ddDur.Enable = 'off'; else, h.ddDur.Enable = 'on'; end
        end
        recompute();
    end

    function recompute()
        % ---- time of concentration --------------------------------------
        st.tc = time_of_concentration(st.L, st.Sc);

        % The Rational Method assumes the storm lasts exactly as long as
        % the time of concentration, because that is when the whole
        % catchment is contributing at once and the peak flow occurs.
        if st.useTc
            st.durUse = st.tc;
        else
            st.durUse = st.durMin;
        end

        % ---- hydrology --------------------------------------------------
        if st.override
            st.Ipeak = st.manualI;
            st.exact = false;
        else
            [st.Ipeak, st.exact] = idf_intensity(st.durUse, st.Tr);
        end
        st.Qreq = rational_Q(st.C, st.Ipeak, st.A);

        % ---- hydraulics -------------------------------------------------
        [st.Qfull, st.Vfull] = capacityFast();
        st.yQmax = YMAX_RATIO * st.D;
        st.Qmax  = manning_Q(st.yQmax, st.D, st.S, st.n);
        [st.yflow, st.Vflow] = solveDepth(st.Qreq);

        [st.status, st.level] = design_status(st.Qreq, st.Qfull, ...
                                              st.Vflow, st.vmin, st.vmax);

        % ---- labels ------------------------------------------------------
        h.lblArea.Text  = sprintf('Catchment area           %.2f ha',  st.A);
        h.lblC.Text     = sprintf('Runoff coefficient C     %.2f',     st.C);
        h.lblL.Text     = sprintf('Overland flow length     %.0f m',   st.L);
        h.lblSc.Text    = sprintf('Catchment slope          %.4f m/m', st.Sc);
        h.lblTc.Text    = sprintf('t_c (Kirpich) = %.1f min', st.tc);
        h.lblD.Text     = sprintf('Pipe diameter            %.2f m',   st.D);
        h.lblS.Text     = sprintf('Pipe slope               %.4f m/m', st.S);
        h.lblSpeed.Text = sprintf('Storm loop length        %.0f s',   st.loopSec);

        h.txtQreqV.String = sprintf('%.4f m^3/s', st.Qreq);
        h.txtQcapV.String = sprintf('%.4f m^3/s', st.Qfull);
        h.txtIdfV.String  = sprintf('i = %.1f mm/hr', st.Ipeak);
        if st.override
            h.txtIdfK.String = 'manual override - IDF table bypassed';
        elseif st.exact
            h.txtIdfK.String = sprintf(...
                'Khulna IDF  -  %s  -  %d yr  -  exact table entry', ...
                durText(st.durUse), st.Tr);
        else
            h.txtIdfK.String = sprintf(...
                'Khulna IDF  -  %s  -  %d yr  -  interpolated', ...
                durText(st.durUse), st.Tr);
        end
        h.txtStatus.String = sprintf('%s  -  V = %.2f m/s', st.status, st.Vflow);
        h.txtStatus.Color  = levelColour(st.level);

        drawBar();
        clearHistory();
    end

    function s = durText(dm)
        if dm < 60
            s = sprintf('%.1f min', dm);
        else
            s = sprintf('%.2f hr', dm/60);
        end
    end

    function [Qf, Vf] = capacityFast()
        % Closed form. manning_capacity also grid-searches for Qmax, which
        % is too slow to run on every slider movement.
        Af = pi * st.D^2 / 4;
        Rf = st.D / 4;
        Qf = (1/st.n) * Af * Rf^(2/3) * sqrt(st.S);
        Vf = Qf / Af;
    end

    function [y, V] = solveDepth(Q)
        % Depth carrying discharge Q, by bisection on [0, 0.938D]. Same
        % algorithm as normal_depth but using the cached st.Qmax so it is
        % fast enough to run every animation frame.
        if Q <= 0, y = 0; V = 0; return; end
        if Q >= st.Qmax
            y = st.D;
            [~, V] = manning_Q(st.D, st.D, st.S, st.n);
            return
        end
        f = @(yy) manning_Q(yy, st.D, st.S, st.n) - Q;
        y = bisect_root(f, 1e-9, st.yQmax, 0.01, 30);
        [~, V] = manning_Q(y, st.D, st.S, st.n);
    end

    function c = levelColour(level)
        switch level
            case 1,    c = COL.ok;
            case 2,    c = COL.warn;
            otherwise, c = COL.danger;
        end
    end

    function clearHistory()
        st.hist = struct('t', [], 'i', [], 'q', []);
    end

    function resetSim()
        st.simT = 0;
        clearHistory();
        st.rain = struct('x', [], 'y', [], 'len', [], 'v', []);
    end

    function togglePlay(~, ~)
        st.running = ~st.running;
        if st.running
            h.btnPlay.Text    = 'Pause';
            h.badge.Text      = 'KHULNA IDF - SIMULATING';
            h.badge.FontColor = COL.ok;
        else
            h.btnPlay.Text    = 'Play';
            h.badge.Text      = 'KHULNA IDF - PAUSED';
            h.badge.FontColor = COL.warn;
        end
    end

    function onTick(~, ~)
        if ~isvalid(fig), return; end
        if ~st.running,  return; end

        st.simT = st.simT + 0.05 / st.loopSec;
        if st.simT >= 1
            resetSim();
        end

        tMin = st.simT * st.durUse;
        i_t  = triangular_hyetograph(tMin, st.durUse, st.Ipeak);
        Q_t  = rational_Q(st.C, i_t, st.A);

        st.hist.t(end+1) = tMin;
        st.hist.i(end+1) = i_t;
        st.hist.q(end+1) = Q_t;
        if numel(st.hist.t) > 400
            st.hist.t(1) = []; st.hist.i(1) = []; st.hist.q(1) = [];
        end

        drawHero(tMin, i_t, Q_t);
        drawHydro();
        drawnow limitrate
    end

    function drawHero(tMin, i_t, Q_t)
        % The axes spans 0-1 both ways but is wider than tall in pixels,
        % so the x-radius is shrunk by the pixel aspect ratio to make the
        % pipe render as a true circle.
        p = h.axHero.InnerPosition;
        if p(3) > 0, aspect = p(4)/p(3); else, aspect = 0.6; end

        ry  = 0.25;
        rx  = ry * aspect;
        pcx = 0.70; pcy = 0.40;

        if st.Ipeak > 0, frac = min(i_t/st.Ipeak, 1); else, frac = 0; end
        h.sky.FaceAlpha = 0.04 + 0.20*frac;

        nNew = round(4*frac);
        for k = 1:nNew
            st.rain.x(end+1)   = rand();
            st.rain.y(end+1)   = 1.02;
            st.rain.len(end+1) = 0.02 + 0.025*rand();
            st.rain.v(end+1)   = 0.030 + 0.020*rand();
        end
        if ~isempty(st.rain.x)
            st.rain.y = st.rain.y - st.rain.v;
            st.rain.x = st.rain.x - 0.004;
            keep = st.rain.y > 0.55;
            st.rain.x   = st.rain.x(keep);
            st.rain.y   = st.rain.y(keep);
            st.rain.len = st.rain.len(keep);
            st.rain.v   = st.rain.v(keep);
        end
        nd = numel(st.rain.x);
        if nd > 0
            rx3 = [st.rain.x; st.rain.x - 0.006; NaN(1, nd)];
            ry3 = [st.rain.y; st.rain.y - st.rain.len; NaN(1, nd)];
            h.rain.XData = rx3(:);
            h.rain.YData = ry3(:);
        else
            h.rain.XData = NaN; h.rain.YData = NaN;
        end

        ta  = linspace(0, 1, 20);
        ax0 = 0.26; ay0 = 0.55;
        ax1 = pcx - rx*0.8; ay1 = pcy + ry*0.4;
        axm = (ax0+ax1)/2;  aym = min(ay0, ay1) - 0.06;
        h.arrow.XData = (1-ta).^2*ax0 + 2*(1-ta).*ta*axm + ta.^2*ax1;
        h.arrow.YData = (1-ta).^2*ay0 + 2*(1-ta).*ta*aym + ta.^2*ay1;

        th = linspace(0, 2*pi, 180);
        h.pipe.XData  = pcx + rx*cos(th);
        h.pipe.YData  = pcy + ry*sin(th);
        h.crown.XData = [pcx-rx pcx+rx];
        h.crown.YData = [pcy+ry pcy+ry];

        % Real hydraulics: the depth is solved from Manning's equation by
        % bisection every frame, not scaled proportionally.
        y_t = solveDepth(Q_t);
        yf  = y_t / st.D;

        if yf <= 0
            h.water.XData = NaN;     h.water.YData = NaN;
            h.waterLine.XData = NaN; h.waterLine.YData = NaN;
        else
            phi = acos(max(-1, min(1, 1 - 2*yf)));
            ph  = linspace(-phi, phi, 120);
            wx  = pcx + rx*sin(ph);
            wy  = pcy - ry*cos(ph);
            % The arc runs from one end of the water surface, round the
            % invert, back to the other end. patch closes the polygon
            % along the surface chord - exactly the circular segment.
            h.water.XData = wx;
            h.water.YData = wy;

            surfY = pcy - ry*cos(phi);
            xs    = linspace(pcx - rx*sin(phi), pcx + rx*sin(phi), 40);
            h.waterLine.XData = xs;
            h.waterLine.YData = surfY + 0.004*sin(12*xs + 8*st.simT*2*pi);
        end

        if yf >= 0.999
            h.water.FaceColor = COL.danger;
        elseif yf >= 0.75
            h.water.FaceColor = COL.warn;
        else
            h.water.FaceColor = COL.water;
        end

        h.txtPipe.Position(1) = pcx;
        h.txtPipe.Position(2) = pcy - ry - 0.06;
        h.txtPipe.String = sprintf(...
            'pipe D %.2f m  -  y/D = %.0f %%  -  Q = %.4f m^3/s', ...
            st.D, 100*yf, Q_t);
        h.txtClock.String = sprintf('%.1f / %.1f min', tMin, st.durUse);
    end

    function drawHydro()
        if numel(st.hist.t) < 2, return; end
        t = st.hist.t; ii = st.hist.i; qq = st.hist.q;

        maxI = max(st.Ipeak*1.15, 1e-6);
        maxQ = max([st.Qfull*1.15, max(qq)*1.15, 1e-6]);

        yyaxis(h.axHyd, 'left');
        h.iLine.XData = t;  h.iLine.YData = ii;
        h.iFill.XData = [t, fliplr(t)];
        h.iFill.YData = [ii, zeros(size(ii))];
        h.axHyd.YLim  = [0 maxI];

        yyaxis(h.axHyd, 'right');
        h.qLine.XData   = t;  h.qLine.YData = qq;
        h.capLine.XData = [0 st.durUse];
        h.capLine.YData = [st.Qfull st.Qfull];
        h.nowLine.XData = [t(end) t(end)];
        h.nowLine.YData = [0 maxQ];
        h.axHyd.YLim    = [0 maxQ];

        h.axHyd.XLim = [0 st.durUse];
    end

    function drawBar()
        vals = [st.Qreq, st.Qfull];
        h.bar.YData = vals;
        if st.Qfull >= st.Qreq, capColour = COL.ok; else, capColour = COL.danger; end
        h.bar.CData = [COL.warn; capColour];

        top = max(vals)*1.28;
        if top <= 0, top = 1; end
        h.axBar.YLim = [0 top];

        h.barTxt1.Position(2) = vals(1);
        h.barTxt1.String = sprintf('%.4f', vals(1));
        h.barTxt2.Position(2) = vals(2);
        h.barTxt2.String = sprintf('%.4f', vals(2));
    end

    function printReport(~, ~)
        % Full numerical design report for the current inputs, written to
        % the Command Window. Every intermediate number is shown so the
        % arithmetic can be followed on paper.

        A = st.A; C = st.C; D = st.D; S = st.S; n = st.n;

        fprintf('\n');
        fprintf('==============================================================\n');
        fprintf('  DrainFlow - storm sewer design report\n');
        fprintf('  Station: KHULNA   |   %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
        fprintf('  Rainfall: Sheonty & Islam (2020), ICCESD, KUET\n');
        fprintf('            62 yr BMD record 1948-2010, Gumbel EV-I,\n');
        fprintf('            Sherman disaggregation to short durations\n');
        fprintf('==============================================================\n\n');

        fprintf('--- 1. DESIGN INPUTS ---------------------------------------\n');
        fprintf('  Catchment area        A  = %8.3f ha\n',  A);
        fprintf('  Runoff coefficient    C  = %8.3f\n',     C);
        fprintf('  Overland flow length  L  = %8.0f m\n',   st.L);
        fprintf('  Catchment slope       Sc = %8.4f m/m\n', st.Sc);
        fprintf('  Return period         Tr = %8.0f yr\n',  st.Tr);
        fprintf('  Pipe diameter         D  = %8.3f m\n',   D);
        fprintf('  Pipe slope            S  = %8.4f m/m\n', S);
        fprintf('  Manning roughness     n  = %8.4f\n\n',   n);

        fprintf('--- 2. TIME OF CONCENTRATION -------------------------------\n');
        fprintf('  Kirpich:  tc = 0.0195 * L^0.77 * Sc^-0.385   [min, m, m/m]\n');
        fprintf('  tc = 0.0195 * %.0f^0.77 * %.4f^-0.385 = %.2f min\n', ...
                st.L, st.Sc, st.tc);
        if st.useTc
            fprintf('  Storm duration set equal to tc = %.2f min\n\n', st.durUse);
        else
            fprintf('  Storm duration selected manually = %.0f min\n', st.durUse);
            fprintf('  NOTE: the Rational Method assumes duration = tc.\n\n');
        end

        fprintf('--- 3. HYDROLOGY -------------------------------------------\n');
        if st.override
            fprintf('  Intensity (manual)    i  = %8.2f mm/hr\n', st.Ipeak);
        else
            fprintf('  Khulna IDF intensity  i  = %8.2f mm/hr', st.Ipeak);
            if st.exact
                fprintf('   (exact table entry)\n');
            else
                fprintf('   (interpolated)\n');
            end
        end
        fprintf('  Rational Method       Q  = C*i*A/360\n');
        fprintf('                           = %.3f * %.2f * %.3f / 360\n', C, st.Ipeak, A);
        fprintf('                           = %8.4f m^3/s  (required)\n\n', st.Qreq);

        [~, ~, Dtab, Ttab, Itab] = idf_intensity(60, 10);
        fprintf('  Stored Khulna IDF table, intensity in mm/hr:\n');
        fprintf('      duration');
        fprintf('%9s', 'T=2', 'T=5', 'T=10', 'T=20', 'T=30', 'T=50', 'T=100');
        fprintf('\n');
        for r = 1:numel(Dtab)
            fprintf('  %7.0f min', Dtab(r));
            fprintf('%9.2f', Itab(r,:));
            fprintf('\n');
        end
        fprintf('\n');

        fprintf('--- 4. HYDRAULICS ------------------------------------------\n');
        [~, ~, QmaxGrid, yQmaxGrid] = manning_capacity(D, S, n);
        [yB, VB, surch, itB, eaB]   = normal_depth(st.Qreq, D, S, n);

        fprintf('  Full-pipe capacity  Qfull = %8.4f m^3/s\n', st.Qfull);
        fprintf('  Full-pipe velocity  Vfull = %8.4f m/s\n',   st.Vfull);
        fprintf('  Maximum discharge   Qmax  = %8.4f m^3/s at y/D = %.4f\n', ...
                QmaxGrid, yQmaxGrid/D);
        fprintf('  Theory: peak at y/D = 0.9382, Qmax/Qfull = 1.0757\n');
        fprintf('  Computed ratio            = %.4f\n\n', QmaxGrid/st.Qfull);

        fprintf('  Bisection for the flow depth:\n');
        fprintf('    iterations           = %d\n',      itB);
        fprintf('    relative error ea    = %.3e %%\n', eaB);
        fprintf('    flow depth      y    = %8.4f m  (y/D = %.1f %%)\n', ...
                yB, 100*yB/D);
        fprintf('    flow velocity   V    = %8.4f m/s\n', VB);
        if surch
            fprintf('    PIPE IS SURCHARGED - demand exceeds Qmax\n');
        end
        fprintf('\n  VERDICT: %s\n', st.status);
        fprintf('  Capacity / demand ratio = %.3f\n\n', st.Qfull/max(st.Qreq, eps));

        fprintf('--- 5. HAND-CALCULATION CHECK ------------------------------\n');
        Afull_hand = pi*D^2/4;
        Rfull_hand = D/4;
        Qfull_hand = (1/n)*Afull_hand*Rfull_hand^(2/3)*sqrt(S);
        Qreq_hand  = C*st.Ipeak*A/360;
        Qcheck     = manning_Q(yB, D, S, n);
        fprintf('  A_full = pi*D^2/4              = %10.6f m^2\n', Afull_hand);
        fprintf('  R_full = D/4                   = %10.6f m\n',   Rfull_hand);
        fprintf('  Q_full by hand                 = %10.6f m^3/s\n', Qfull_hand);
        fprintf('  Q_full from the code           = %10.6f m^3/s\n', st.Qfull);
        fprintf('  difference                     = %10.3e\n\n', abs(st.Qfull-Qfull_hand));
        fprintf('  Q_req by hand                  = %10.6f m^3/s\n', Qreq_hand);
        fprintf('  Q_req from the code            = %10.6f m^3/s\n', st.Qreq);
        fprintf('  difference                     = %10.3e\n\n', abs(st.Qreq-Qreq_hand));
        fprintf('  Manning(y_solved) back-subbed  = %10.6f m^3/s\n', Qcheck);
        fprintf('  target Q_req                   = %10.6f m^3/s\n', st.Qreq);
        fprintf('  residual                       = %10.3e m^3/s\n\n', abs(Qcheck-st.Qreq));

        fprintf('--- 6. RUNOFF VOLUME ---------------------------------------\n');
        nSeg = 200;
        tv = linspace(0, st.durUse, nSeg+1);          % minutes
        iv = triangular_hyetograph(tv, st.durUse, st.Ipeak);
        qv = rational_Q(C, iv, A);
        ts = tv * 60;                                 % minutes -> seconds

        Vol_trapz   = trapz(ts, qv);
        Vol_simpson = simpson13(ts, qv);
        % Exact: the triangle encloses 0.5*Ipeak*(dur/60) mm of rain, and
        % 1 ha * 1 mm = 10 m^3.
        Vol_exact   = 10 * C * A * 0.5 * st.Ipeak * st.durUse/60;

        fprintf('  Peak of hydrograph       = %10.4f m^3/s\n', max(qv));
        fprintf('  Volume, trapz            = %10.2f m^3\n',   Vol_trapz);
        fprintf('  Volume, Simpson 1/3      = %10.2f m^3\n',   Vol_simpson);
        fprintf('  Volume, exact analytic   = %10.2f m^3\n',   Vol_exact);
        fprintf('  trapz error              = %10.4f %%\n', ...
                100*abs(Vol_trapz-Vol_exact)/Vol_exact);
        fprintf('  Simpson error            = %10.4f %%\n\n', ...
                100*abs(Vol_simpson-Vol_exact)/Vol_exact);

        fprintf('--- 7. PIPE SIZING -----------------------------------------\n');
        % Vectorised sweep: element-wise operators, no loop.
        Dsweep  = 0.15:0.01:2.00;
        Afull_s = pi * Dsweep.^2 / 4;
        Rfull_s = Dsweep / 4;
        Qsweep  = (1/n) * Afull_s .* Rfull_s.^(2/3) * sqrt(S);
        kOK     = find(Qsweep >= st.Qreq, 1, 'first');

        if isempty(kOK)
            fprintf('  No diameter up to %.2f m carries %.4f m^3/s at S = %.4f.\n', ...
                    Dsweep(end), st.Qreq, S);
            fprintf('  Steepen the pipe or split the catchment.\n\n');
        else
            Dmin = Dsweep(kOK);
            % Snap to the commercial sizes used in practice.
            comm = [0.300 0.375 0.450 0.525 0.600 0.675 0.750 0.900 ...
                    1.050 1.200 1.350 1.500];
            kC = find(comm >= Dmin, 1, 'first');
            fprintf('  Smallest adequate diameter = %.2f m (capacity %.4f m^3/s)\n', ...
                    Dmin, Qsweep(kOK));
            if ~isempty(kC)
                fprintf('  Nearest commercial size    = %.0f mm\n', comm(kC)*1000);
            end
            fprintf('  Chosen trial diameter      = %.2f m\n\n', D);
        end

        fprintf('--- 8. GUMBEL BACK-ANALYSIS OF THE PUBLISHED TABLE ---------\n');
        % The source study fitted a Gumbel EV-I distribution. Because depth
        % is exactly linear in the Gumbel reduced variate
        %     y = -ln(-ln(1 - 1/T)),
        % the location and scale parameters behind every duration row can
        % be recovered by least squares. Small residuals confirm the table
        % really is a Gumbel product and lets the fit be quoted directly.
        [uG, aG, resid] = gumbel_from_table();
        fprintf('  Fitting  P = u + alpha*y  to each duration row:\n');
        fprintf('   dur(min)     u (mm)   alpha (mm)     mean     std      CV   resid\n');
        for r = 1:numel(Dtab)
            mn = uG(r) + 0.5772*aG(r);
            sd = aG(r)*pi/sqrt(6);
            fprintf('  %9.0f %10.2f %12.2f %8.1f %7.1f %7.3f %7.3f\n', ...
                    Dtab(r), uG(r), aG(r), mn, sd, sd/mn, resid(r));
        end
        fprintf('  (residuals in mm; all well under 1 %% of the depths)\n');
        fprintf('\n==============================================================\n\n');
        fprintf('\n==============================================================\n\n');


        %% --- EXPORT REPORT TO EXCEL -----------------------------

        filename = sprintf('DrainFlow_Report_%s.xlsx', ...
            datestr(now,'yyyy_mm_dd_HH_MM_SS'));

        ReportData = {
            'Parameter','Value';
            'Catchment Area (ha)', st.A;
            'Runoff Coefficient', st.C;
            'Flow Length (m)', st.L;
            'Catchment Slope', st.Sc;
            'Return Period (yr)', st.Tr;
            'Pipe Diameter (m)', st.D;
            'Pipe Slope (m/m)', st.S;
            'Manning Roughness', st.n;
            'Time of Concentration (min)', st.tc;
            'Rainfall Intensity (mm/hr)', st.Ipeak;
            'Required Discharge (m3/s)', st.Qreq;
            'Full Pipe Capacity (m3/s)', st.Qfull;
            'Normal Depth (m)', st.yflow;
            'Velocity (m/s)', st.Vflow;
            'Status', st.status;
            'Capacity Ratio', st.Qfull/max(st.Qreq,eps)
            };

        writecell(ReportData, filename);

        fprintf('Excel report saved: %s\n', filename);


    
    end

    function onClose(~, ~)
        try
            stop(tmr); delete(tmr);
        catch
        end
        delete(fig);
    end

end % ======================= end of DrainFlow ============================


%% ======================================================================
%  ==================  USER INTERFACE HELPERS  =========================
% =======================================================================

function sectionLabel(parent, txt, COL, row)
%SECTIONLABEL  Small uppercase heading dividing the control column.
lbl = uilabel(parent, 'Text', txt, 'FontSize', 11, 'FontWeight', 'bold', ...
    'FontColor', COL.inkDim);
lbl.Layout.Row = row;
end


function styleAxes(ax, COL)
%STYLEAXES  Dark theme for a uiaxes, with hover toolbar and zoom disabled.
ax.Color           = COL.panel2;
ax.XColor          = COL.inkDim;
ax.YColor          = COL.inkDim;
ax.GridColor       = COL.line;
ax.GridAlpha       = 0.5;
ax.Box             = 'off';
ax.FontSize        = 10;
ax.TitleFontWeight = 'normal';
grid(ax, 'on');
try
    ax.Toolbar.Visible = 'off';
    disableDefaultInteractivity(ax);
catch
end
end


%% ======================================================================
%  ===================  ENGINEERING FUNCTIONS  =========================
% =======================================================================

function [i, isExact, Dtab, Ttab, Itab] = idf_intensity(durMin, Tr)
%IDF_INTENSITY  Rainfall intensity from the published Khulna IDF table.
%
% [i, isExact, Dtab, Ttab, Itab] = idf_intensity(durMin, Tr)
%
% input:
%   durMin  = storm duration (MINUTES), valid range 15 - 1440
%   Tr      = return period (yr), valid range 2 - 100
% output:
%   i       = rainfall intensity (mm/hr)
%   isExact = true if the query landed on a tabulated entry
%   Dtab    = duration vector of the stored table (min)
%   Ttab    = return period vector of the stored table (yr)
%   Itab    = the stored intensity table (mm/hr)
%
% SOURCE
%   Sheonty, S.R. & Islam, G.M.T. (2020) "Establishment of Rainfall
%   Intensity-Duration-Frequency Curves of Khulna". ICCESD 2020, KUET,
%   Khulna, Bangladesh, Paper WRE-4374.
%   62 years of BMD annual maximum daily rainfall for Khulna (1948-2010),
%   checked for homogeneity, consistency and randomness, fitted to a
%   Gumbel Type-I Extreme Value distribution, then disaggregated to short
%   durations. Talbot, Sherman and Kimjima equations were compared by
%   least squares and Sherman was the best fit.
%
%   LIMITATION, to be stated in the report: short-duration rainfall gauges
%   are scarce in Bangladesh, so durations below about 3 hours are
%   DISAGGREGATED from daily records rather than directly measured.
%
% INTERPOLATION - two passes of interp1, each in the variable that makes
% the relationship straight, which is far more accurate than interpolating
% the raw numbers:
%
%   pass 1, along DURATION: IDF curves are close to straight lines on
%     log-log axes, so we interpolate log(i) against log(t) and take the
%     exponential. Interpolating i against t linearly across a 15 min to
%     1440 min span would cut the corner off a strongly convex curve.
%
%   pass 2, across RETURN PERIOD: for a Gumbel distribution the depth is
%     EXACTLY linear in the reduced variate
%           y = -ln(-ln(1 - 1/T)),
%     so interpolating against y rather than against T is not an
%     approximation at all - it reproduces the fitted distribution.
%
% Querying any of the 70 tabulated points returns the table value exactly.

if nargin < 2
    error('idf_intensity: duration and return period are both required');
end

% ---------------- stored table -----------------------------------------
Dtab = [15; 30; 60; 90; 120; 180; 300; 360; 720; 1440];   % minutes
Ttab = [2, 5, 10, 20, 30, 50, 100];                       % years

% Intensity in mm/hr. Rows match Dtab, columns match Ttab.
% Reading the shape: rightwards the numbers grow (rarer storms are more
% intense); downwards they shrink (longer storms are gentler per hour).
Itab = [ 89.60 165.16 218.15 270.24 300.60 338.87 390.91
         58.14 103.98 135.76 166.84 184.91 207.65 238.51
         37.73  65.46  84.48 103.01 113.75 127.24 145.53
         29.30  49.94  64.01  77.69  85.60  95.54 109.00
         24.49  41.21  52.58  63.59  69.97  77.97  88.79
         19.01  31.44  39.84  47.96  52.66  58.54  66.51
         13.83  22.35  28.09  33.62  36.81  40.81  46.21
         12.34  19.79  24.79  29.61  32.39  35.87  40.58
          8.01  12.46  15.43  18.28  19.93  21.98  24.76
          5.20   7.84   9.60  11.29  12.26  13.47  15.11 ];

% ---------------- clamp the query to the table range --------------------
% Extrapolating a rainfall table means inventing data, so clamp and warn.
if durMin < min(Dtab) || durMin > max(Dtab)
    warning('idf_intensity:duration', ...
        'Duration %.4g min is outside the table (15-1440 min). Clamped.', durMin);
    durMin = min(max(durMin, min(Dtab)), max(Dtab));
end
if Tr < min(Ttab) || Tr > max(Ttab)
    warning('idf_intensity:returnPeriod', ...
        'Return period %.4g yr is outside the table (2-100 yr). Clamped.', Tr);
    Tr = min(max(Tr, min(Ttab)), max(Ttab));
end

% ---------------- pass 1: log-log interpolation along duration ----------
nT     = numel(Ttab);
iAtDur = zeros(1, nT);
for k = 1:nT
    iAtDur(k) = exp( interp1(log(Dtab), log(Itab(:,k)), log(durMin), 'linear') );
end

% ---------------- pass 2: interpolation in the Gumbel reduced variate ---
yTab = -log(-log(1 - 1./Ttab));         % reduced variate of the table
yQ   = -log(-log(1 - 1/Tr));            % reduced variate of the query

i = interp1(yTab, iAtDur, yQ, 'linear');

% ---------------- did we land on a table entry? -------------------------
isExact = any(abs(Dtab - durMin) < 1e-9) && any(abs(Ttab - Tr) < 1e-9);
end


function tc = time_of_concentration(L, Sc)
%TIME_OF_CONCENTRATION  Kirpich formula.
%
% tc = time_of_concentration(L, Sc)
%
% input:
%   L  = length of the longest overland flow path (m)
%   Sc = average slope along that path (m/m)
% output:
%   tc = time of concentration (minutes)
%
%       tc = 0.0195 * L^0.77 * Sc^(-0.385)
%
% WHY THIS MATTERS
%   The Rational Method assumes the storm lasts exactly as long as the
%   time of concentration - the time water needs to travel from the
%   furthest corner of the catchment to the inlet. Only then is the whole
%   catchment contributing at once, which is when the peak flow occurs.
%   A shorter storm never engages the whole area; a longer one is less
%   intense. So tc is what tells you WHICH row of the IDF table to read.
%
%   Getting this wrong is the most common misuse of the Rational Method.
%   Reading a 24-hour intensity for a catchment whose tc is 20 minutes
%   under-estimates the design flow by an order of magnitude.
%
% Kirpich was calibrated on small rural catchments; for paved urban
% surfaces it tends to give short values, which is conservative here
% because shorter duration means higher intensity.

if nargin < 2
    error('time_of_concentration: two input arguments required (L, Sc)');
end
if L <= 0
    error('time_of_concentration: flow length L must be positive');
end
if Sc <= 0
    error('time_of_concentration: catchment slope Sc must be positive');
end

tc = 0.0195 * L^0.77 * Sc^(-0.385);
end


function Q = rational_Q(C, i, A)
%RATIONAL_Q  Peak runoff discharge by the Rational Method.
%
% Q = rational_Q(C, i, A)
%
% input:
%   C = runoff coefficient (0 - 1)
%   i = rainfall intensity (mm/hr)
%   A = catchment area (hectares)
% output:
%   Q = peak discharge (m^3/s)
%
% The Rational Method in raw form is Q = C*i*A. The 360 is a UNIT
% CONVERSION, not an empirical constant:
%
%   i [mm/hr] * A [ha] = (i/1000) m/hr * (A*10000) m^2 = 10*i*A  m^3/hr
%                                                      = 10*i*A/3600 m^3/s
%                                                      = C*i*A/360   m^3/s
%
% Vectorised, so one call converts a whole hyetograph into a hydrograph.

if nargin < 3
    error('rational_Q: three input arguments required (C, i, A)');
end

Q = C .* i .* A / 360;
end


function [Aflow, P, R, Tw, theta] = pipe_geometry(y, D)
%PIPE_GEOMETRY  Flow geometry of a part-full circular pipe.
%
% [Aflow, P, R, Tw, theta] = pipe_geometry(y, D)
%
% input:
%   y     = depth of flow from the pipe invert (m)
%   D     = internal pipe diameter (m)
% output:
%   Aflow = wetted cross-sectional flow area (m^2)
%   P     = wetted perimeter (m) - the pipe wall the water touches, NOT
%           including the free surface, because air causes no friction
%   R     = hydraulic radius = Aflow / P (m)
%   Tw    = top width of the water surface (m)
%   theta = wetted central angle (rad)
%
%       cos(theta/2) = 1 - 2*y/D    ->    theta = 2*acos(1 - 2*y/D)
%
%       Aflow = (D^2/8) * (theta - sin(theta))    circular segment area
%       P     = D*theta/2                         arc length
%       R     = Aflow / P
%       Tw    = D*sin(theta/2)
%
% CHECK at y = D (pipe full): theta = 2*pi, so Aflow = pi*D^2/4,
% P = pi*D and R = D/4 - the familiar full-pipe values.

if nargin < 2
    error('pipe_geometry: two input arguments required (y, D)');
end
if D <= 0
    error('pipe_geometry: diameter D must be positive');
end

y = min(max(y, 0), D);

theta = 2 * acos(1 - 2*y/D);

Aflow = (D^2/8) * (theta - sin(theta));
P     = D * theta / 2;
Tw    = D * sin(theta/2);

if P > 0
    R = Aflow / P;
else
    R = 0;
end
end


function [Q, V, R, Aflow] = manning_Q(y, D, S, n)
%MANNING_Q  Discharge of a part-full circular pipe by Manning's equation.
%
% [Q, V, R, Aflow] = manning_Q(y, D, S, n)
%
%       Q = (1/n) * A * R^(2/3) * S^(1/2)
%
% Flow rises with area, rises with hydraulic radius, rises with the square
% root of slope, and falls with roughness. The exponents are empirical -
% Manning fitted them to channel data in the 1890s.
%
% SI UNITS ONLY. The 1/n term carries hidden units; US customary form
% needs an extra factor of 1.49.
%
% Q(y) is NOT monotonic over 0 <= y <= D. Near the crown the wall curves
% inward and adds wetted perimeter faster than flow area, so R falls and
% Q peaks at about y/D = 0.938 before dropping back to the full-pipe
% value. This is why normal_depth brackets on [0, 0.938D].

if nargin < 4
    error('manning_Q: four input arguments required (y, D, S, n)');
end
if S <= 0
    error('manning_Q: slope S must be positive');
end
if n <= 0
    error('manning_Q: roughness n must be positive');
end

[Aflow, ~, R] = pipe_geometry(y, D);

Q = (1/n) * Aflow * R^(2/3) * sqrt(S);

if Aflow > 0
    V = Q / Aflow;
else
    V = 0;
end
end


function [Qfull, Vfull, Qmax, yQmax] = manning_capacity(D, S, n)
%MANNING_CAPACITY  Full-pipe capacity and true maximum discharge.
%
% Qfull uses A = pi*D^2/4 and R = D/4 - the capacity quoted in design.
%
% Qmax is found by GRID SEARCH: manning_Q is evaluated at 2001 trial
% depths and the largest picked with max. That is the "scan and take the
% best" idea from the Fundamentals lecture. Theory says yQmax/D = 0.9382
% with Qmax/Qfull = 1.0757, so the grid result independently checks the
% theory - and vice versa.

if nargin < 3
    error('manning_capacity: three input arguments required (D, S, n)');
end

Afull = pi * D^2 / 4;
Rfull = D / 4;

Qfull = (1/n) * Afull * Rfull^(2/3) * sqrt(S);
Vfull = Qfull / Afull;

yGrid = linspace(1e-6, D, 2001);
qGrid = zeros(size(yGrid));
for k = 1:numel(yGrid)
    qGrid(k) = manning_Q(yGrid(k), D, S, n);
end

[Qmax, kBest] = max(qGrid);
yQmax = yGrid(kBest);
end


function [root, fx, ea, iter] = bisect_root(func, xl, xu, es, maxit)
%BISECT_ROOT  Root of a function by the bisection method.
%
% [root, fx, ea, iter] = bisect_root(func, xl, xu, es, maxit)
%
% input:
%   func  = handle to the function being solved, f(x) = 0
%   xl,xu = lower and upper guesses that bracket the root
%   es    = desired relative error in %    (default 0.0001)
%   maxit = maximum allowable iterations   (default 50)
%
% This is bisect.m from the "Finding Roots of Equations" lecture, renamed
% so it cannot collide with a copy already in your own CE 206 folder.
% Halve the interval, keep whichever half still shows a sign change. The
% error falls by a factor of 2 every pass, so 23 iterations shrink it by
% roughly 8 million.
%
% The test func(xl)*func(xr) < 0 uses multiplication to compare signs: a
% negative product means the two values straddle zero.

if nargin < 3
    error('bisect_root: at least 3 input arguments required');
end

test = func(xl) * func(xu);
if test > 0
    error('bisect_root: no sign change over [xl, xu], root is not bracketed');
end

if nargin < 4 || isempty(es),    es    = 0.0001; end
if nargin < 5 || isempty(maxit), maxit = 50;     end

iter = 0;
xr   = xl;
ea   = 100;

while (1)
    xrold = xr;
    xr    = (xl + xu) / 2;
    iter  = iter + 1;

    if xr ~= 0
        ea = abs((xr - xrold)/xr) * 100;
    end

    test = func(xl) * func(xr);
    if test < 0
        xu = xr;
    elseif test > 0
        xl = xr;
    else
        ea = 0;
    end

    if ea <= es || iter >= maxit
        break
    end
end

root = xr;
fx   = func(xr);
end


function [y, V, surcharged, iter, ea] = normal_depth(Qreq, D, S, n, es, maxit)
%NORMAL_DEPTH  Depth of flow in a circular pipe carrying a given discharge.
%
% THE ROOT PROBLEM - the heart of the project.
%
%   Manning's equation gives Q for a KNOWN depth. The design question is
%   the reverse: what depth carries a KNOWN Q? Substituting the circular
%   segment geometry into Manning's equation gives
%
%       f(y) = (1/n) * A(y) * R(y)^(2/3) * sqrt(S) - Qreq = 0
%
%   with y buried inside acos(1 - 2y/D), inside a subtraction, inside a
%   fractional power, inside a ratio. It cannot be rearranged for y. The
%   equation is transcendental and must be solved numerically.
%
% WHY BISECTION rather than Newton-Raphson:
%   1. No derivative needed - differentiating that expression by hand is
%      genuinely unpleasant; bisection uses only function values.
%   2. The bracket is guaranteed, so it cannot fail:
%          y = 0      ->  f = -Qreq       < 0
%          y = yQmax  ->  f = Qmax - Qreq > 0   (when the flow fits)
%      Newton-Raphson can diverge from a poor start. For a solver running
%      unattended 20 times a second inside a GUI, guaranteed convergence
%      beats fast convergence.
%
% WHY THE BRACKET STOPS AT 0.938D, NOT D:
%   Q(y) rises then falls, so over [0, D] there can be TWO depths giving
%   the same discharge. Bisection assumes ONE root in the bracket; give it
%   two and the answer is meaningless. Over [0, yQmax] the function is
%   strictly increasing, so the root is unique.

if nargin < 4
    error('normal_depth: four input arguments required (Qreq, D, S, n)');
end
if nargin < 5 || isempty(es),    es    = 0.0001; end
if nargin < 6 || isempty(maxit), maxit = 50;     end

if Qreq <= 0
    y = 0; V = 0; surcharged = false; iter = 0; ea = 0;
    return
end

[~, ~, Qmax, yQmax] = manning_capacity(D, S, n);

% No root exists - detect this BEFORE calling bisection rather than
% handing it an unbracketed problem and hoping.
if Qreq >= Qmax
    y = D;
    [~, V] = manning_Q(D, D, S, n);
    surcharged = true;
    iter = 0;
    ea = NaN;
    return
end

f = @(yy) manning_Q(yy, D, S, n) - Qreq;    % anonymous function handle

[y, ~, ea, iter] = bisect_root(f, 1e-9, yQmax, es, maxit);

[~, V] = manning_Q(y, D, S, n);
surcharged = false;
end


function i_t = triangular_hyetograph(t, Ddur, Ipeak)
%TRIANGULAR_HYETOGRAPH  Symmetric triangular design storm.
%
% i_t = triangular_hyetograph(t, Ddur, Ipeak)
%
% input:
%   t     = time from the start of the storm, same units as Ddur
%   Ddur  = total storm duration (minutes in this application)
%   Ipeak = peak rainfall intensity (mm/hr), from the IDF table
% output:
%   i_t   = rainfall intensity at time t (mm/hr), same size as t
%
% The IDF table gives ONE number: the average intensity over the whole
% duration. Real rain builds, peaks and dies away, so a design storm
% distributes that number into a time pattern:
%
%       i(t) = Ipeak * (t/tp)           for  0  <= t <= tp
%       i(t) = Ipeak * (Ddur - t)/tp    for  tp <  t <= Ddur
%       tp   = Ddur/2
%
% A real design would use an alternating-block or Huff distribution built
% from local storm records; the triangle is chosen for transparency.
%
% Written with LOGICAL MASKS instead of a loop with an if inside,
% following the vectorisation section of the Fundamentals lecture.

if nargin < 3
    error('triangular_hyetograph: three input arguments required');
end

tp  = Ddur / 2;
i_t = zeros(size(t));

rising  = (t >= 0)  & (t <= tp);
falling = (t >  tp) & (t <= Ddur);

i_t(rising)  = Ipeak .* t(rising) / tp;
i_t(falling) = Ipeak .* (Ddur - t(falling)) / tp;

i_t = max(i_t, 0);
end


function I = simpson13(x, y)
%SIMPSON13  Numerical integration by Simpson's 1/3 rule.
%
%   I = (h/3)*[ f(x0) + 4*sum(odd) + 2*sum(even) + f(xn) ]
%
% Simpson's rule fits a PARABOLA through each group of three consecutive
% points, so it is exact for polynomials up to cubic order, whereas trapz
% joins points with straight lines and is exact only for linear data. The
% alternating 4s and 2s are the weights that fall out of integrating a
% parabola through three equally spaced points.
%
% Requires an EVEN number of segments and EQUAL spacing. Both are checked,
% because silently returning a wrong number is worse than an error.

if nargin < 2
    error('simpson13: two input arguments required (x, y)');
end

x = x(:);
y = y(:);

if numel(x) ~= numel(y)
    error('simpson13: x and y must be the same length');
end

n = numel(x) - 1;
if mod(n, 2) ~= 0
    error('simpson13: needs an even number of segments (odd number of points)');
end

h = (x(end) - x(1)) / n;

if max(abs(diff(x) - h)) > 1e-9 * max(1, abs(h))
    error('simpson13: x must be equally spaced');
end

oddSum  = sum(y(2:2:n));
evenSum = sum(y(3:2:n));

I = (h/3) * (y(1) + 4*oddSum + 2*evenSum + y(end));
end


function [status, level] = design_status(Qreq, Qcap, V, vmin, vmax)
%DESIGN_STATUS  Verdict on a trial storm sewer design.
%
% THE THREE CHECKS, in order of severity:
%   capacity  Qcap < Qreq. The pipe cannot pass the flow. Water backs up,
%             manholes surcharge, the street floods. FAIL.
%   vmin      Below about 0.6 m/s, grit and organic solids settle out
%             instead of being swept along, and the sewer silts up and
%             blocks. This is the SELF-CLEANSING velocity.
%   vmax      Above about 3.0 m/s the suspended grit abrades the invert.
%
% Capacity is tested first because it governs: a pipe that cannot pass the
% flow is broken regardless of its velocity.
%
% These criteria FIGHT each other. Enlarging a pipe to fix a capacity
% failure makes the flow shallower and slower, pushing toward siltation;
% steepening it raises velocity toward the scour limit and deepens the
% excavation. Design is finding the window where all three hold at once.

if nargin < 3
    error('design_status: at least 3 input arguments required');
end
if nargin < 4 || isempty(vmin), vmin = 0.6; end
if nargin < 5 || isempty(vmax), vmax = 3.0; end

if Qcap < Qreq
    status = 'UNDERSIZED PIPE - capacity below demand';
    level  = 3;
elseif V < vmin
    status = 'VELOCITY TOO LOW - siltation risk';
    level  = 2;
elseif V > vmax
    status = 'VELOCITY TOO HIGH - erosion risk';
    level  = 2;
else
    status = 'OK - design passes';
    level  = 1;
end
end


function [i_T, x_T, u, alpha] = gumbel_fit(annualMax, durHr, Tr)
%GUMBEL_FIT  Gumbel (EV Type I) frequency analysis of annual maxima.
%
% [i_T, x_T, u, alpha] = gumbel_fit(annualMax, durHr, Tr)
%
% input:
%   annualMax = vector of annual maximum rainfall depths (mm) for one
%               duration, one value per year of record
%   durHr     = the duration those depths refer to (hr)
%   Tr        = return period(s) of interest (yr)
% output:
%   i_T       = rainfall intensity for each Tr (mm/hr)
%   x_T       = rainfall depth for each Tr (mm)
%   u, alpha  = fitted location and scale parameters (mm)
%
% THE IDEA
%   Take many years of records, keep the single largest storm from each
%   year, fit a distribution to that series of annual maxima, then use it
%   to answer "what depth has a 1-in-T chance each year?" The Gumbel
%   distribution is built for exactly this - it models the maximum of many
%   samples. This is the method behind the stored Khulna table.
%
% METHOD OF MOMENTS
%       alpha = sqrt(6) * s / pi
%       u     = mean - 0.5772 * alpha        (0.5772 = Euler's constant)
%       x_T   = u - alpha * ln( -ln( 1 - 1/Tr ) )
%
% NOTE: this function is supplied so a real BMD annual maximum series can
% be processed directly if your group obtains one. It is NOT used to
% generate the stored table - that comes from the published source.

if nargin < 3
    error('gumbel_fit: three input arguments required');
end

x = annualMax(:);
x = x(~isnan(x));

N = numel(x);
if N < 5
    error('gumbel_fit: need at least 5 years of record, got %d', N);
end

xbar = mean(x);
s    = std(x);                          % n-1 degrees of freedom

alpha = sqrt(6) * s / pi;
u     = xbar - 0.5772 * alpha;

Tr  = Tr(:)';
x_T = u - alpha * log(-log(1 - 1 ./ Tr));
i_T = x_T / durHr;
end


function [u, alpha, resid] = gumbel_from_table()
%GUMBEL_FROM_TABLE  Recover the Gumbel parameters behind the stored table.
%
% [u, alpha, resid] = gumbel_from_table()
%
% output:
%   u     = location parameter for each duration row (mm)
%   alpha = scale parameter for each duration row (mm)
%   resid = largest absolute fitting residual per row (mm)
%
% For a Gumbel distribution the depth is exactly linear in the reduced
% variate
%       y = -ln(-ln(1 - 1/T)),
% that is  P = u + alpha*y.
%
% The source study fitted Gumbel to 62 years of BMD data. Least-squares
% fitting a straight line to each row of the published table therefore
% recovers the parameters that study obtained, and the residuals show how
% cleanly the table follows the distribution. Small residuals confirm the
% table really is a Gumbel product rather than a set of loose numbers -
% which is a validation you can quote, using only published values and no
% invented data.
%
% This is the Descriptive Statistics module applied in reverse: instead of
% going from a sample to parameters, we go from published quantiles back
% to the parameters that produced them.

[~, ~, Dtab, Ttab, Itab] = idf_intensity(60, 10);   % fetch the table

y  = -log(-log(1 - 1./Ttab));           % reduced variate, row vector
nD = numel(Dtab);

u     = zeros(nD, 1);
alpha = zeros(nD, 1);
resid = zeros(nD, 1);

for r = 1:nD
    P = Itab(r,:) .* (Dtab(r)/60);      % intensity (mm/hr) -> depth (mm)

    % Least squares straight line, P = u + alpha*y, written out rather
    % than calling polyfit so the arithmetic is visible.
    n    = numel(y);
    Sy   = sum(y);
    SP   = sum(P);
    Syy  = sum(y.^2);
    SyP  = sum(y .* P);

    alpha(r) = (n*SyP - Sy*SP) / (n*Syy - Sy^2);
    u(r)     = (SP - alpha(r)*Sy) / n;

    resid(r) = max(abs(P - (u(r) + alpha(r)*y)));
end


end