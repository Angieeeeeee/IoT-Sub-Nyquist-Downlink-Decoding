%% ========================================================================
% Packet Detection + Embedded Preamble Check (10xA) + Symbol Correlation
%
%   - Using Schmidl-Cox Packet Detection algorithm identify all packets
%     (L-STF check)
%   - For each packet, find the exact sampling offset using L-LTF
%     correlation
%   - Slice first 11 Symbols from Payload (Skip first since its the MAC
%     Header) - 10 Symbol Embedded Preamble (AAA...)
%   - First Packet to pass Correlation threshold (0.6 for now) of embedded preamble is
%     chosen
%   - Run symbol to symbol correlation
%   - Estimated CFO (not sure if correct)
%
% Assumptions:
%   - fs = 20 MHz
%   - Non-HT CBW20
%   - Embedded Preamble: 10 identical "A" OFDM data symbols starting at DATA #2
%% ========================================================================

clear; clc;

%% ---------------- USER SETTINGS ----------------
filename      = "/Users/angelina/Desktop/Research/DualUSRPExp4_15_1.dat";  % Put your file Path here to run
fs            = 20e6;                 % capture sample rate (Hz)
fileFormat    = "fc32iq";             % IQ sample dat file
chunkSamples  = 2e6;                  % stream chunk size (complex samples)
overlapSamps  = 6000;                 % overlap between chunks

% Schmidl–Cox parameters (CBW20)
D = 16;                               % repetition spacing in samples
L = 64;                               % sliding sum length

% Refinement window around candidate
refinePre_us  = 200;                  % read this much before candidate idx
refineSpan_us = 2000;                 % read this much total window length

% Embedded preamble check
thrA   = 0.60;                        % threshold on min corr across 10 A symbols
M_A    = 10;                           % number of A symbols in marker
startSymA = 2;                        % marker begins at DATA symbol #2
KneedA = startSymA + M_A - 1;         % need DATA symbols #1..#11

% How many L-STF candidates to test
NkeepGlobal = 50;                     % store top-N peaks across the whole file
NlocPerChunk = 5;                     % pick top-N peaks per chunk

% Packet plot / slicing
Kplot_data_syms  = 50;                % how many DATA symbols to include in plot
Kslice_data_syms = 50;                % how many DATA symbols to slice for analysis

% Known A/B pattern for symbol-to-symbol correlation (comparisons start at DATA #2)
pattern = 'AAAAAAAAAAABBAAABBAAABBAABABABABAABABBAAABBAAABBAA';

%% ---------------- OFDM CONSTANTS (Non-HT, CBW20) ----------------
Nfft = 64;
Ncp  = 16;
Ns   = Nfft + Ncp;                    % 80 samples per OFDM symbol

LSTF_len = 160;
LLTF_len = 160;
LSIG_len = 80;

assert(abs(fs-20e6) < 1, "This script assumes fs=20 MHz.");

%% ---------------- WLAN REFERENCE (L-LTF) ----------------
cfg = wlanNonHTConfig('ChannelBandwidth','CBW20','MCS',0,'PSDULength',200);
lltf_ref = wlanLLTF(cfg);
lltf_ref = lltf_ref(:);

%% ---------------- FILE FORMAT ----------------
switch fileFormat
    case "fc32iq"
        bytesPerComplex = 8;
        readChunk = @(fid, N) localReadFC32IQ(fid, N);
    case "sc16iq"
        bytesPerComplex = 4;
        readChunk = @(fid, N) localReadSC16IQ(fid, N);
    otherwise
        error("fileFormat must be 'fc32iq' or 'sc16iq'.");
end

normcorr = @(a,b) abs((a(:)' * b(:)) / sqrt((a(:)'*a(:)) * (b(:)'*b(:)) + 1e-12));

%% ========================================================================
% PASS 1: STREAM THROUGH FILE AND COLLECT TOP-N SCHMIDL–COX CANDIDATES
% (L-STF correlation)
%% ========================================================================
fid = fopen(filename,"rb");
assert(fid>0, "Couldn't open file.");

candM   = -inf(1, NkeepGlobal);
candIdx = -ones(1, NkeepGlobal);

globalIndex = 0;
tail = complex(single([]), single([]));

while true
    x = readChunk(fid, chunkSamples);
    if isempty(x), break; end

    x = [tail; x];
    N = length(x);

    if N < (L + D + 2)
        tail = x;
        continue;
    end

    a = x(1:end-D);
    b = x(1+D:end);

    prod_ab = a .* conj(b);
    pow_b   = abs(b).^2;

    P = filter(ones(L,1,'single'), 1, prod_ab);
    R = filter(ones(L,1,'single'), 1, pow_b);

    P = P(L:end);
    R = R(L:end);

    M = abs(P) ./ (R + 1e-12);

    % take top peaks from this chunk
    [vals, locs] = maxk(M, min(NlocPerChunk, numel(M)));

    for ii = 1:numel(locs)
        Mpk = vals(ii);
        n_in_x = locs(ii) + (L-1);
        abs_idx = globalIndex + (n_in_x - 1);

        [worstVal, worstPos] = min(candM);
        if Mpk > worstVal
            candM(worstPos)   = Mpk;
            candIdx(worstPos) = abs_idx;
        end
    end

    % advance with overlap
    keep = min(overlapSamps, N);
    tail = x(end-keep+1:end);
    globalIndex = globalIndex + (N - keep);
end

fclose(fid);

% sort candidates best-first + remove invalid
[~, order] = sort(candM, 'descend');
candM   = candM(order);
candIdx = candIdx(order);

valid = isfinite(candM) & candIdx >= 0;
candM   = candM(valid);
candIdx = candIdx(valid);

assert(~isempty(candIdx), "No S&C candidates found (unexpected).");


%% ========================================================================
% PASS 2: FOR EACH CANDIDATE, REFINE WITH L-LTF + RUN EMBEDDED AAAAA CHECK
%% ========================================================================
refPre  = round(refinePre_us  * 1e-6 * fs);
refSpan = round(refineSpan_us * 1e-6 * fs);

found = false;

% variables we will set when the correct packet is found
idx_lstf = -1; idx_lltf = -1; idx_lsig = -1; idx_data = -1;
peakVal_best = NaN;
rhoA_best = [];

for ci = 1:numel(candIdx)

    % --- Read refine window around candidate ---
    refStart = max(0, candIdx(ci) - refPre);
    refLen   = refSpan;

    fid = fopen(filename,"rb"); assert(fid>0);
    fseek(fid, refStart * bytesPerComplex, "bof");
    x_ref = readChunk(fid, refLen);
    fclose(fid);

    if length(x_ref) < LLTF_len + 20
        continue;
    end

    % --- L-LTF correlation refine ---
    c = conv(x_ref, flipud(conj(lltf_ref)), "valid");
    [peakVal, kmax] = max(abs(c));

    idx_lltf_try = refStart + (kmax - 1);
    idx_lstf_try = idx_lltf_try - LSTF_len;
    idx_lsig_try = idx_lltf_try + LLTF_len;
    idx_data_try = idx_lsig_try + LSIG_len;

    if idx_lstf_try < 0
        continue;
    end

    % --- Read enough DATA symbols for embedded check ---
    dataNeededA = KneedA * Ns;

    fid = fopen(filename,"rb"); assert(fid>0);
    fseek(fid, idx_data_try * bytesPerComplex, "bof");
    x_dataA = readChunk(fid, dataNeededA);
    fclose(fid);

    if length(x_dataA) < dataNeededA
        continue;
    end

    symA_80 = reshape(x_dataA, Ns, KneedA);

    % --- Embedded AAAAA score ---
    [scoreA, rhoA] = embeddedA10_waveform(symA_80, startSymA, M_A);

    fprintf("Cand %2d: S&C M=%.3f | LLTFpeak=%.3g | AAAAA cor=%.3f | each sym cor =[%s]\n", ...
        ci, candM(ci), peakVal, scoreA, sprintf("%.3f ", rhoA));

    if scoreA >= thrA
        fprintf("\n>>> MATCH: Found packet\n");
        fprintf("    idx_lstf=%d | idx_lltf=%d | idx_lsig=%d | idx_data=%d\n\n", ...
            idx_lstf_try, idx_lltf_try, idx_lsig_try, idx_data_try);

        found = true;

        idx_lstf = idx_lstf_try;
        idx_lltf = idx_lltf_try;
        idx_lsig = idx_lsig_try;
        idx_data = idx_data_try;

        peakVal_best = peakVal;
        rhoA_best = rhoA;

        break;
    end
end

assert(found, "No packet matched the embedded preamble (AAAAA) with threshold %.2f.", thrA);

%% ========================================================================
% PASS 3: READ ONLY THE SINGLE PACKET (PREAMBLE+HEADER+K DATA SYMBOLS)
%% ========================================================================
pkt_start = idx_lstf;
pkt_end   = idx_data + Kplot_data_syms*Ns - 1;
pkt_len   = pkt_end - pkt_start + 1;

fid = fopen(filename,"rb"); assert(fid>0);
fseek(fid, pkt_start * bytesPerComplex, "bof");
xpkt = readChunk(fid, pkt_len);
fclose(fid);

tp = (0:length(xpkt)-1).' / fs;

t_lstf = 0;
t_lltf = (idx_lltf - pkt_start)/fs;
t_lsig = (idx_lsig - pkt_start)/fs;
t_data = (idx_data - pkt_start)/fs;

%% ========================================================================
% PLOTS: WAVEFORM (I/Q) + MAGNITUDE, MARKERS
%% ========================================================================
figure;
subplot(2,1,1);
plot(tp*1e6, real(xpkt)); grid on; hold on;
plot(tp*1e6, imag(xpkt)); grid on;
ylabel("I / Q");
title("Selected Packet (I/Q)");
xline(t_lstf*1e6,'--','L-STF','LineWidth',1.5);
xline(t_lltf*1e6,'--','L-LTF','LineWidth',1.5);
xline(t_lsig*1e6,'--','L-SIG','LineWidth',1.5);
xline(t_data*1e6,'--','DATA','LineWidth',1.5);
legend("I","Q","Location","best");

subplot(2,1,2);
plot(tp*1e6, abs(xpkt)); grid on;
xlabel("Time (\mus)");
ylabel("|x[n]|");
title("Selected Packet Magnitude");
xline(t_lstf*1e6,'--','L-STF','LineWidth',1.5);
xline(t_lltf*1e6,'--','L-LTF','LineWidth',1.5);
xline(t_lsig*1e6,'--','L-SIG','LineWidth',1.5);
xline(t_data*1e6,'--','DATA','LineWidth',1.5);


%% ========================================================================
% PASS 4: SLICE K DATA SYMBOLS (CP included + CP removed) FOR ANALYSIS
%% ========================================================================
K = Kslice_data_syms;
dataNeeded = K * Ns;

fid = fopen(filename,"rb"); assert(fid>0);
fseek(fid, idx_data * bytesPerComplex, "bof");
x_data = readChunk(fid, dataNeeded);
fclose(fid);

assert(length(x_data) == dataNeeded, "Not enough samples after idx_data in file.");

sym_td   = reshape(x_data, Ns, K);      % 80 x K
sym_noCP = sym_td(Ncp+1:end, :);        % 64 x K

%% ========================================================================
% SNR ESTIMATE (prints only): Noise from pre-packet with guard padding
%% ========================================================================

noiseWin_us = 25;                          % noise window length (us)
noiseWin    = round(noiseWin_us*1e-6*fs);  % samples
padSyms     = 0;                           % symbols to ignore before packet
padSamps    = padSyms * Ns;                % samples to skip
epsP        = 1e-12;

% Define noise window END safely before packet
noiseEnd   = max(0, idx_lstf - padSamps);
noiseStart = max(0, noiseEnd - noiseWin);
noiseLen   = noiseEnd - noiseStart;

assert(noiseLen > 0, "Noise window collapsed. Reduce padding or increase noiseWin_us.");

fid = fopen(filename,"rb"); assert(fid>0);
fseek(fid, noiseStart * bytesPerComplex, "bof");
x_noise = readChunk(fid, noiseLen);
fclose(fid);

assert(~isempty(x_noise), "Noise window empty. Check padding/window sizes.");

% Noise power
Pn = mean(abs(x_noise).^2);

% Signal+noise power from DATA region (already read)
Prx = mean(abs(x_data).^2);

% SNR
Ps      = max(Prx - Pn, epsP);
snr_lin = Ps / max(Pn, epsP);
snr_db  = 10*log10(snr_lin);

fprintf("\n=== SNR estimate (power-based, padded) ===\n");
fprintf("Ignored %d symbols (%d samples) before idx_lstf\n", padSyms, padSamps);
fprintf("Noise window: %d samples (%.1f us)\n", noiseLen, noiseLen/fs*1e6);
fprintf("Pn  (noise power)      = %.3e\n", Pn);
fprintf("Prx (signal+noise pow) = %.3e\n", Prx);
fprintf("SNR = %.2f dB\n\n", snr_db);

fprintf("mean|x_noise| = %.3e,  max|x_noise| = %.3e\n", ...
        mean(abs(x_noise)), max(abs(x_noise)));
fprintf("mean|x_data|  = %.3e,  max|x_data|  = %.3e\n", ...
        mean(abs(x_data)),  max(abs(x_data)));

figure; plot(abs(x_noise)); grid on;
title("Noise window magnitude (with padding)");

figure; plot(abs(x_data)); grid on;
title("DATA window magnitude");

%% ========================================================================
% PASS 5: REFERENCE-SYMBOL CORRELATION + A/B DECISION (CP INCLUDED 80)
% - Save DATA symbol #1 as reference
% - Correlate every DATA symbol k with reference
% - If corr >= thrAB => label 'A' else 'B'
% - Plot correlation values vs symbol index
%% ========================================================================

% --- user knob for A/B decision ---
thrAB = 0.6;   

% Use DATA #1 as reference
refSym = sym_td(:, 2);

rho_ref = zeros(K,1);
labelsAB = repmat('B', K, 1);  % default B

for k = 1:K
    rho_ref(k) = normcorr(refSym, sym_td(:, k));
    if rho_ref(k) >= thrAB
        labelsAB(k) = 'A';
    end
end

% Print a compact result line
fprintf("=== Reference-based A/B classification (ref = DATA #1, thr=%.2f) ===\n", thrAB);
fprintf("Symbols 1..%d labels:\n%s\n\n", K, string(labelsAB.'));

% Plot correlation vs symbol index
figure;
stem(1:K, rho_ref, 'filled', 'LineWidth', 1.5); grid on;
ylim([0 1.05]);
xlabel("DATA symbol index k");
ylabel("|ρ(ref, sym k)|");
title("Correlation of Each DATA Symbol vs Reference (DATA #1) [CP included]");
xline(1,'--','Reference (DATA #1)','LineWidth',1.5);
yline(thrAB,'--',sprintf("thrAB=%.2f",thrAB),'LineWidth',1.5);

% Optional: show A/B labels on top of points
hold on;
for k = 1:K
    text(k, min(1.03, rho_ref(k)+0.03), labelsAB(k), 'HorizontalAlignment','center');
end
hold off;

%% ========================================================================
% SIMPLE CFO ESTIMATE USING THE TWO L-LTF SYMBOLS
%% ========================================================================
lltf_pkt_start = (idx_lltf - pkt_start) + 1;

s1 = xpkt(lltf_pkt_start + Ncp : lltf_pkt_start + Ncp + (Nfft-1));
s2 = xpkt(lltf_pkt_start + Ns  + Ncp : lltf_pkt_start + Ns  + Ncp + (Nfft-1));

phi = angle(s1' * conj(s2));
cfo_cycles_per_sample = phi / (2*pi*Ns);
cfo_hz = cfo_cycles_per_sample * fs;

fprintf("\nL-LTF CFO estimate:\n");
fprintf("  phi = %.6f rad\n", phi);
fprintf("  CFO = %.6e cycles/sample\n", cfo_cycles_per_sample);
fprintf("  CFO = %.3f Hz\n", cfo_hz);

%% ========================================================================
% LOCAL HELPERS
%% ========================================================================
function [score, rho_vec] = embeddedA10_waveform(sym80, startSym, M_A)
% Waveform-only, CP-included (80) embedded Ax check.
% Correlates DATA startSym vs startSym+1 ... startSym+M_A-1
% Uses the same complex normalized correlation as symbol to symbol.

    normcorr = @(a,b) abs((a(:)' * b(:)) / sqrt((a(:)'*a(:)) * (b(:)'*b(:)) + 1e-12));

    lastNeeded = startSym + M_A - 1;
    if size(sym80,2) < lastNeeded
        rho_vec = [];
        score = -Inf;
        return;
    end

    ref = sym80(:, startSym);

    rho_vec = zeros(M_A-1,1);
    for m = 1:(M_A-1)
        ksym = startSym + m;
        rho_vec(m) = normcorr(ref, sym80(:, ksym));
    end

    score = min(rho_vec);   % worst of the 9 comparisons (for M_A=10)
end



function x = localReadFC32IQ(fid, N)
    raw = fread(fid, 2*N, "float32=>single");
    if isempty(raw)
        x = complex(single([]), single([]));
        return;
    end
    raw = reshape(raw, 2, []);
    x = complex(raw(1,:), raw(2,:)).';
end

function x = localReadSC16IQ(fid, N)
    raw = fread(fid, 2*N, "int16=>single");
    if isempty(raw)
        x = complex(single([]), single([]));
        return;
    end
    raw = reshape(raw, 2, []);
    x = (complex(raw(1,:), raw(2,:)).') / 32768;
end

%{
%% ========================================================================
% PASS 5: SYMBOL-TO-SYMBOL CORRELATION (PATTERN) — CP INCLUDED (80)
% comparisons: DATA #2->#3 ... based on 'pattern'
%% ========================================================================
Np  = length(pattern);
Ncmp = Np - 1;

numSymsNeeded = Np + 1;      % need DATA #1..#(Np+1)
dataNeeded2   = numSymsNeeded * Ns;

fid = fopen(filename,"rb"); assert(fid>0);
fseek(fid, idx_data * bytesPerComplex, "bof");
x_data2 = readChunk(fid, dataNeeded2);
fclose(fid);

if length(x_data2) < dataNeeded2
    warning('File ended early: got %d samples, need %d. Shortening pattern.', length(x_data2), dataNeeded2);
    numSymsAvailable = floor(length(x_data2) / Ns);
    if numSymsAvailable < 3
        error('Not enough DATA symbols in capture. Need at least 3.');
    end
    pattern = pattern(1:min(length(pattern), numSymsAvailable-1));
    Np  = length(pattern);
    Ncmp = Np - 1;
    numSymsNeeded = Np + 1;
    x_data2 = x_data2(1:numSymsNeeded*Ns);
end

sym2 = reshape(x_data2, Ns, numSymsNeeded);

rho    = zeros(Ncmp,1);
labels = strings(Ncmp,1);

for k = 1:Ncmp
    s1 = sym2(:, k+1);               % DATA #2, #3, ...
    s2 = sym2(:, k+2);               % DATA #3, #4, ...
    rho(k) = normcorr(s1, s2);
    labels(k) = pattern(k) + "→" + pattern(k+1) + " (sym " + (k+1) + "→" + (k+2) + ")";
end

figure;
stem(rho, 'filled', 'LineWidth', 1.5); grid on;
ylim([0 1.05]);
xlabel("Comparison index (starting at DATA symbol 2)");
ylabel("|ρ| (normalized complex correlation, CP included)");
title("Symbol-to-Symbol Correlation (80-sample symbols, CP included)");
xticks(1:Ncmp);
xticklabels(labels);
xtickangle(45);
%}

%{
%% ========================================================================
% PASS 5: REFERENCE-SYMBOL CORRELATION + A/B DECISION (CP INCLUDED 80)
% - Save DATA symbol #1 as reference
% - Correlate every DATA symbol k with reference
% - If corr >= thrAB => label 'A' else 'B'
% - Plot correlation values vs symbol index
%% ========================================================================

% --- user knob for A/B decision ---
thrAB = 0.1;   

% Use DATA #1 as reference
refSym = sym_td(:, 1);

rho_ref = zeros(K,1);
labelsAB = repmat('B', K, 1);  % default B

for k = 1:K
    rho_ref(k) = normcorr(refSym, sym_td(:, k));
    if rho_ref(k) >= thrAB
        labelsAB(k) = 'A';
    end
end

% Print a compact result line
fprintf("=== Reference-based A/B classification (ref = DATA #1, thr=%.2f) ===\n", thrAB);
fprintf("Symbols 1..%d labels:\n%s\n\n", K, string(labelsAB.'));

% Plot correlation vs symbol index
figure;
stem(1:K, rho_ref, 'filled', 'LineWidth', 1.5); grid on;
ylim([0 1.05]);
xlabel("DATA symbol index k");
ylabel("|ρ(ref, sym k)|");
title("Correlation of Each DATA Symbol vs Reference (DATA #1) [CP included]");
xline(1,'--','Reference (DATA #1)','LineWidth',1.5);
yline(thrAB,'--',sprintf("thrAB=%.2f",thrAB),'LineWidth',1.5);

% Optional: show A/B labels on top of points
hold on;
for k = 1:K
    text(k, min(1.03, rho_ref(k)+0.03), labelsAB(k), 'HorizontalAlignment','center');
end
hold off;
%}





%{
%% ========================================================================
% PASS 5: Dual-reference A/B classification using DATA #1 and DATA #6 as 
%         A templates
%% ========================================================================

thrAB = 0.20;     % tune this
K = size(sym_td, 2);

ref1_idx = 2;
ref6_idx = 6;

assert(K >= ref6_idx, "Need at least %d symbols, but K=%d.", ref6_idx, K);

ref1 = sym_td(:, ref1_idx);
ref6 = sym_td(:, ref6_idx);

rho1 = zeros(K,1);
rho6 = zeros(K,1);
labelsAB = repmat('B', K, 1);   % default B

for k = 1:K
    rho1(k) = normcorr(ref1, sym_td(:,k));
    rho6(k) = normcorr(ref6, sym_td(:,k));

    if (rho1(k) >= thrAB) || (rho6(k) >= thrAB)
        labelsAB(k) = 'A';
    end
end

detectedPattern = string(labelsAB.');

fprintf("=== Dual-ref A/B classification (ref=DATA#%d and DATA#%d, thr=%.2f) ===\n", ...
        ref1_idx, ref6_idx, thrAB);
fprintf("Detected pattern (1..%d):\n%s\n\n", K, detectedPattern);

% ---------- "Correlation tables" printout ----------
T = table((1:K).', rho1, rho6, labelsAB, ...
          'VariableNames', {'SymIdx','CorrToRef1','CorrToRef6','Label'});
disp(T);

% ---------- Plot correlations ----------
figure;
stem(1:K, rho1, 'filled', 'LineWidth', 1.2); grid on;
ylim([0 1.05]);
xlabel("DATA symbol index k");
ylabel("|ρ(DATA#1, sym k)|");
title("Correlation vs Reference DATA #1 [CP included]");
xline(ref1_idx,'--','Ref #1','LineWidth',1.5);
yline(thrAB,'--',sprintf("thrAB=%.2f",thrAB),'LineWidth',1.5);

figure;
stem(1:K, rho6, 'filled', 'LineWidth', 1.2); grid on;
ylim([0 1.05]);
xlabel("DATA symbol index k");
ylabel("|ρ(DATA#6, sym k)|");
title("Correlation vs Reference DATA #6 [CP included]");
xline(ref6_idx,'--','Ref #6','LineWidth',1.5);
yline(thrAB,'--',sprintf("thrAB=%.2f",thrAB),'LineWidth',1.5);

% ---------- Optional: combined view (easier to see "A if either") ----------
figure;
plot(1:K, rho1, '-o', 1:K, rho6, '-o', 'LineWidth', 1.2); grid on;
ylim([0 1.05]);
xlabel("DATA symbol index k");
ylabel("Correlation magnitude");
title("Correlation to DATA #1 and DATA #6 (A if either ≥ thr)");
yline(thrAB,'--',sprintf("thrAB=%.2f",thrAB),'LineWidth',1.5);
legend("|ρ(ref1)|","|ρ(ref6)|","Location","best");

% Annotate labels on combined plot
hold on;
for k = 1:K
    y = max(rho1(k), rho6(k));
    text(k, min(1.03, y+0.03), labelsAB(k), 'HorizontalAlignment','center');
end
hold off;
%}