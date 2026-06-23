clc; clear; close all;
try
    release(tx);
    release(rx);
catch
end

clear tx rx

load("x_corrupt.mat");

x = x_corrupt(:);
fs_in = 48000;

%% Low-pass filter and downsample to 8 kHz

fs = 8000;

lpFilt = designfilt('lowpassfir', ...
    'PassbandFrequency', 3500, ...
    'StopbandFrequency', 4000, ...
    'PassbandRipple', 1, ...
    'StopbandAttenuation', 60, ...
    'SampleRate', fs_in);

x_filt = filtfilt(lpFilt, x);       % Anti-aliasing filter
x8 = resample(x_filt, fs, fs_in);   % Downsample to 8 kHz

% sound(x8, fs); % Optional

%% Frame settings

frameLen = 160;      % 20 ms at 8 kHz
overlapLen = 40;     
hopLen = frameLen - overlapLen;
overlapPercentage = overlapLen/frameLen*100;
fprintf("The overlap percentage: %d Percent \n", overlapPercentage);
p = 10;              % LPC order

frames = buffer(x8, frameLen, overlapLen, 'nodelay');
numFrames = size(frames, 2);

win = hamming(frameLen);

disp("Number of frames:");
disp(numFrames);

LSF = zeros(numFrames, p);
Gain = zeros(numFrames, 1);
voicedFlag = zeros(numFrames, 1);
pitch = zeros(numFrames, 1);

%% LPC and LSF extraction + voiced/unvoiced + pitch

for k = 1:numFrames

    frame = frames(:, k) .* win;

    [A, E] = lpc(frame, p);
    LSF(k,:) = poly2lsf(A);
    Gain(k) = sqrt(E);

    frameEnergy = sum(frame.^2) / length(frame);

    R = xcorr(frame);
    R = R(length(frame):end);

    % Valid pitch lag range for 50-400 Hz
    minLag = round(fs/500); % 8k/500=16
    maxLag = round(fs/50); % 8k/50=160, 16<=lag<=160

    % Only search inside valid lag range
    R_pitch = R(minLag:maxLag);
    
    [peakVal, peakLoc] = max(R_pitch);

    % Convert local index back to actual lag
    pitchLag = peakLoc + minLag - 1;
    
    if frameEnergy > 1e-4 && peakVal > 0.3*max(R)
        voicedFlag(k) = 1;
        pitch(k) = fs / pitchLag; % pitch estimation
    else
        voicedFlag(k) = 0;
        pitch(k) = 0;
    end
end



%% Frequency spectrum comparison

N = length(x);
f = linspace(-fs_in/2, fs_in/2, N);

X_original = fftshift(fft(x));
X_filtered = fftshift(fft(x_filt));

figure;

subplot(2,2,1);
plot(x);
xlim([1 length(x)]);
title("Original Corrupted Speech Signal");
xlabel("Sample");
ylabel("Amplitude");

subplot(2,2,2);
plot(x8);
xlim([1 length(x8)]);
title("Prepared Speech Signal at 8 kHz");
xlabel("Sample");
ylabel("Amplitude");

subplot(2,2,3);
plot(f, abs(X_original));
title("Spectrum of Original Corrupted Speech");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

subplot(2,2,4);
plot(f, abs(X_filtered));
title("Spectrum After Low-pass Filtering");
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

%% Data size comparison

raw_bits_per_sec = fs * 16;

lsfBits = 10;
gainBits = 8;
pitchBits = 8;
voicedBits = 1;

bits_per_frame = p*lsfBits + gainBits + pitchBits + voicedBits;
coded_bits_per_sec = bits_per_frame * (fs/frameLen);

fprintf("\n--- Data Size Comparison ---\n");
fprintf("Raw speech bitrate: %.2f kbps\n", raw_bits_per_sec/1000);
fprintf("LPC/LSF parameter bitrate: %.2f kbps\n", coded_bits_per_sec/1000);
fprintf("Compression ratio: %.2f\n", raw_bits_per_sec/coded_bits_per_sec);

%% =========================================================
%% QUANTIZATION + PACKETIZATION + BPSK + USRP TX/RX
%% =========================================================

%% Packet settings

preamble = repmat([1;0;1;0;1;1;0;0], 4, 1); % 32-bit preamble
frameIDBits = 6;                             % 5 is enough until 32 frames, 6 is for 64 frames

pitchMin = 50;
pitchMax = 400;

gainMin = min(Gain);
gainMax = max(Gain);

%% Build TOTAL BITSTREAM

totalBitstream = [];
bits = [];

for k = 1:numFrames

    %% Quantize LSF

    lsf_q = round( ...
        (LSF(k,:) - 0) / pi * (2^lsfBits - 1) ...
    );

    lsf_q = max(0, min(2^lsfBits - 1, lsf_q));

    %% Quantize Gain

    if gainMax == gainMin
        gain_q = 0;
    else
        gain_q = round( ...
            (Gain(k) - gainMin) / (gainMax - gainMin) ...
            * (2^gainBits - 1) ...
        );
    end

    gain_q = max(0, min(2^gainBits - 1, gain_q));

    %% Quantize Pitch

    if voicedFlag(k) == 1

        pitch_q = round( ...
            (pitch(k) - pitchMin) / (pitchMax - pitchMin) ...
            * (2^pitchBits - 1) ...
        );

        pitch_q = max(0, min(2^pitchBits - 1, pitch_q));

    else
        pitch_q = 0;
    end

    %% Convert to bits

    lsf_bits = de2bi(lsf_q, lsfBits, 'left-msb');
    gain_bits = de2bi(gain_q, gainBits, 'left-msb');
    pitch_bits = de2bi(pitch_q, pitchBits, 'left-msb');

    voiced_bits = voicedFlag(k);

    %% Create 91-bit payload

    frameBits = [
        reshape(lsf_bits.', [], 1);
        gain_bits(:);
        pitch_bits(:);
        voiced_bits
    ];

    bits = [bits; frameBits];

    %% Frame ID bits

    frameID = de2bi(k-1, frameIDBits, 'left-msb').';

    %% Create packet

    packetBits = [
        preamble;
        frameID;
        frameBits
    ];

    %% Repeat same packet 5 times

    repeatedPacket = repmat(packetBits, 3, 1);

    %% Add to total stream

    totalBitstream = [totalBitstream; repeatedPacket];

end

fprintf("\n--- TOTAL BITSTREAM ---\n");
fprintf("Total transmitted bits = %d\n", length(totalBitstream));

%% =========================================================
%% BPSK MODULATION + RRC
%% =========================================================

sps = 8;
rolloff = 0.35;
span = 10;

rootraisedcos = rcosdesign(rolloff, span, sps, "sqrt");

symbols = 2*totalBitstream - 1;


%% TX Constellation Diagram
figure;
plot(real(symbols), imag(symbols), 'y.', 'MarkerSize', 8);
title('TX Constellation Diagram (BPSK - After Modulation)');
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
xlim([-2 2]); ylim([-2 2]); grid on; axis square;
xline(0,'k--','LineWidth',1); yline(0,'k--','LineWidth',1);
figure;

subplot(2,1,1)
stem(real(symbols(1:100)), 'filled')
title("TX BPSK Symbols")
xlabel("Symbol Index")
ylabel("Amplitude")
grid on

symbols = complex(symbols);

txSig = upfirdn(symbols, rootraisedcos, sps, 1);

txSig = txSig / max(abs(txSig));

subplot(2,1,2)
plot(real(txSig(1:1000)))
title("TX Baseband Signal After RRC ")
xlabel("Sample Index"), ylabel("Amplitude"), grid on
fprintf("Transmit samples = %d\n", length(txSig));


%% =========================================================
%% USRP SETTINGS
%% =========================================================

radioID = "31AE1A0";

fc = 0.433e9;

platformType = "B200";

masterClockRate = 20e6;
interpDecim = 500;

Fs_usrp = masterClockRate / interpDecim;
SamplesPerFrame = 10000;

%% TX object

tx = comm.SDRuTransmitter( ...
    Platform = platformType, ...
    SerialNum = radioID, ...
    CenterFrequency = fc, ...
    Gain = 15, ...
    MasterClockRate = masterClockRate, ...
    InterpolationFactor = interpDecim);

%% RX object

rx = comm.SDRuReceiver( ...
    Platform = platformType, ...
    SerialNum = radioID, ...
    CenterFrequency = fc, ...
    Gain = 30, ...
    MasterClockRate = masterClockRate, ...
    DecimationFactor = interpDecim, ...
    SamplesPerFrame = SamplesPerFrame, ...
    OutputDataType = "double");

%% =========================================================
%% RX CONTINUOUS LISTEN
%% =========================================================

rxSig = [];

fprintf("\n--- RX STARTED ---\n");

% Start RX buffer
for k = 1:10
    [data, len] = rx();

    if len > 0
        rxSig = [rxSig; data(:)];
    end
end

%% =========================================================
%% TRANSMIT BIG SIGNAL
%% =========================================================

fprintf("\n--- TX STARTED ---\n");

rxDuration = 8; % seconds
numRxLoops = 30;
%2*ceil(rxDuration * Fs_usrp / SamplesPerFrame);

for k = 1:4

    fprintf("Transmission repeat %d/4\n", k);
    
    tx(txSig);

    % Keep listening while TX bursts happen

    for r = 1:numRxLoops

        [data, len] = rx();

        if len > 0
            rxSig = [rxSig; data(:)];
        end
    end
end

fprintf("\nTransmission complete.\n");
fprintf("Received samples = %d\n", length(rxSig));

%% =========================================================
%% CLEANUP
%% =========================================================

release(tx);
release(rx);

%% =========================================================
%% RX DEMODULATION + PACKET DETECTION + BER EVALUATION
%% =========================================================

fprintf("\n--- RX DEMODULATION STARTED ---\n");

rxSig = rxSig(:);

% DC removal + normalization (std ile — ani gürültü tepelerine dayanıklı)
rxSig = rxSig - mean(rxSig);
% rxSig = rxSig / (std(rxSig) + eps);

% Matched filter
rxFiltered = filter(rootraisedcos, 1, rxSig);

totalDelay = span * sps;

bestScore = -inf;
bestOffset = 1;
bestBits = [];

%% CFO (Carrier Frequency Offset) Estimation and Correction
% BPSK sinyalinin karesi alınınca CFO frekansı 2x'e katlanır,
% FFT ile tepe noktası bulunur ve düzeltme yapılır.

%rxSq = rxFiltered .^ 2;
%N_fft = 2^nextpow2(length(rxSq));
%RxSq_fft = fft(rxSq, N_fft);
%[~, cfo_idx] = max(abs(RxSq_fft(1:N_fft/2)));
%cfo_est_hz = (cfo_idx - 1) / N_fft * Fs_usrp;  % tek taraflı, Hz
%cfo_est_hz = cfo_est_hz / 2;                    % BPSK'da /2 gerekli

%fprintf("Estimated CFO = %.2f Hz\n", cfo_est_hz);

%t_rx = (0:length(rxFiltered)-1).' / Fs_usrp;
%rxFiltered = rxFiltered .* exp(-1j * 2 * pi * cfo_est_hz .* t_rx);

%% Try all symbol timing offsets
bestRxSamples = [];

for offset = 1:sps

    startIndex = totalDelay + offset;

    if startIndex > length(rxFiltered)
        continue;
    end

    rxSamples = rxFiltered(startIndex:sps:end);

    

    if isempty(rxSamples)
        continue;
    end

    % max() yerine std() kullan — ani gürültü tepelerine karşı dayanıklı
    rxSamples = rxSamples / (std(rxSamples) + eps);

    % BPSK phase correction
    phaseOffset = angle(mean(rxSamples.^2)) / 2;
    rxSamples = rxSamples * exp(-1j*phaseOffset);

    rxBits_candidate = real(rxSamples) > 0;

    % --- Polarity correction: preamble ile normal vs ters karşılaştır ---
    if length(rxBits_candidate) >= length(preamble)
        score_normal = sum(rxBits_candidate(1:length(preamble)) == preamble(:));
        score_inv    = sum(~rxBits_candidate(1:length(preamble)) == preamble(:));
        if score_inv > score_normal
            rxBits_candidate = ~rxBits_candidate;
        end
    end

    rxBits = rxBits_candidate;

    % Score this offset by total preamble matches
    scoreSum = 0;
    searchLimit = length(rxBits) - length(preamble) + 1;

    for idx = 1:searchLimit

        candidate = rxBits(idx:idx+length(preamble)-1);
        score = sum(candidate(:) == preamble(:));

        if score >= 0.8 * length(preamble)
            scoreSum = scoreSum + score;
        end
    end

    if scoreSum > bestScore
        bestScore = scoreSum;
        bestOffset = offset;
        bestBits = rxBits;
        bestRxSamples = rxSamples;
    end
end

figure;

subplot(1,1,1)
plot(real(bestRxSamples(1e5:1e5+1000)))
title("RX Sampled Symbols")
xlabel("Symbol Index")
ylabel("Amplitude")
grid on



fprintf("Best timing offset = %d\n", bestOffset);
fprintf("Recovered hard bits = %d\n", length(bestBits));
fprintf("Preamble score sum = %.2f\n", bestScore);

%% Packet detection

packetLen = length(preamble) + frameIDBits + bits_per_frame;

bestPayloadPerFrame = cell(numFrames, 1);
bestScorePerFrame = -inf(numFrames, 1);

detectedPackets = 0;

searchLimit = length(bestBits) - packetLen + 1;

for idx = 1:searchLimit

    candidatePreamble = bestBits(idx : idx + length(preamble) - 1);
    score = sum(candidatePreamble(:) == preamble(:));

    if score >= 0.8 * length(preamble)

        frameIDStart = idx + length(preamble);
        frameIDEnd = frameIDStart + frameIDBits - 1;

        payloadStart = frameIDEnd + 1;
        payloadEnd = payloadStart + bits_per_frame - 1;

        if payloadEnd <= length(bestBits)

            frameID_rx = bestBits(frameIDStart:frameIDEnd).';
            frameNumber = bi2de(frameID_rx, 'left-msb') + 1;

            if frameNumber >= 1 && frameNumber <= numFrames

                payload = bestBits(payloadStart:payloadEnd);
                detectedPackets = detectedPackets + 1;

                if score > bestScorePerFrame(frameNumber)
                    bestScorePerFrame(frameNumber) = score;
                    bestPayloadPerFrame{frameNumber} = payload(:);
                end
            end
        end
    end
end

%% Reconstruct received payload stream in correct frame order

rxPayloadBits = [];
missingFrames = [];

for k = 1:numFrames

    if isempty(bestPayloadPerFrame{k})
        missingFrames = [missingFrames; k];
    else
        rxPayloadBits = [rxPayloadBits; bestPayloadPerFrame{k}];
    end
end

fprintf("\n--- PACKET DETECTION RESULTS ---\n");
fprintf("Detected packet candidates = %d\n", detectedPackets);
fprintf("Recovered frames = %d / %d\n", numFrames - length(missingFrames), numFrames);

if ~isempty(missingFrames)
    fprintf("Missing frames: ");
    fprintf("%d ", missingFrames);
    fprintf("\n");
end

fprintf("Recovered payload bits = %d\n", length(rxPayloadBits));
fprintf("Expected payload bits = %d\n", length(bits));

%% BER evaluation

numUseful = min(length(rxPayloadBits), length(bits));

if numUseful == 0
    error("No payload bits recovered. BER cannot be calculated.");
end

BER_payload = sum(rxPayloadBits(1:numUseful) ~= bits(1:numUseful)) / numUseful;

fprintf("\n--- BER RESULT ---\n");
fprintf("Compared bits = %d\n", numUseful);
fprintf("Payload BER = %.6f\n", BER_payload);

%% Optional inverted BER check

BER_inverted = sum(~rxPayloadBits(1:numUseful) ~= bits(1:numUseful)) / numUseful;

fprintf("Inverted BER = %.6f\n", BER_inverted);

if BER_inverted < BER_payload
    fprintf("NOTE: Inverted BER is lower. Bit polarity may be flipped.\n");
end

%% RX Constellation Diagram (paket içi semboller, Gardner timing)
c_allSamples = [];
c_searchLimit = length(bestBits) - packetLen + 1;

for c_idx = 1:c_searchLimit
    c_cand  = bestBits(c_idx : c_idx + length(preamble) - 1);
    c_score = sum(c_cand(:) == preamble(:));

    if c_score >= 0.8 * length(preamble)
        c_sampleEnd = c_idx + packetLen - 1;
        if c_sampleEnd <= length(bestRxSamples)
            c_chunk = bestRxSamples(c_idx : c_sampleEnd);
            c_ph    = angle(mean(c_chunk.^2)) / 2;
            c_chunk = c_chunk .* exp(-1j * c_ph);
            c_allSamples = [c_allSamples; c_chunk];
        end
    end
end

figure;
scatter(real(c_allSamples), imag(c_allSamples), 10, 'filled', 'MarkerFaceColor', '#77AC30');
title(sprintf('RX Constellation Diagram - BPSK (BER = %.4f)', BER_payload));
xlabel('In-Phase (I)'); ylabel('Quadrature (Q)');
xline(0,'--k'); yline(0,'--k');
axis([-2 2 -2 2]); grid on;

clearvars c_allSamples c_searchLimit c_idx c_cand c_score c_sampleEnd c_chunk c_ph


%% Diagnostic plots

figure;

subplot(3,1,1)
plot(real(rxSig))
title("Received Signal - Real Part")
xlabel("Sample Index")
ylabel("Amplitude")
grid on

subplot(3,1,2)
stem(bestBits(1:min(200,length(bestBits))), 'filled')
title("First Recovered Hard Bits")
xlabel("Bit Index")
ylabel("Bit")
ylim([-0.2 1.2])
grid on

subplot(3,1,3)
bar(bestScorePerFrame)
title("Best Preamble Score per Frame")
xlabel("Frame Index")
ylabel("Preamble Score")
grid on

%% =========================================================
%% SPEECH RECONSTRUCTION FROM RECEIVED PAYLOADS
%% =========================================================

fprintf("\n--- SPEECH RECONSTRUCTION STARTED ---\n");

% Decode each recovered frame separately.
% If a frame is missing, reuse previous valid parameters.
LSF_rec = zeros(numFrames, p);
Gain_rec = zeros(numFrames, 1);
pitch_rec = zeros(numFrames, 1);
voiced_rec = zeros(numFrames, 1);

lastLSF = LSF(1,:);
lastGain = mean(Gain);
lastPitch = 0;
lastVoiced = 0;

for k = 1:numFrames

    if isempty(bestPayloadPerFrame{k})
        % Missing frame handling
        LSF_rec(k,:) = lastLSF;
        Gain_rec(k) = 0.2 * lastGain;   % reduce energy for missing frame
        pitch_rec(k) = lastPitch;
        voiced_rec(k) = 0;
        continue;
    end

    frameBits_rx = bestPayloadPerFrame{k}(:);

    if length(frameBits_rx) ~= bits_per_frame
        warning("Frame %d has wrong bit length. Using previous parameters.", k);
        LSF_rec(k,:) = lastLSF;
        Gain_rec(k) = 0.2 * lastGain;
        pitch_rec(k) = lastPitch;
        voiced_rec(k) = 0;
        continue;
    end

    %% Split 91-bit frame payload
    lsf_bits_rx = frameBits_rx(1 : p*lsfBits);
    gain_bits_rx = frameBits_rx(p*lsfBits + 1 : p*lsfBits + gainBits);
    pitch_bits_rx = frameBits_rx(p*lsfBits + gainBits + 1 : p*lsfBits + gainBits + pitchBits);
    voiced_bit_rx = frameBits_rx(end);

    %% Bits -> quantized integers
    lsf_bits_mat = reshape(lsf_bits_rx, lsfBits, p).';
    lsf_q_rx = bi2de(lsf_bits_mat, 'left-msb');

    gain_q_rx = bi2de(gain_bits_rx(:).', 'left-msb');
    pitch_q_rx = bi2de(pitch_bits_rx(:).', 'left-msb');

    %% Dequantization
    LSF_rec(k,:) = lsf_q_rx / (2^lsfBits - 1) * pi;

    if gainMax == gainMin
        Gain_rec(k) = gainMin;
    else
        Gain_rec(k) = gain_q_rx / (2^gainBits - 1) * (gainMax - gainMin) + gainMin;
    end

    pitch_rec(k) = pitch_q_rx / (2^pitchBits - 1) * (pitchMax - pitchMin) + pitchMin;
    voiced_rec(k) = voiced_bit_rx;

    if voiced_rec(k) == 0
        pitch_rec(k) = 0;
    end

    %% Update fallback parameters
    lastLSF = LSF_rec(k,:);
    lastGain = Gain_rec(k);
    lastPitch = pitch_rec(k);
    lastVoiced = voiced_rec(k);
end

%% LPC vocoder synthesis

reconstructedSpeech = zeros(numFrames * frameLen, 1);

for k = 1:numFrames

    currentLSF = LSF_rec(k,:);

    % Safety for lsf2poly
    currentLSF = sort(currentLSF);
    currentLSF = max(0.01, min(pi - 0.01, currentLSF));

    A_rec = lsf2poly(currentLSF);

    if voiced_rec(k) == 1 && pitch_rec(k) > 0

        pitchPeriod = max(1, round(fs / pitch_rec(k)));

        excitation = zeros(frameLen, 1);
        excitation(1:pitchPeriod:end) = 1;

    else

        excitation = randn(frameLen, 1);

    end

    excitation = excitation * Gain_rec(k);

    speechFrame = filter(1, A_rec, excitation);

    idx1 = (k-1)*hopLen + 1; 
    idx2 = idx1 + frameLen - 1;
    
    reconstructedSpeech(idx1:idx2) = speechFrame;
end

%% Normalize, listen, plot, save

reconstructedSpeech = reconstructedSpeech - mean(reconstructedSpeech);
reconstructedSpeech = reconstructedSpeech / (max(abs(reconstructedSpeech)) + eps);

x8_plot = x8(1:min(length(x8), length(reconstructedSpeech)));
rec_plot = reconstructedSpeech(1:min(length(x8), length(reconstructedSpeech)));

fprintf("Reconstructed speech samples = %d\n", length(reconstructedSpeech));


%% Plots for comparison of DSP part
frameAxis = 1:numFrames;

figure;
subplot(2,1,1)
plot(x8_plot), title("Original Prepared Speech x8")
xlabel("Sample Index"), ylabel("Amplitude"), grid on

subplot(2,1,2)
plot(rec_plot), title("Reconstructed Speech")
xlabel("Sample Index"), ylabel("Amplitude"), grid on

figure;
subplot(2,1,2); plot(frameAxis, LSF_rec); title("Recovered LSF"); grid on
xlabel("Frame Index"); ylabel("Recovered LSF Value");
subplot(2,1,1), plot(frameAxis, LSF); title("LSF Parameters");
xlabel("Frame Index"); ylabel("LSF Value");

figure; 
subplot(3,1,1), plot(frameAxis, Gain, frameAxis, Gain_rec, "*r");
title("Gain Values Comparison"); xlabel("Frame Index"); ylabel("Gain");
legend("Before Transmission", "Reconstructed"), grid on

subplot(3,1,2), plot(frameAxis, voicedFlag, frameAxis, voiced_rec, "--r");
title("Voiced/Unvoiced Flag comparison"); xlabel("Frame Index"); 
ylabel("Voiced (1) / Unvoiced (0)"); legend("Before Transmission", "Reconstructed");

subplot(3,1,3), plot(frameAxis, pitch, frameAxis,pitch_rec);
title("Estimated Pitch for Voiced Frames"); legend("Before Transmission", "Reconstructed");
xlabel("Frame Index"); ylabel("Pitch (Hz)"), grid on;

 sound(reconstructedSpeech, fs);

audiowrite("reconstructed_speech.wav", reconstructedSpeech, fs);
fprintf("Saved reconstructed_speech.wav\n");
audiowrite("reconstructed.mp3", reconstructedSpeech, 44100); % mp3 olarak kaydetmek için

%% =========================================================
%% RECONSTRUCTION QUALITY METRICS
%% =========================================================

minLen = min(length(x8), length(reconstructedSpeech));

x_ref = x8(1:minLen);
x_rec = reconstructedSpeech(1:minLen);

%% Overall SNR

noiseSignal = x_ref - x_rec;

SNR_rec = 10 * log10( ...
    sum(x_ref.^2) / (sum(noiseSignal.^2) + eps) ...
);

fprintf("\n--- RECONSTRUCTION QUALITY ---\n");
fprintf("Overall Reconstruction SNR = %.2f dB\n", SNR_rec);

%% Segmental SNR

segFrameLen = 160;
numSeg = floor(minLen / segFrameLen);

segSNR = zeros(numSeg,1);

for k = 1:numSeg

    idx1 = (k-1)*segFrameLen + 1;
    idx2 = idx1 + segFrameLen - 1;

    cleanSeg = x_ref(idx1:idx2);
    recSeg = x_rec(idx1:idx2);

    noiseSeg = cleanSeg - recSeg;

    segSNR(k) = 10 * log10( ...
        sum(cleanSeg.^2) / (sum(noiseSeg.^2) + eps) ...
    );
end

meanSegSNR = mean(segSNR);

fprintf("Mean Segmental SNR = %.2f dB\n", meanSegSNR);



%% NOTLAR

% her framin başına preamble koy
% Her frame + preamble'ı 5 er kez gönder
% en çok bit hangi paketteyse onu seçtirme kullandırma kodu

% kurcalanabilir parametreler
    % frameLen
    % p
    % overlap ???
    % gainler
    % parametre bit sayısı
    % 