function [d_preprocessed,sfx_preprocessed] = preprocess_clip(d,sfx)

    %zero center the data
    choffset = median(d, 1, 'omitnan');
    choffset(isnan(choffset)) = 0; %avoid nans
    d_zero = bsxfun(@minus, d, choffset);

    %downsample
    freqs = [1 255]/(sfx/2);
    %anti-aliasing filter
    [b,a]=butter(2,freqs,'bandpass');
    d_anti = zeros(size(d_zero));
    for i=1:size(d_zero,2)
        d_anti(:,i)=filtfilt(b,a,d_zero(:,i));
    end
    d_down = zeros(size(d_anti));
    sfx_preprocessed = 512;
    for i=1:size(d_anti,2)
       d_down(:,i) = resample(d_anti(:,i), sfx_preprocessed, round(sfx));
    end
    %notch filter
    notch_freq = 60; %base freq
    nns = ~isnan(d_down(1,:)); %is this necessary?
    d_preprocessed = d_down; 
    while notch_freq < sfx_preprocessed/2
        [b,a]=fir2(1000,[0 notch_freq-1 notch_freq-.5 notch_freq+.5 notch_freq+1 sfx/2]/(sfx/2),[1 1 0 0 1 1 ]);
        d_preprocessed(:, nns) = filtfilt(b,a,d_preprocessed(:, nns));
        notch_freq = notch_freq + 60;
    end
end