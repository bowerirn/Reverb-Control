function [panel, ref, reverb] = clipped_irs(max_len, source_to_ref_cm)
    [panel, ref, reverb, fs] = irs();
    source_to_ref_samples = source_to_ref_cm * fs / (34300); % speed of sound in cm/s

    [~, source_idx] = max(ref);
    start = source_idx - round(source_to_ref_samples);
    cut = start + max_len - 1;

    panel = panel(start:cut);
    ref = ref(start:cut);
    reverb = reverb(start:cut);

end