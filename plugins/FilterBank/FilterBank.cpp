#include "SC_PlugIn.hpp"
#include "FilterBank.hpp"

static InterfaceTable *ft;

FilterBank::FilterBank() : m_sampleRate(static_cast<float>(sampleRate()))
{
    // Initialize parameter cache
    freqPast = in0(Freq);
    spreadPast = in0(Spread);
    warpPast = in0(Warp);
    resonancePast = in0(Resonance);
   
    mCalcFunc = make_calc_function<FilterBank, &FilterBank::next>();
    next(1);
}

void FilterBank::next(int nSamples) {
    // Audio-rate parameters
    const float* input = in(Input);
   
    // Control-rate parameters with smooth interpolation
    auto slopedFreq = makeSlope(sc_clip(in0(Freq), 20.0f, m_sampleRate * 0.49f), freqPast);
    auto slopedSpread = makeSlope(sc_clip(in0(Spread), 0.0f, 2.0f), spreadPast);
    auto slopedWarp = makeSlope(sc_clip(in0(Warp), -1.0f, 1.0f), warpPast);
    auto slopedResonance = makeSlope(sc_clip(in0(Resonance), 0.0f, 0.99f), resonancePast);
       
    // Output pointers
    float* outbuf = out(Out);
   
    // Process audio
    for (int i = 0; i < nSamples; ++i) {
        // Process filter bank
        outbuf[i] = filterBank.process(
            input[i],
            slopedFreq.consume(),
            slopedSpread.consume(),
            slopedWarp.consume(),
            slopedResonance.consume(),
            m_sampleRate
        );
    }
   
    // Update parameter cache
    freqPast = slopedFreq.value;
    spreadPast = slopedSpread.value;
    warpPast = slopedWarp.value;
    resonancePast = slopedResonance.value;
}

PluginLoad(FilterBankUGens) {
    ft = inTable;
    registerUnit<FilterBank>(ft, "FilterBank", false);
}