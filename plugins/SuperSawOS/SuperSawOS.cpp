#include "SC_PlugIn.hpp"
#include "SuperSawOS.hpp"

static InterfaceTable *ft;

SuperSawOS::SuperSawOS() : m_sampleRate(static_cast<float>(sampleRate()))
{
    RGen& rgen = *mParent->mRGen;
    m_oversampling.reset(m_sampleRate);
    m_superSaw.reset(rgen);
    mCalcFunc = make_calc_function<SuperSawOS, &SuperSawOS::next_aa>();
    next_aa(1);
}

void SuperSawOS::next_aa(int nSamples) {

    // Audio-rate parameters
    const float* freqIn = in(Freq);

    // Control-rate parameters
    const float mix = sc_clip(in0(Mix), 0.0f, 1.0f);       
    const float detune = sc_clip(in0(Detune), 0.0f, 1.0f); 
    const int oversampleIndex = sc_clip(static_cast<int>(in0(Oversample)), 0, 4);
    
    // Output buffer
    float* outbuf = out(Out);
    
    // Set oversampling factor
    m_oversampling.setOversamplingIndex(oversampleIndex);
    
    if (oversampleIndex == 0) {
        // No oversampling - direct processing
        for (int i = 0; i < nSamples; ++i) {
            float slope = freqIn[i] / m_sampleRate;
            
            // Process SuperSaw
            float supersawOut = m_superSaw.process(slope, mix, detune);
            
            // Apply pitch-tracking highpass filter
            outbuf[i] = m_pitchTrackingHPF.processHighpass(supersawOut, slope);
        }
    } else {
        // Oversampling enabled
        const int osRatio = m_oversampling.getOversamplingRatio();
        const float osRate = m_sampleRate * osRatio;
        
        for (int i = 0; i < nSamples; ++i) {
            // Prepare oversampling buffer
            m_oversampling.upsample(0.0f);
            float* osBuffer = m_oversampling.getOSBuffer();
            
            for (int k = 0; k < osRatio; k++) {
                float osSlope = freqIn[i] / osRate;
                
                // Process SuperSaw at oversampled rate
                float supersawOut = m_superSaw.process(osSlope, mix, detune);
                
                // Apply pitch-tracking highpass filter
                osBuffer[k] = m_pitchTrackingHPF.processHighpass(supersawOut, osSlope);
            }
            
            outbuf[i] = m_oversampling.downsample();
        }
    }
}

PluginLoad(SuperSawOSUGens)
{
    ft = inTable;
    registerUnit<SuperSawOS>(ft, "SuperSawOS", false);
}