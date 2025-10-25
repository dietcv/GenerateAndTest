#pragma once
#include "SC_PlugIn.hpp"
#include <array>
#include <cmath>

namespace Utils {

// ===== CONSTANTS =====

inline constexpr float PI = 3.14159265358979323846f;
inline constexpr float TWO_PI = 6.28318530717958647692f;

// ===== FAST APPROXIMATIONS =====

inline float tanApprox(float x) {
    float x2 = x * x;
    return x * (0.999999492001f + x2 * -0.096524608111f) / 
                (1.0f + x2 * (-0.429867256894f + x2 * 0.009981877999f));
}

// ===== BASIC MATH UTILITIES =====

inline float lerp(float a, float b, float t) {
    return a + t * (b - a);
}

// ===== ONE POLE FILTERS =====

struct OnePoleHz {
    float m_state{0.0f};
   
    float processLowpass(float input, float cutoffHz, float sampleRate) {

        // Clip slope to Nyquist range and take absolute value
        float slope = cutoffHz / sampleRate;
        float safeSlope = std::abs(sc_clip(slope, -0.5f, 0.5f));
        
        // Calculate coefficient: b = exp(-2π * slope)
        float coeff = std::exp(-TWO_PI * safeSlope);
        
        // OnePole formula: y[n] = x[n] * (1-b) + y[n-1] * b
        m_state = input * (1.0f - coeff) + m_state * coeff;
        
        return m_state;
    }
   
    float processHighpass(float input, float cutoffHz, float sampleRate) {
        float lowpassed = processLowpass(input, cutoffHz, sampleRate);
        return input - lowpassed;
    }

    void reset() {
        m_state = 0.0f;
    }
};

struct OnePoleSlope {
    float m_state{0.0f};
      
    float processLowpass(float input, float slope) {

        // Clip slope to Nyquist range and take absolute value
        float safeSlope = std::abs(sc_clip(slope, -0.5f, 0.5f));
       
        // Calculate coefficient: b = exp(-2π * slope)
        float coeff = std::exp(-2.0f * Utils::PI * safeSlope);
       
        // OnePole formula: y[n] = x[n] * (1-b) + y[n-1] * b
        m_state = input * (1.0f - coeff) + m_state * coeff;
       
        return m_state;
    }
   
    float processHighpass(float input, float slope) {
        float lowpassed = processLowpass(input, slope);
        return input - lowpassed;
    }

    void reset() {
        m_state = 0.0f;
    }
};

// ===== ÉMILIE GILLET SVF =====

struct EmilieSVF {
    float y0{0.0f};
    float y1{0.0f};
    
    float m_g{0.0f};
    float m_r{0.0f};
    float m_h{0.0f};
    float m_rpg{0.0f};
    
    void updateCoefficients(float cutoff, float q, float sampleRate) {

        m_g = std::tan(cutoff * Utils::PI / sampleRate);
        m_r = 1.0f / q;
        m_h = 1.0f / (1.0f + m_r * m_g + m_g * m_g);
        m_rpg = m_r + m_g;
        
        // Flush denormals from coefficients
        m_g = zapgremlins(m_g);
        m_r = zapgremlins(m_r);
        m_h = zapgremlins(m_h);
        m_rpg = zapgremlins(m_rpg);
    }
    
    inline float processBandpass(float xin) {
        float hp = (xin - m_rpg * y0 - y1) * m_h;
        float bp = m_g * hp + y0;
        y0 = zapgremlins(m_g * hp + bp);
        float lp = m_g * bp + y1;
        y1 = zapgremlins(m_g * bp + lp);
        
        return bp;
    }
    
    inline float processLowpass(float xin) {
        float hp = (xin - m_rpg * y0 - y1) * m_h;
        float bp = m_g * hp + y0;
        y0 = zapgremlins(m_g * hp + bp);
        float lp = m_g * bp + y1;
        y1 = zapgremlins(m_g * bp + lp);
        
        return lp;
    }
    
    inline float processHighpass(float xin) {
        float hp = (xin - m_rpg * y0 - y1) * m_h;
        float bp = m_g * hp + y0;
        y0 = zapgremlins(m_g * hp + bp);
        float lp = m_g * bp + y1;
        y1 = zapgremlins(m_g * bp + lp);
        
        return hp;
    }
    
    void reset() {
        y0 = y1 = 0.0f;
    }
};

// ===== BIQUAD COEFFICIENTS =====

struct BiquadCoefficients {
    float a1, a2;
    float b0, b2;  // For bandpass (b1 always 0)
    
    // RBJ Audio EQ Cookbook allpass formula
    static BiquadCoefficients allpass(float freq, float q, float sampleRate) {
        BiquadCoefficients c;
        
        float w0 = Utils::TWO_PI * freq / sampleRate;
        float cosw0 = std::cos(w0);
        float sinw0 = std::sin(w0);
        float alpha = sinw0 / (2.0f * q);
        
        float a0 = 1.0f + alpha;
        
        c.a1 = -2.0f * cosw0 / a0;
        c.a2 = (1.0f - alpha) / a0;
        c.b0 = 0.0f;  // Not used for allpass
        c.b2 = 0.0f;  // Not used for allpass
        
        return c;
    }

    // RBJ Audio EQ Cookbook bandpass (constant skirt gain, peak gain = Q)
    static BiquadCoefficients bandpass(float freq, float q, float sampleRate) {
        BiquadCoefficients c;
        
        float w0 = Utils::TWO_PI * freq / sampleRate;
        float cosw0 = std::cos(w0);
        float sinw0 = std::sin(w0);
        float alpha = sinw0 / (2.0f * q);
        
        float a0 = 1.0f + alpha;
        
        c.a1 = -2.0f * cosw0 / a0;
        c.a2 = (1.0f - alpha) / a0;
        c.b0 = (sinw0 / 2.0f) / a0;   // Constant skirt gain: peak gain = Q
        c.b2 = -(sinw0 / 2.0f) / a0;
        
        return c;
    }
};

// ===== BIQUAD ALLPASS FILTERS =====

// Direct Form I (4 states)
struct BiquadAllpass_DF1 {
    float x1{0.0f}, x2{0.0f};
    float y1{0.0f}, y2{0.0f};
    
    inline float process(float x0, const BiquadCoefficients& c) {
        // Allpass: y[n] = a2*x[n] + a1*x[n-1] + x[n-2] - a1*y[n-1] - a2*y[n-2]
        float y0 = c.a2 * x0 + c.a1 * x1 + x2 - c.a1 * y1 - c.a2 * y2;
        
        y2 = zapgremlins(y1);
        x2 = zapgremlins(x1);
        x1 = zapgremlins(x0);
        y1 = zapgremlins(y0);
        
        return y0;
    }
    
    void reset() {
        x1 = x2 = y1 = y2 = 0.0f;
    }
};

// Transposed Direct Form II (2 states)
struct BiquadAllpass_TDF2 {
    float z1{0.0f}, z2{0.0f};
    
    inline float process(float x, const BiquadCoefficients& c) {
        float y = c.a2 * x + z1;
        z1 = c.a1 * x - c.a1 * y + z2;
        z2 = x - c.a2 * y;
        
        z1 = zapgremlins(z1);
        z2 = zapgremlins(z2);
        
        return y;
    }
    
    void reset() {
        z1 = z2 = 0.0f;
    }
};

// ===== BIQUAD BANDPASS FILTERS =====

// Direct Form I (4 states)
struct BiquadBandpass_DF1 {
    float x1{0.0f}, x2{0.0f};
    float y1{0.0f}, y2{0.0f};
    
    inline float process(float x0, const BiquadCoefficients& c) {
        // BPF: y[n] = b0*x[n] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
        // (b1 = 0, so x[n-1] term is omitted)
        float y0 = c.b0 * x0 + c.b2 * x2 - c.a1 * y1 - c.a2 * y2;
        
        y2 = zapgremlins(y1);
        x2 = zapgremlins(x1);
        x1 = zapgremlins(x0);
        y1 = zapgremlins(y0);
        
        return y0;
    }
    
    void reset() {
        x1 = x2 = y1 = y2 = 0.0f;
    }
};

// Transposed Direct Form II (2 states)
struct BiquadBandpass_TDF2 {
    float z1{0.0f}, z2{0.0f};
    
    inline float process(float x, const BiquadCoefficients& c) {
        float y = c.b0 * x + z1;
        z1 = -c.a1 * y + z2;  // b1 = 0, so b1*x term omitted
        z2 = c.b2 * x - c.a2 * y;
        
        z1 = zapgremlins(z1);
        z2 = zapgremlins(z2);
        
        return y;
    }
    
    void reset() {
        z1 = z2 = 0.0f;
    }
};

// ===== DISPERSER =====

template<int NumAllpasses>
struct Disperser {
    std::array<BiquadAllpass_TDF2, NumAllpasses> allpasses;

    Disperser() = default;

    // Process audio through cascaded allpass filters
    inline float process(float input, float freq, float resonance, float sampleRate) {

        // Convert resonance (0 - 1) to Q (1.0 - 0.01) (inverted for allpass) 
        float q = sc_clip(1.0f - resonance, 0.01f, 1.0f);
        
        // Calculate coefficients
        auto coeffs = BiquadCoefficients::allpass(freq, q, sampleRate);
        
        // Cascade allpass filters
        float processed = input;
        for (int i = 0; i < NumAllpasses; ++i) {
            processed = allpasses[i].process(processed, coeffs);
        }
        
        return processed;
    }
};

// ===== FILTERBANK =====

template<int NumBands>
struct FilterBank {
    std::array<EmilieSVF, NumBands> filters;
    std::array<float, NumBands> frequencies;

    FilterBank() = default;

    // Calculate frequency distribution
    void calculateFrequencies(float baseFreq, float spread, float warp) {
        
        for (int i = 0; i < NumBands; ++i) {

            // 1-based indexing like Fors Opal
            float x = static_cast<float>(i + 1);  
            float len = static_cast<float>(NumBands);
            
            // Fors Opal formula: x = pow(x, spread * pow(len / x, warp))
            float exponent = spread * std::pow(len / x, warp);
            float multiplier = std::pow(x, exponent);
            
            frequencies[i] = baseFreq * multiplier;
        }
    }
    
    // Process audio through filter bank
    inline float process(float input, float baseFreq, float spread, float warp, 
                        float resonance, float sampleRate) {
        
        // Convert resonance (0 - 1) to Q (0.5 - 25.0)
        float q = 0.5f + std::sqrt(resonance) * 24.5f;
        
        // Calculate frequency distribution
        calculateFrequencies(baseFreq, spread, warp);

        // Pre-calculate frequency limits
        const float nyquistLimit = sampleRate * 0.49f;
        const float fadeWidth = 2000.0f;
        const float fadeStart = nyquistLimit - fadeWidth;
                
        float output = 0.0f;
        
        // Update and process each filter
        for (int i = 0; i < NumBands; ++i) {

            float freq = frequencies[i];
            
            // Check if frequency is valid and update filter
            if (freq > 20.0f && freq < nyquistLimit) {
                filters[i].updateCoefficients(freq, q, sampleRate);
                
                float fadeOut = 1.0f;
                if (freq > fadeStart) {
                    fadeOut = 1.0f - (freq - fadeStart) / fadeWidth;
                }
                
                float filtered = filters[i].processBandpass(input);
                output += filtered * fadeOut;
            }
        }
        
        return output * (1.0f / std::sqrt(static_cast<float>(NumBands)));
    }
};

// ===== SUPERSAW OSCILLATOR =====

struct SuperSawOsc {
    // Individual oscillator states for the 7 sawtooth waves
    struct SawState {
        float m_phase{0.0f};
        
        float next(float phaseInc) {
            // Update internal phase (standard 0-1 range)
            m_phase += phaseInc;
            
            // Wrap phase to [0, 1] range
            if (m_phase >= 1.0f)
                m_phase -= 1.0f;
            else if (m_phase < 0.0f)
                m_phase += 1.0f;
            
            // Generate sawtooth output: ramp from -1 to +1
            return (m_phase * 2.0f) - 1.0f;
        }

        void reset(bool isCenter, RGen& rgen) {
            if (isCenter) {
                m_phase = 0.0f;  // Center oscillator always starts at phase 0
            } else {
                // Side oscillators get random phase
                m_phase = rgen.frand() * 2.0f - 1.0f;
            }
        }
    };
    
    // 7 individual sawtooth oscillators (1 center + 6 sides)
    SawState m_centerOsc;
    std::array<SawState, 6> m_sideOscs;
    
    // JP-8000 polynomial curves for detune and gain compensation
    float detuneCurve(float x) const {
        return (10028.7312891634f * static_cast<float>(std::pow(x, 11))) -
            (50818.8652045924f * static_cast<float>(std::pow(x, 10))) +
            (111363.4808729368f * static_cast<float>(std::pow(x, 9))) -
            (138150.6761080548f * static_cast<float>(std::pow(x, 8))) +
            (106649.6679158292f * static_cast<float>(std::pow(x, 7))) -
            (53046.9642751875f * static_cast<float>(std::pow(x, 6))) +
            (17019.9518580080f * static_cast<float>(std::pow(x, 5))) -
            (3425.0836591318f * static_cast<float>(std::pow(x, 4))) +
            (404.2703938388f * static_cast<float>(std::pow(x, 3))) -
            (24.1878824391f * static_cast<float>(std::pow(x, 2))) +
            (0.6717417634f * x) +
            0.0030115596f;
    }
    
    float centerGain(float x) const {
        return (-0.55366f * x) + 0.99785f;
    }
    
    float sideGain(float x) const {
        return (-0.73764f * x * x) + (1.2841f * x) + 0.044372f;
    }
    
    void reset(RGen& rgen) {
        m_centerOsc.reset(true, rgen);   // Center: fixed phase at 0
        for (auto& osc : m_sideOscs) {
            osc.reset(false, rgen);      // Sides: random phase offsets
        }
    }
    
    float process(float slope, float mix, float detune) {

        // Calculate detune ratios using JP-8000 curve
        const float detuneAmount = detuneCurve(detune);
        const std::array<float, 6> detuneRatios = {
            -0.11002313f, 
            -0.06288439f,
            -0.01952356f,
            0.01991221f,  
            0.06216538f,
            0.10745242f
        };
        
        // Generate center oscillator
        float centerOut = m_centerOsc.next(slope);
        
        // Generate side oscillators with detune
        float sideOut = 0.0f;
        for (int i = 0; i < 6; ++i) {
            float detuneSlope = slope * (1.0f + (detuneAmount * detuneRatios[i]));
            sideOut += m_sideOscs[i].next(detuneSlope);
        }
        
        // Apply JP-8000 gain compensation curves
        float centerGained = centerOut * centerGain(mix);
        float sideGained = sideOut * sideGain(mix);
        
        // Sum all oscillators
        return centerGained + sideGained;
    }
};

} // namespace Utils