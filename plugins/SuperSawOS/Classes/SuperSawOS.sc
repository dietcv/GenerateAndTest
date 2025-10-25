SuperSawOS : UGen {
    *ar { |freq = 440, mix, detune, oversample = 0|
		if(freq.rate != 'audio') { freq = K2A.ar(freq) };
        ^this.multiNew('audio', freq, mix, detune, oversample)
    }
}
