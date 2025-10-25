FilterBank : UGen {
	*ar { |input, freq, spread, warp, q|
		^this.multiNew('audio', input, freq, spread, warp, q)
	}

	checkInputs {
		^this.checkValidInputs
	}
}
