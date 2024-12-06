PV_Freezish : PV_ChainUGen {
	*new { |buffer, freeze = 0|
		^this.multiNew('control', buffer, freeze)
	}
}

PV_Lag : PV_ChainUGen {
	*new { |buffer, attack = 0.01, decay = 1, hop = 0.5|
		^this.multiNew('control', buffer,
			attack, decay,
			SampleRate.ir / BufFrames.kr(buffer) / hop
		)
	}
}
