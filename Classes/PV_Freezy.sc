PV_Freezish : PV_ChainUGen {
	*new { |buffer, atkCoeff = 0, dcyCoeff = 0|
		^this.multiNew('control', buffer, atkCoeff, dcyCoeff)
	}

	// exp(log001 / numDecay)
	// numDecay = time * (sr / framesize / hop)
	// exp(log001 / (time * sr / framesize / hop)))
	// exp(log001 / time * sd * framesize * hop))
	*lag { |buffer, atkTime = 0.01, dcyTime = 0.01, hop = 0.5|
		var sdur = SampleDur.ir;
		var frames = if(buffer.rate == \control) {
			Latch.kr(BufFrames.kr(buffer), buffer >= 0)
		} {
			buffer.numFrames
		};
		var factor = log(0.001) * hop * frames * sdur;
		var coeff = { |time|
			exp(factor / max(0.001, time))
		};
		^this.multiNew('control', buffer, coeff.(atkTime), coeff.(dcyTime))
	}
}
