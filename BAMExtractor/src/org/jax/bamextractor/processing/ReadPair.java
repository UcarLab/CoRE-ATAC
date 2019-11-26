package org.jax.bamextractor.processing;

import htsjdk.samtools.SAMRecord;

public class ReadPair {

	private SAMRecord _fr;
	private SAMRecord _rr;
	
	public ReadPair(SAMRecord fr, SAMRecord rr) {
		_fr = fr;
		_rr = rr;
	}
	
	public SAMRecord getForwardRead() {
		return _fr;
	}
	
	public SAMRecord getReverseRead() {
		return _rr;
	}
	
}
