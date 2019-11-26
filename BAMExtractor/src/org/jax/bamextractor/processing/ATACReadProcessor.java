package org.jax.bamextractor.processing;

import java.util.LinkedList;

import org.jax.bamextractor.util.Location;

import htsjdk.samtools.SAMRecord;

public interface ATACReadProcessor {

	public void processReads(Location l, String referenceseq, LinkedList<SAMRecord> forwardreads, LinkedList<SAMRecord> reversereads);
	
	public void close();
}
