package org.jax.bamextractor.processing;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;

import org.jax.bamextractor.util.Location;

import htsjdk.samtools.SAMRecord;

public class ReadCounts implements ATACReadProcessor{

	private BufferedWriter _bw;
	
	public ReadCounts(String outpath) {
		try {
			_bw = new BufferedWriter(new FileWriter(outpath));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	@Override
	public void close() {
		try {
			_bw.flush();
			_bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void processReads(Location l, String referenceseq, LinkedList<SAMRecord> forwardreads, LinkedList<SAMRecord> reversereads) {
		try {
			ATACProcessorFunctions f = new ATACProcessorFunctions();
			f.writeReadCounts(_bw, l, f.getReadPairs(forwardreads, reversereads).size());
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
