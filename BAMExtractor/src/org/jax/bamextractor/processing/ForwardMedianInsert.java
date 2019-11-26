package org.jax.bamextractor.processing;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;

import org.jax.bamextractor.util.Location;

import htsjdk.samtools.SAMRecord;

public class ForwardMedianInsert implements ATACReadProcessor{

	private BufferedWriter _bw;
	
	public ForwardMedianInsert(String outpath) {
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
			f.writeMedianInserts(_bw, l, f.getForwardMedianInserts(l, forwardreads));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
