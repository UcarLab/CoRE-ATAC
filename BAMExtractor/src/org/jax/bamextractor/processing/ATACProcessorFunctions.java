package org.jax.bamextractor.processing;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import org.jax.bamextractor.util.Location;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ATACProcessorFunctions {

	public void writeBaseFrequencies(BufferedWriter bw, Location l, String refseq, double[][] freqmat) throws IOException {
		bw.write(l.getChr()+"\t"+l.getStart()+"\t"+l.getEnd()+"\t"+l.getId()+"\n");
		bw.write(refseq+"\n");
		
		for(int i = 0; i < 4; i++) {
			bw.write(Double.toString(freqmat[i][0]));
			for(int j = 1; j < freqmat.length; j++) {
				bw.write(","+freqmat[j][i]);
			}
			bw.write("\n");
		}
	}
	
	public void writeReadCounts(BufferedWriter bw, Location l, int count) throws IOException {
		bw.write(l.getChr()+"\t"+l.getStart()+"\t"+l.getEnd()+"\t"+Integer.toString(count)+"\n");
	}
	
	
	public void writePileups(BufferedWriter bw, Location l, int[] cuts) throws IOException {
		bw.write(l.getChr()+"\t"+l.getStart()+"\t"+l.getEnd()+"\t"+l.getId()+"\n");
		bw.write(Integer.toString(cuts[0]));
		for(int i = 1; i < cuts.length; i++) {
			bw.write(","+cuts[i]);
		}
		bw.write("\n");
	}
	
	public void writeMedianInserts(BufferedWriter bw, Location l, double[] cuts) throws IOException {
		bw.write(l.getChr()+"\t"+l.getStart()+"\t"+l.getEnd()+"\t"+l.getId()+"\n");
		bw.write(Double.toString(cuts[0]));
		for(int i = 1; i < cuts.length; i++) {
			bw.write(","+cuts[i]);
		}
		bw.write("\n");
	}
	
	
	public void writeCutMatrixList(BufferedWriter bw, Location l, int[][] matrix) throws IOException {
		int ml = matrix.length;
		int totalcuts = 0;
		for(int i = 0; i < ml; i++) {
			for(int j = (i+1); j < ml; j++) {
				if(matrix[i][j] > 0) {
					totalcuts++;
				}
			}
		}
		
		String metadata = l.getChr()+"\t"+l.getStart()+"\t"+l.getEnd()+"\t"+l.getId()+"\t"+ml;
		bw.write(metadata+"\n");
		bw.write(totalcuts+"\n");
		
		for(int i = 0; i < ml; i++) {
			for(int j = (i+1); j < ml; j++) {
				if(matrix[i][j] > 0) {
					bw.write(i+","+j+"\n");
				}
			}
		}
		
	}
	
	public void writeCutMatrix(FileOutputStream bw, Location l, int[][] matrix) throws IOException {
		int mlen = matrix.length;
		
		int numentries = matrix.length*matrix.length;
		int trailingbits = numentries%8;
		int numbytes = numentries/8;
		
		ByteBuffer bufsize = ByteBuffer.allocate(4);
		String metadata = l.getChr()+"\t"+l.getStart()+"\t"+l.getEnd()+"\t"+l.getId()+"\t"+mlen;
		byte[] mbytes = metadata.getBytes();
		bw.write(bufsize.putInt(0,mbytes.length).array());
		bw.write(mbytes);
		
		byte cur = 0;
		int position = 0;
		int sum = 0;
		byte[] allbytes = new byte[numbytes+(trailingbits == 0 ? 0 : 1)];
		int idx = 0;
		for(int i = 0; i < matrix.length; i++) {
			for(int j = 0; j < matrix.length; j++) {
				if(matrix[i][j] > 0) {
					sum++;
					cur = (byte) (cur | (1 << position));
				}
				position++;
				if(position == 8) {
					allbytes[idx++] = cur;
					position = 0;
					cur = 0;
				}
			}
		}
		if(position != 0) {
			allbytes[idx++] = cur;
		}
		
		bw.write(allbytes);
		bw.write(bufsize.putInt(0,sum).array());
	}
	
	public void writeInsertSizeDistribution(BufferedWriter bw, Location l, LinkedList<SAMRecord> reads) throws IOException {
		bw.write(l.getChr()+"\t"+l.getStart()+"\t"+l.getEnd()+"\t"+l.getId()+"\n");
		
		for(Iterator<SAMRecord> it = reads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			bw.write(next.getInferredInsertSize()+"\n");
		}
	}
	
	public int[][] getCutMatrix(Location l, LinkedList<ReadPair> pairs){
		int s = l.getStart();
		int len = l.getEnd()-s+1;
		int[][] rv = new int[len][len];
		for(Iterator<ReadPair> it = pairs.iterator(); it.hasNext();) {
			ReadPair next = it.next();
			int p1 = next.getForwardRead().getStart()-s;
			int p2 = next.getReverseRead().getEnd()-s;
			if(p1 > -1 && p1 < len && p2 > -1 && p2 < len) {
				rv[p1][p2] += 1;
			}
		}
		return rv;
	}
	
	public int[] getCutPileupForward(Location l, LinkedList<SAMRecord> reads){
		int s = l.getStart();
		int len = l.getEnd()-s+1;
		int[] rv = new int[len];
		
		for(Iterator<SAMRecord> it = reads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			int p = next.getStart()-s;
			if(p > -1 && p < len) {
				rv[p] += 1;
			}
		}
		
		return rv;
	}
	
	public int[] getCutPileupReverse(Location l, LinkedList<SAMRecord> reads){
		int s = l.getStart();
		int len = l.getEnd()-s+1;
		int[] rv = new int[len];
		
		for(Iterator<SAMRecord> it = reads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			int p = next.getEnd()-s;
			if(p > -1 && p < len) {
				rv[p] += 1;
			}
		}
		
		return rv;
	}
	
	public int[] getCutPileup(Location l, LinkedList<SAMRecord> forwardreads, LinkedList<SAMRecord> reversereads){
		int s = l.getStart();
		int len = l.getEnd()-s+1;
		int[] rv = new int[len];
		
		for(Iterator<SAMRecord> it = forwardreads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			int p = next.getStart()-s;
			if(p > -1 && p < len) {
				rv[p] += 1;
			}
		}
		
		for(Iterator<SAMRecord> it = reversereads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			int p = next.getEnd()-s;
			if(p > -1 && p < len) {
				rv[p] += 1;
			}
		}
		return rv;
	}
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public double[] getReverseMedianInserts(Location l,LinkedList<SAMRecord> reads){
		int s = l.getStart();
		int len = l.getEnd()-s+1;
		double[] rv = new double[len];
		//NOTE: The method implemented here only works because it is assumed that duplicate paired end reads are removed.
		//Otherwise need to make a linked, sort the list, and pull the median value.
		LinkedList[] mtree = new LinkedList[len]; 
		for(int i = 0; i < mtree.length; i++) {
			mtree[i] = new LinkedList<Integer>();
		}
		for(Iterator<SAMRecord> it = reads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			int p = next.getEnd()-s;
			if(p > -1 && p < len) {
				mtree[p].add(Math.abs(next.getInferredInsertSize()));
			}
		}
		
		for(int i = 0; i < mtree.length; i++) {
			Collections.sort((LinkedList<Integer>)mtree[i]);
			Integer[] pvals = ((LinkedList<Integer>)mtree[i]).toArray(new Integer[0]);
			if(pvals.length % 2 == 0) {
				if(pvals.length > 0) {
					rv[i] = (pvals[pvals.length/2]+pvals[(pvals.length/2)-1])/2.0;
				}
				else {
					rv[i] = 0;
				}
			}
			else {
				rv[i] = pvals[pvals.length/2];
			}

		}
		return rv;
	}
	
	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public double[] getForwardMedianInserts(Location l,LinkedList<SAMRecord> reads){
		int s = l.getStart();
		int len = l.getEnd()-s+1;
		double[] rv = new double[len];
		//NOTE: The method implemented here only works because it is assumed that duplicate paired end reads are removed.
		//Otherwise need to make a linked, sort the list, and pull the median value.
		LinkedList[] mtree = new LinkedList[len]; 
		for(int i = 0; i < mtree.length; i++) {
			mtree[i] = new LinkedList<Integer>();
		}
		for(Iterator<SAMRecord> it = reads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			int p = next.getStart()-s;
			if(p > -1 && p < len) {
				mtree[p].add(Math.abs(next.getInferredInsertSize()));
			}
		}
		
		for(int i = 0; i < mtree.length; i++) {
			Collections.sort((LinkedList<Integer>)mtree[i]);
			Integer[] pvals = ((LinkedList<Integer>)mtree[i]).toArray(new Integer[0]);
			if(pvals.length % 2 == 0) {
				if(pvals.length > 0) {
					rv[i] = (pvals[pvals.length/2]+pvals[(pvals.length/2)-1])/2.0;
				}
				else {
					rv[i] = 0;
				}
			}
			else {
				rv[i] = pvals[pvals.length/2];
			}

		}
		return rv;
	}
	

	
	
	
	public int[] getInsertPileup(Location l, LinkedList<ReadPair> pairs){
		int s = l.getStart();
		int len = l.getEnd()-s+1;
		int[] rv = new int[len];
		
		for(Iterator<ReadPair> it = pairs.iterator(); it.hasNext();) {
			ReadPair next = it.next();
			int start = Math.max(next.getForwardRead().getStart()-s, 0);
			int end = Math.max(Math.min(next.getReverseRead().getEnd()-s, len-1), start);
			for(int i = start; i <= end; i++) {
				rv[i] += 1;
			}
		}
		
		return rv;
	}
	

	public double[][] getBaseFrequency(Location l, LinkedList<SAMRecord> forwardreads, LinkedList<SAMRecord> reversereads){
		int s = l.getStart();
		int len = l.getEnd()-s+1;
		double[][] rv = new double[len][4];
		setBaseFrequencies(rv, s, len, forwardreads);
		setBaseFrequencies(rv, s, len, reversereads);			
		return rv;
	}
	
	private void setBaseFrequencies(double[][] fmat, int s, int len, LinkedList<SAMRecord> reads) {
		for(Iterator<SAMRecord> it = reads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			
			int fmatpos = next.getAlignmentStart()-s;
			int basepos = 0;
			char[] sequence = (new String(next.getReadBases())).toCharArray();
			
			Cigar c = next.getCigar();
			List<CigarElement> cel = c.getCigarElements();
			for(Iterator<CigarElement> it2 = cel.iterator(); it2.hasNext();) {
				CigarElement ce = it2.next();
				CigarOperator op = ce.getOperator();
				int l = ce.getLength();
				if(op.consumesReadBases() && op.consumesReferenceBases()) {
					for(int i = basepos; i < basepos+l; i++) {
						int curfmatpos = fmatpos++;
						if(curfmatpos > -1 && curfmatpos < len) {
							setFrequencyMatBase(fmat, curfmatpos, sequence[i]);
						}
					}
					basepos += l;
				}
				else if(op.consumesReadBases()) {
					basepos += l;
				}
				else if(op.consumesReferenceBases()) {
					fmatpos += l;
				}
			}
		}
	}
	
	private void setFrequencyMatBase(double[][] fmat, int curfmatpos, char base) {
		int[] intbases = BaseToInt.getIntBases(base);
		if(intbases != null) {
			double value = (1.0)/intbases.length;
			for(int i = 0; i < intbases.length; i++) {
				fmat[curfmatpos][intbases[i]] += value;
			}
		}
	}
	
	public LinkedList<ReadPair> getReadPairs(LinkedList<SAMRecord> forwardreads, LinkedList<SAMRecord> reversereads){
		LinkedList<ReadPair> rv = new LinkedList<ReadPair>();
		TreeMap<String, SAMRecord> frmap = new TreeMap<String, SAMRecord>();
		for(Iterator<SAMRecord> it = forwardreads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			frmap.put(next.getReadName(), next);
		}
		for(Iterator<SAMRecord> it = reversereads.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			if(frmap.containsKey(next.getReadName())) {
				rv.add(new ReadPair(frmap.get(next.getReadName()), next));
			}
		}
		return rv;
	}
}
