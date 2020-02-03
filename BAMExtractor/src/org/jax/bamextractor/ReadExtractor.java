package org.jax.bamextractor;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.TreeMap;

import org.jax.FASTAReader.FASTAReader;
import org.jax.FASTAReader.FASTAReaderFactory;
import org.jax.bamextractor.processing.ATACReadProcessor;
import org.jax.bamextractor.processing.CutMatrices;
import org.jax.bamextractor.processing.CutPileups;
import org.jax.bamextractor.processing.ForwardCutPileups;
import org.jax.bamextractor.processing.ForwardMedianInsert;
import org.jax.bamextractor.processing.InsertPileups;
import org.jax.bamextractor.processing.InsertSizeDistribution;
import org.jax.bamextractor.processing.ReverseCutPileups;
import org.jax.bamextractor.processing.ReverseMedianInsert;
import org.jax.bamextractor.processing.SequenceFrequencies;
import org.jax.bamextractor.util.Location;
import org.jax.bamextractor.util.Util;

public class ReadExtractor {
	
	public static void main(String[] args){
		boolean useduplicateflag = true;
		if(args.length >= 6) {
			
			if(args[args.length-1].equals("--keepduplicates")) {
				useduplicateflag = false;
				System.out.println("Keeping duplicates.");
			}
			
			String bamfile = args[0];
			String locifile = args[1];
			Map<String, String> fastafiles = null;
			try {
				fastafiles = FileMap.getFileMap(args[2]);
			} catch (IOException e1) {
				e1.printStackTrace();
				System.exit(1);
			}
			int charperline = Integer.parseInt(args[3]);
			String outdir = args[4];
			if(outdir.endsWith("/") || outdir.endsWith("\\")) {
				outdir = outdir.substring(0, outdir.length()-1);
			}
			
			String prefix = args[5];
			//Processor in the sense they are processing reads to generate different representations of the data.
			ATACReadProcessor[] processors = new ATACReadProcessor[9];
			processors[0] = new CutMatrices(outdir+"/"+prefix+"_cutmatrices.txt");
			processors[1] = new CutPileups(outdir+"/"+prefix+"_cutpileups.txt");
			processors[2] = new ForwardCutPileups(outdir+"/"+prefix+"_forwardcutpileups.txt");
			processors[3] = new ReverseCutPileups(outdir+"/"+prefix+"_reversecutpileups.txt");
			processors[4] = new ForwardMedianInsert(outdir+"/"+prefix+"_forwardmedianinsert.txt");
			processors[5] = new ReverseMedianInsert(outdir+"/"+prefix+"_reversemedianinsert.txt");
			processors[6] = new InsertPileups(outdir+"/"+prefix+"_insertpileups.txt");
			processors[7] = new InsertSizeDistribution(outdir+"/"+prefix+"_insertsizes.txt");
			processors[8] = new SequenceFrequencies(outdir+"/"+prefix+"_sequencefreq.txt");
			
			int thresh = Integer.MAX_VALUE;
			if(args.length > 6 && !args[5].trim().equals("")) {
				thresh = Integer.parseInt(args[6]);
			}
			else {
				InsertSizeThreshold ist = new InsertSizeThreshold();
				try {
					thresh = ist.getInsertSizeThreshold(bamfile, useduplicateflag);
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}
			
			System.out.println("Insert size threshold: "+thresh);
			
			try {
				ReadExtractor t = new ReadExtractor();
				t.extractFeatures(bamfile, locifile, fastafiles, charperline, processors, thresh, useduplicateflag);
			} catch (IOException e) {
				e.printStackTrace();
			}
			for(int i = 0; i < processors.length; i++) {
				processors[i].close();
			}
		}
		
	}

	public void extractFeatures(String chrbamfile, String peakfile, Map<String, String> fastafiles, int charperline, ATACReadProcessor[] processors, int insertthresh, boolean useduplicateflag) throws IOException{
		
		SamReaderFactory factory = SamReaderFactory.makeDefault()
	              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	              .validationStringency(ValidationStringency.SILENT);
		
		final SamReader reader = factory.open(new File(chrbamfile));
		
		if(!(reader.getFileHeader().getSortOrder().getComparatorInstance() instanceof SAMRecordCoordinateComparator)) {
			System.out.println("The input BAM file must be coordinate sorted.");
			reader.close();
			System.exit(0);
		}
		
		FASTAReader fr = null;
		
		SAMRecordIterator it = reader.iterator();
		it = it.assertSorted(SAMFileHeader.SortOrder.coordinate);
		
		Util util = new Util();
		TreeMap<String, Location[]> m = util.getChrStartSorted((util.readLocations(peakfile)));
		Location[] sortedlocations = null;
		LinkedList<SAMRecord> forwardcurrent = null;
		LinkedList<SAMRecord> forwardnext = null;
		LinkedList<SAMRecord> reversecurrent = null;
		LinkedList<SAMRecord> reversenext = null;
		int cli = 0; //Current location index
		
		int positiveinserts = 0;
		int zeroinserts = 0;
		int negativeinserts = 0;
		int invalidreads = 0;
		int totalreads = 0;
		boolean first = true;
		String curchr = null;
		while(it.hasNext()){
			SAMRecord next = it.next();
			
			if(!fastafiles.containsKey(next.getReferenceName())){
				continue;
			}

			if(first || !curchr.equals(next.getReferenceName())) {
				
				if(!first && !curchr.equals(next.getReferenceName())){
					//process remaining peaks
					while(cli < sortedlocations.length){
						fr.setRegion(sortedlocations[cli].getStart(), sortedlocations[cli].getEnd());
						String refseq = fr.getFASTASequence(sortedlocations[cli].getStart(), sortedlocations[cli].getEnd());
						refseq = fixRefSeq(sortedlocations[cli], refseq);
						for(int i = 0; i < processors.length; i++) {
							processors[i].processReads(sortedlocations[cli], refseq, forwardcurrent, reversecurrent);
						}
						forwardcurrent = forwardnext;
						reversecurrent = reversenext;
						forwardnext = new LinkedList<SAMRecord>();
						reversenext = new LinkedList<SAMRecord>();
						cli++;
						if((cli+1) < sortedlocations.length) {
							Location nextlocation = sortedlocations[cli+1];
							addCurrentToNext(forwardcurrent, forwardnext, nextlocation);
							addCurrentToNext(reversecurrent, reversenext, nextlocation);
						}
					}
				}
				
				curchr = next.getReferenceName();
				first = false;
				forwardcurrent = new LinkedList<SAMRecord>();
				forwardnext = new LinkedList<SAMRecord>();
				reversecurrent = new LinkedList<SAMRecord>();
				reversenext = new LinkedList<SAMRecord>();
				String curfasta = fastafiles.get(curchr);
				fr = FASTAReaderFactory.getForwardFASTAReader(curfasta, charperline);
				sortedlocations = m.get(curchr);
				if(sortedlocations == null) {
					sortedlocations = new Location[0];
				}
				cli = 0;
			}
			
			totalreads++;
			boolean isnegative = next.getReadNegativeStrandFlag();

			if(!isValidRead(next, insertthresh, useduplicateflag)){
				invalidreads++;
				continue;
			}
			
			int insertsize = next.getInferredInsertSize();
			if(isnegative) {
				insertsize = insertsize*-1;
			}
			if(insertsize < 0) {
				negativeinserts++;
				continue;
			}
			
			if(insertsize == 0) {
				zeroinserts++;
				continue;
			}
			positiveinserts++;

			if(cli >= sortedlocations.length) {
				continue;
			}
			
			int start = -1;
			int end = -1;
			int curstart = next.getAlignmentStart();
			if(isnegative) {
				end = next.getAlignmentEnd();
				start = end-insertsize+1;
			}
			else{
				start = next.getAlignmentStart();
				end = start+insertsize-1;
			}
			
			
			Location curlocation = sortedlocations[cli];
			if(curlocation.getStart() <= end && start <= curlocation.getEnd()) {
				if(isnegative) {
					reversecurrent.add(next);
				}
				else {
					forwardcurrent.add(next);
				}
			}
			
			if((cli+1) < sortedlocations.length) {
				Location nextlocation = sortedlocations[cli+1];
				if(nextlocation.getStart() <= end && start <= nextlocation.getEnd()) {
					if(isnegative) {
						reversenext.add(next);
					}
					else {
						forwardnext.add(next);
					}
				}
			}
			

			if(curstart > sortedlocations[cli].getEnd()) {
				fr.setRegion(sortedlocations[cli].getStart(), sortedlocations[cli].getEnd());
				String refseq = fr.getFASTASequence(sortedlocations[cli].getStart(), sortedlocations[cli].getEnd());
				refseq = fixRefSeq(sortedlocations[cli], refseq);
				for(int i = 0; i < processors.length; i++) {
					processors[i].processReads(sortedlocations[cli], refseq, forwardcurrent, reversecurrent);
				}
				forwardcurrent = forwardnext;
				reversecurrent = reversenext;
				forwardnext = new LinkedList<SAMRecord>();
				reversenext = new LinkedList<SAMRecord>();
				cli++;
				//add current peaks to next peaks if they overlap
				if((cli+1) < sortedlocations.length) {
					Location nextlocation = sortedlocations[cli+1];
					addCurrentToNext(forwardcurrent, forwardnext, nextlocation);
					addCurrentToNext(reversecurrent, reversenext, nextlocation);
				}
			}
		}

		//Close out the last peak?
		while(cli < sortedlocations.length){
			fr.setRegion(sortedlocations[cli].getStart(), sortedlocations[cli].getEnd());
			String refseq = fr.getFASTASequence(sortedlocations[cli].getStart(), sortedlocations[cli].getEnd());
			refseq = fixRefSeq(sortedlocations[cli], refseq);
			for(int i = 0; i < processors.length; i++) {
				processors[i].processReads(sortedlocations[cli], refseq, forwardcurrent, reversecurrent);
			}
			forwardcurrent = forwardnext;
			reversecurrent = reversenext;
			forwardnext = new LinkedList<SAMRecord>();
			reversenext = new LinkedList<SAMRecord>();
			cli++;
			if((cli+1) < sortedlocations.length) {
				Location nextlocation = sortedlocations[cli+1];
				addCurrentToNext(forwardcurrent, forwardnext, nextlocation);
				addCurrentToNext(reversecurrent, reversenext, nextlocation);
			}
		}
		

		it.close();
		reader.close();
		fr.close();
		
		System.out.println("Invalid Reads: "+invalidreads);
		System.out.println("Positive Inserts: "+positiveinserts);
		System.out.println("Negative Inserts: "+negativeinserts);
		System.out.println("Zero Inserts: "+zeroinserts);
		System.out.println("Total Reads: "+totalreads);

	}
	
	private String fixRefSeq(Location l, String refseq) {
		int peaklength = l.getEnd()-l.getStart()+1;
		int reflength = refseq.length();
		
		//for instances where the peak is over the end position, add N
		if(reflength < peaklength) {
			StringBuilder sb = new StringBuilder();
			for(int i = reflength; i <= peaklength; i++) {
				sb.append("N");
			}
			return refseq+sb.toString();
		}
		return refseq;
	}
	
	private void addCurrentToNext(LinkedList<SAMRecord> currentlist, LinkedList<SAMRecord> nextlist, Location nextpeak) {
		for(Iterator<SAMRecord> it = currentlist.iterator(); it.hasNext();) {
			SAMRecord next = it.next();
			boolean isnegative = next.getReadNegativeStrandFlag();
			int start = -1;
			int end = -1;
			
			int insertsize = next.getInferredInsertSize();
			if(isnegative) {
				insertsize = insertsize*-1;
				end = next.getAlignmentEnd();
				start = end-insertsize+1;
			}
			else{
				start = next.getAlignmentStart();
				end = start+insertsize-1;
			}
			
			if(nextpeak.getStart() <= end && start <= nextpeak.getEnd()) {
				if(isnegative) {
					nextlist.add(next);
				}
				else {
					nextlist.add(next);
				}
			}
		}
	}
	
	private boolean isValidRead(SAMRecord record, int thresh, boolean useduplicateflag){
		if ((!record.getReadPairedFlag() ||
                record.getReadUnmappedFlag() ||
                record.getMateUnmappedFlag() ||
                record.isSecondaryOrSupplementary() ||
                (record.getDuplicateReadFlag() && useduplicateflag) )) {
			return false;
		}
		
		if(record.getReferenceIndex() != record.getMateReferenceIndex()){
			return false;
		}
		
		if(Math.abs(record.getInferredInsertSize()) > thresh){
			return false;
		}
		return true;
	}
	
}
