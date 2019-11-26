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
import java.util.TreeMap;

import org.jax.bamextractor.processing.ATACReadProcessor;
import org.jax.bamextractor.processing.ReadCounts;
import org.jax.bamextractor.util.Location;
import org.jax.bamextractor.util.Util;

public class ReadCountExtractor {
	
	
	//TODO make the code reusable instead of copying, only copying due to time limitations
	
	public static void main(String[] args){
		
		if(args.length >= 3) {
			
			String bamfile = args[0];
			String locifile = args[1];
			String outfile = args[2];
			
			//Processor in the sense they are processing reads to generate different representations of the data.
			ATACReadProcessor[] processors = new ATACReadProcessor[1];
			processors[0] = new ReadCounts(outfile);

			
			int thresh = Integer.MAX_VALUE;
			if(args.length > 3) {
				thresh = Integer.parseInt(args[3]);
			}
			else {
				InsertSizeThreshold ist = new InsertSizeThreshold();
				try {
					thresh = ist.getInsertSizeThreshold(bamfile);
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(1);
				}
			}
			
			System.out.println("Insert size threshold: "+thresh);
			
			try {
				ReadCountExtractor t = new ReadCountExtractor();
				t.extractReadCounts(bamfile, locifile, processors, thresh);
			} catch (IOException e) {
				e.printStackTrace();
			}
			for(int i = 0; i < processors.length; i++) {
				processors[i].close();
			}
		}
		else {
			System.out.println("Usage: ReadCountExtractor.jar <bamfile> <peakfile> <outfile> <insertsize threshold>");
		}
		
	}

	public void extractReadCounts(String chrbamfile, String peakfile, ATACReadProcessor[] processors, int insertthresh) throws IOException{
		
		SamReaderFactory factory = SamReaderFactory.makeDefault()
	              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	              .validationStringency(ValidationStringency.SILENT);
		
		final SamReader reader = factory.open(new File(chrbamfile));
		
		if(!(reader.getFileHeader().getSortOrder().getComparatorInstance() instanceof SAMRecordCoordinateComparator)) {
			System.out.println("The input BAM file must be coordinate sorted.");
			reader.close();
			System.exit(0);
		}
		
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

			if(first || !curchr.equals(next.getReferenceName())) {
				
				if(!first && !curchr.equals(next.getReferenceName())){
					//process remaining peaks
					while(cli < sortedlocations.length){
						String refseq = "";
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
				sortedlocations = m.get(curchr);
				if(sortedlocations == null) {
					sortedlocations = new Location[0];
				}
				cli = 0;
			}
			
			totalreads++;
			boolean isnegative = next.getReadNegativeStrandFlag();

			if(!isValidRead(next, insertthresh)){
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
				String refseq = "";
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
			String refseq = "";
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
		
		System.out.println("Invalid Reads: "+invalidreads);
		System.out.println("Positive Inserts: "+positiveinserts);
		System.out.println("Negative Inserts: "+negativeinserts);
		System.out.println("Zero Inserts: "+zeroinserts);
		System.out.println("Total Reads: "+totalreads);

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
	
	private boolean isValidRead(SAMRecord record, int thresh){
		if ((!record.getReadPairedFlag() ||
                record.getReadUnmappedFlag() ||
                record.getMateUnmappedFlag() ||
                record.isSecondaryOrSupplementary() ||
                record.getDuplicateReadFlag())) {
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
