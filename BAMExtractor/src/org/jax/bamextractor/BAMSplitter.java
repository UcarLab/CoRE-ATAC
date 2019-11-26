package org.jax.bamextractor;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.TreeMap;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class BAMSplitter {

	public static void main(String[] args) {
		if(args.length == 3) {
			BAMSplitter b = new BAMSplitter();
			String bamfile = args[0];
			String chromosomefile = args[1];
			String outdir = args[2];
			try {
				b.writeThreshold(outdir, b.getThresholdAndSplitBAM(bamfile, outdir, b.readChromosomes(chromosomefile)));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private void writeThreshold(String outdir, int thresh) throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outdir+"/threshold.txt"));
		bw.write(Integer.toString(thresh));
		bw.flush();
		bw.close();
	}
	
	private String[] readChromosomes(String file) throws IOException {
		LinkedList<String> rv = new LinkedList<String>();
		BufferedReader br  = new BufferedReader(new FileReader(file));
		while(br.ready()) {
			rv.add(br.readLine());
		}
		br.close();
		return rv.toArray(new String[0]);
	}
	
	private int getThresholdAndSplitBAM(String bamfile, String outdir, String[] chromosomes){
		SamReaderFactory factory = SamReaderFactory.makeDefault()
	              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
	              .validationStringency(ValidationStringency.SILENT);
		final SamReader reader = factory.open(new File(bamfile));
		
		TreeMap<String, SAMFileWriter> writers = new TreeMap<String, SAMFileWriter>();
		
		for(int i = 0; i < chromosomes.length; i++){
			String key = chromosomes[i];
			writers.put(key, new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), 
				true, new File(outdir+"/"+key+".bam")));
		}
		
		SAMRecordIterator it = reader.iterator();
		it = it.assertSorted(SAMFileHeader.SortOrder.coordinate);
		ReadStatistics stats = new ReadStatistics();
		
		while(it.hasNext()){
			SAMRecord next = it.next();
			stats.addRecord(next);
			SAMFileWriter w = writers.get(next.getReferenceName());
			if(w != null){
				w.addAlignment(next);
			}
		}
		
		while(!writers.isEmpty()){
			writers.pollFirstEntry().getValue().close();
		}
		
		int thresh = 20;
		int rv = stats.getInsertSizeWidthThreshold(thresh);
		
		stats.printStatistics(thresh);
		stats.trim(rv);
		it.close();
		
		return rv;
	}
	
}
