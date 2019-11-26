package org.jax.bamextractor;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordCoordinateComparator;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

public class InsertSizeThreshold {

	
	public int getInsertSizeThreshold(String chrbamfile) throws IOException{
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
		ReadStatistics stats = new ReadStatistics();
		
		while(it.hasNext()){
			SAMRecord next = it.next();
			stats.addRecord(next);
		}
		
		int thresh = 20;
		int rv = stats.getInsertSizeWidthThreshold(thresh);
		
		it.close();
		reader.close();
		
		return rv;
	}
	
}
