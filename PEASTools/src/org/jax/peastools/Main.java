package org.jax.peastools;
import org.jax.peastools.annotate.AnnotatePeaks;
import org.jax.peastools.conservation.Conservation;
import org.jax.peastools.filter.BAMFilter;
import org.jax.peastools.filter.FilterRegions;
import org.jax.peastools.insertmetrics.PeakInsertMetrics;
import org.jax.peastools.insertmetrics.PeakThreshold;
import org.jax.peastools.merge.Merge;
import org.jax.peastools.merge.MergeMACSFree;


public class Main {

	public static void main(String[] args){
		if(args.length == 0){
			System.out.println("PEAS Tools Options:");
			System.out.println("insertsizethresh");
			System.out.println("insertmetrics");
			System.out.println("conservation");
			System.out.println("merge");
			System.out.println("mergeMACSFree");
			System.out.println("filter");
			System.out.println("bamfilter");
			System.out.println("annotate");
			return;
		}
		String protocol = args[0].toLowerCase();
		
		String[] restargs = new String[args.length-1];
		int j = 0;
		for(int i = 1; i < args.length; i++){
			restargs[j++] = args[i];
		}
		
		if(protocol.equals("insertsizethresh")){
			PeakThreshold.main(restargs);
		}
		if(protocol.equals("insertmetrics")){
			PeakInsertMetrics.main(restargs);
		}
		else if(protocol.equals("conservation")){
			Conservation.main(restargs);
		}
		else if(protocol.equals("merge")){
			Merge.main(restargs);
		}
		else if(protocol.equals("mergedl")){
			System.out.println("GOT HERE");
			System.out.println(restargs[0]);
			MergeMACSFree.main(restargs);
		}
		else if(protocol.equals("filter")){
			FilterRegions.main(restargs);
		}
		else if(protocol.equals("bamfilter")){
			BAMFilter.main(restargs);
		}
		else if(protocol.equals("annotate")){
			AnnotatePeaks.main(restargs);
		}
	}
	
}
