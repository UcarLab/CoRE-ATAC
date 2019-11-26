package org.jax.bamextractor.util;


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;


public class Util {

	public TreeMap<String, Location[]> getChrStartSorted(Location[] pe){
		TreeMap<String, Location[]> chrsorted = getChrSorted(pe);
		TreeMap<String, Location[]> rv = new TreeMap<String, Location[]>();
		while(!chrsorted.isEmpty()){
			Entry<String, Location[]> e = chrsorted.pollFirstEntry();
			rv.put(e.getKey(), getSortedLocations(e.getValue()));
		}	
		return rv;
	}
	
	public TreeMap<String, Location[]> getChrSorted(Location[] pe){
		TreeMap<String, LinkedList<Location>> tm = new TreeMap<String, LinkedList<Location>>();
		for(int i = 0; i < pe.length; i++){
			String chr = pe[i].getChr().toLowerCase();
			if(!tm.containsKey(chr)){
				tm.put(chr, new LinkedList<Location>());
			}
			tm.get(chr).add(pe[i]);
		}
		
		TreeMap<String, Location[]> rv = new TreeMap<String, Location[]>();
		while(!tm.isEmpty()){
			Entry<String, LinkedList<Location>> e = tm.pollFirstEntry();
			rv.put(e.getKey(), e.getValue().toArray(new Location[0]));
		}
		return rv;
	}
	
	private Location[] getSortedLocations(Location[] pe){
		TreeSet<Location> ts = new TreeSet<Location>();
		for(int i = 0; i < pe.length; i++) {
			ts.add(pe[i]);
		}
		return ts.toArray(new Location[0]);
	}

	
	public Location[] readLocations(String file) throws IOException{
		LinkedList<Location> rv = new LinkedList<Location>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		int id = 1;
		while(br.ready()){
			String line = br.readLine();
			String[] split = line.split("\t");
			if(split.length >= 3){
				try{
					String sid = Integer.toString(id);
					if(split.length > 3) {
						sid = split[3];
					}
					rv.add(new Location(sid, split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2])));
				}
				catch(NumberFormatException e){
					System.out.println("Number Format Exception: "+line);
				}
			}
			id++;
		}
		br.close();
		return rv.toArray(new Location[0]);
	}
	
	
/*	private int findFirstIndex(Location[] loci, Location l){
		int start = l.getStart();
		int end = l.getEnd();
		
		//Binary search to get an approximate index
		int s = 0;
		int e = loci.length-1;
		while((e-s) > 1){
			int mi = s+((e-s)/2);
			Location ml = loci[mi];
			int mstart = ml.getStart();
			if(mstart < start){
				s = mi;
			}
			else if(mstart > start){
				e = mi;
			}
			else{
				s = mi;
				e = mi;
				break;
			}
		}
		
		//Loop to make sure we have the location just before the start position goes beyond
		for(int i = s+1; i < loci.length; i++){
			if(loci[i].getStart() > start){
				s = i-1;
				break;
			}
		}
		
		//Choose position just before or just after depending on the overlap. 
		//If neither of them do return -1
		int n = s+1;
		if(s < loci.length && loci[s].getEnd() > start){ //0base
			return s;
		}
		else if(n < loci.length && loci[n].getStart() <= end && loci[n].getEnd() > start){ //0base
			return n;
		}
		else {
			return -1;
		}
	}
*/	

}
