package org.jax.bamextractor.util;

import java.util.Comparator;

public class Location implements Comparable<Location>, Comparator<Location>{
	private String _id;
	protected String _chr;
	protected int _start;
	protected int _end;
	
	public Location(String id, String chr, int start, int end){
		_id = id;
		_chr = chr;
		_start = start;
		_end = end;
	}
	
	public String getId(){
		return _id;
	}
	
	public String getChr(){
		return _chr;
	}
	
	public int getStart(){
		return _start;
	}
	
	public int getEnd(){
		return _end;
	}
	
	public void setChr(String chr){
		_chr = chr;
	}
	
	public void setStart(int s){
		_start = s;
	}
	
	public void setEnd(int e){
		_end = e;
	}

	@Override
	public int compare(Location o1, Location o2) {
		int chrcomp = o1._chr.compareTo(o2._chr);
		if(chrcomp == 0) {
			if(o1._start < o2._start) {
				return -1;
			}
			else if(o1._start > o2._start) {
				return 1;
			}
			else {
				if(o1._end < o2._end) {
					return -1;
				}
				else if(o1._end > o2._end) {
					return 1;
				}
				else {
					return 0;
				}
			}
		}
		return chrcomp;
	}

	@Override
	public int compareTo(Location o) {
		return compare(this, o);
	}
}
