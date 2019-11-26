package org.jax.bamextractor.processing;

import java.util.HashMap;

public class BaseToInt {

	private static HashMap<Character, int[]> _map;
	
	private static void initMap() {
		_map = new HashMap<Character, int[]>(30);
		_map.put('A', new int[]{0});
		_map.put('C', new int[]{1});
		_map.put('M', new int[]{0,1});
		_map.put('G', new int[]{2});
		_map.put('R', new int[]{0,2});
		_map.put('S', new int[]{1,2});
		_map.put('V', new int[]{0,1,2});
		_map.put('T', new int[]{3});
		_map.put('W', new int[]{0,3});
		_map.put('Y', new int[]{1,3});
		_map.put('H', new int[]{0,1,3});
		_map.put('K', new int[]{2,3});
		_map.put('D', new int[]{0,2,3});
		_map.put('B', new int[]{1,2,3});
		_map.put('N', new int[]{0,1,2,3});
	}
	
	public static int[] getIntBases(char base) {
		if(_map == null) {
			initMap();
		}
		return _map.get(base);
	}
	
}
