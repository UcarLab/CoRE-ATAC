package org.jax.bamextractor;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;
import java.util.TreeMap;

public class FileMap {

	public static Map<String, String> getFileMap(String file) throws IOException{
		Map<String, String> rv = new TreeMap<String, String>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		while(br.ready()) {
			String line = br.readLine();
			String[] split = line.split("\t");
			rv.put(split[0], split[1]);
		}
		br.close();
		return rv;
	}
	
}
