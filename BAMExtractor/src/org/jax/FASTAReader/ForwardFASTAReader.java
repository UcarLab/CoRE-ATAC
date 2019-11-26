package org.jax.FASTAReader;


import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;


public class ForwardFASTAReader implements FASTAReader{

	private BufferedReader _br;
	private int _charperline = 50;
	private int _start = 0;
	private int _end = 0;
	private String _fasta = "";
	
	public ForwardFASTAReader(String file, int charperline) throws IOException{
		if(file.endsWith(".gz")){
			_br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		}
		else{
			_br = new BufferedReader(new FileReader(file));
		}
		_br.readLine();
		_charperline = charperline;
	}

	@Override
	public void setRegion(int start, int end) throws IOException {
		start = Math.max(_start, start);
		end = Math.max(_end, end);
		
		if(_start < start){
			int diff = start-_end;
			if(diff > 0){
				_start = _end;
				_fasta = "";

				int skip = (diff/_charperline)-1;
	
				for(int i = 0; i < skip && _br.ready(); i++){
					_br.readLine();
					_start += _charperline;
				}
				_end = _start;
			}
		}
		
		StringBuilder sbuilder = new StringBuilder(_fasta);
		while(_end < end && _br.ready()){
			String line = _br.readLine();
			_end += line.length();
			sbuilder.append(line);
		}
		_fasta = sbuilder.toString();
		
		if (start > _start){
			_fasta = _fasta.substring(start-_start-1);
			_start = start-1;
		}
	}

	@Override
	public String getFASTASequence(int start, int end) {
		int s = Math.max(_start, start)-_start-1;
		int e = (Math.min(_end, end)-_start);
		return _fasta.substring(s,e);
	}

	@Override
	public void close() {
		try {
			_br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
}
