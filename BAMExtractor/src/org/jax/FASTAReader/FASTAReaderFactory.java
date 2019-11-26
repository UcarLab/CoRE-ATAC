package org.jax.FASTAReader;

import java.io.IOException;

public class FASTAReaderFactory {

	public static FASTAReader getDefaultFASTAReader(String src, int charperline) throws IOException{
		return new ForwardFASTAReader(src, charperline);
	}
	public static FASTAReader getForwardFASTAReader(String src, int charperline) throws IOException{
		return new ForwardFASTAReader(src, charperline);
	}
	
}
