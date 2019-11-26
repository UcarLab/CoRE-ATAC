package org.jax.FASTAReader;

import java.io.IOException;

public interface FASTAReader {

	/**
	 * Sets the region of interest to extract the FASTA information.
	 * @param start The start of the region to look at.
	 * @param end The end of the region to look at.
	 */
	public void setRegion(int start, int end) throws IOException;
	
	/**
	 * 
	 * @param start The start position of the returned FASTA string.
	 * @param end The end position of the returned FASTA string.
	 * @return
	 */
	public String getFASTASequence(int start, int end);
	
	public void close();
	
}
