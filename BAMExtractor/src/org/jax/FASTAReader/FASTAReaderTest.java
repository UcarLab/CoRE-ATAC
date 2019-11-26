package org.jax.FASTAReader;

import java.io.IOException;

public class FASTAReaderTest {

	
	public static void main(String[] args) {
		FASTAReader fr;
		try {
			fr = FASTAReaderFactory.getForwardFASTAReader("/Users/athib/Desktop/ATAC Prediction/PipelineScripts/homer/data/genomes/hg19/chr18.fa", 50);
			fr.setRegion(3497677, 34971477);
			System.out.println(fr.getFASTASequence(34970127, 34971127));
			fr.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
