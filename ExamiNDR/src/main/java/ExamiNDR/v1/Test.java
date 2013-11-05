package ExamiNDR.v1;

import java.io.File;

public class Test {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		Gap gap = new Gap();
                int total = 0;
                for(int i = 1; i <= 16; i++) {
                    gap.readFile(new File("./TestFiles/" + i +"_occupancy.txt"));
                    gap.gapFinder(-1, 50, 100);
                    gap.writeOutput(new File("./Output/" + i + "_out.txt"));
                    System.out.println("identified " + gap.getNDRCount() + " NDRs in chromosome " + i);
                    total += gap.getNDRCount();
                }
		System.out.println("identified " + total + " total NDRs");
	}

}
