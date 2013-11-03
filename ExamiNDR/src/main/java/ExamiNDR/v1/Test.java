package ExamiNDR.v1;

import java.io.File;

public class Test {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		Gap gap = new Gap();
		gap.readFile(new File("./TestFiles/10_occupancy.txt"));
		gap.gapFinder(1, .2);
	}

}
