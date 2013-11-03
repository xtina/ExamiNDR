package ExamiNDR.v1;

import java.awt.Point;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

public class Gap {
	/*
	 * original matlab code function [gap, gapplot] = gapfinder(Location,
	 * Density, densLim, lengthLim)
	 * 
	 * M = find(Density<densLim); n=1; gapindex(n,1)=M(1);
	 * 
	 * for i=2:length(M) if ~(M(i)-M(i-1)==1) gapindex(n,2)=M(i-1); n=n+1;
	 * gapindex(n,1)=M(i); end end gapindex(n,2)=M(length(M)); gap1 =
	 * Location(gapindex);
	 * 
	 * m=1; gap(m,1)=gap1(1,1);
	 * 
	 * for i=2:size(gap1,1) if gap1(i,1)-gap1(i-1,2)>lengthLim
	 * gap(m,2)=gap1(i-1,2); m=m+1; gap(m,1)=gap1(i,1); end end
	 * gap(m,2)=gap1(size(gap1,1),2);
	 * 
	 * j=1; i=1; k=1; while (j<=size(gap,1)&& i<=length(Location)-1) if
	 * Location(i) >= gap(j,1)&& Location(i) <= gap(j,2) gapplot(k,1)=
	 * Location(i); gapplot(k,2)= densLim; k=k+1; end if Location(i) <=
	 * gap(j,2)&& Location(i+1) > gap(j,2) j=j+1; end i=i+1; end
	 */
	
	private ArrayList<Point2D.Double> NDR = new ArrayList<Point2D.Double>(); //contains locations and densities of a genome sequence
	private static final int LOC = 0; //location indicie in array
	private static final int DENS = 1; //density indicie in array
	
	
	//Reads in a file containing locations and densities of a genome. This file must be tab delineated.
	public void readFile(File inputFile){
		try {
			String line; //line contents
			String[] temp = new String [2]; //parse line into string array
			BufferedReader br = new BufferedReader(new FileReader(inputFile));
			while((line = br.readLine()) != null) {
				temp = line.split("\t");
				NDR.add(new Point2D.Double(Double.parseDouble(temp[LOC]), Double.parseDouble(temp[DENS])));
			}
			
			for(Point2D.Double p : NDR) {
				System.out.print(p.getX() + " ");
				System.out.println(p.getY());
			}
		} 
		catch (FileNotFoundException e) {
			System.out.println("File " + inputFile + " not found.");
			e.printStackTrace();
		}
		catch(IOException g) {
			System.out.println("IO Exception in readFIle");
			g.printStackTrace();
		}
	}
	
	//Finds gaps
	public void gapFinder(int densityLimit, int lengthLimit){
		
	}
}