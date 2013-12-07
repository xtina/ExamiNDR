package ExamiNDR.v1;

import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class Gap {
	/*
	 * original matlab code function [gap, gapplot] = gapfinder(Location,
	 * Density, densLim, lengthLim)
	 * 
	 * M = find(Density<densLim); n=1; gapindex(n,1)=M(1);
	 * 
	 * for i=2:length(M) 
            if ~(M(i)-M(i-1)==1) 
                gapindex(n,2)=M(i-1); 
                n=n+1;
	 *      gapindex(n,1)=M(i); 
            end 
           end 
           gapindex(n,2)=M(length(M)); 
           gap1 = Location(gapindex);
	 * 
	 * m=1; gap(m,1)=gap1(1,1);
	 * 
	 * for i=2:size(gap1,1) if gap1(i,1)-gap1(i-1,2)>lengthLim
	 * gap(m,2)=gap1(i-1,2); m=m+1; gap(m,1)=gap1(i,1); end end
	 * gap(m,2)=gap1(size(gap1,1),2);
	 * 
	 * j=1; i=1; k=1; 
    while (j<=size(gap,1)&& i<=length(Location)-1) 
        if Location(i) >= gap(j,1)&& Location(i) <= gap(j,2) 
            gapplot(k,1)= Location(i); gapplot(k,2)= densLim; k=k+1; 
        end 
        if Location(i) <= gap(j,2)&& Location(i+1) > gap(j,2) j=j+1; 
        end 
        i=i+1; 
    end
	 */
	

	private ArrayList<Point2D.Double> densities, metLimit;// = new ArrayList<Point2D.Double>(); //contains locations and densities of a genome sequence
	//private ArrayList<Point2D.Double> metLimit// = new ArrayList<Point2D.Double>();
        private static final int LOC = 0; //location indicie in array
	private static final int DENS = 1; //density indicie in array
	
	//Reads in a file containing locations and densities of a genome and parses values into ArrayList NDR. This file must be tab delineated.
	public void readFile(File inputFile){
            try {
                densities = new ArrayList<Point2D.Double>();
                metLimit = new ArrayList<Point2D.Double>();
                String line; //line contents
                String[] temp = new String[2]; //parse line into string array

                //read in lines and parse
                BufferedReader br = new BufferedReader(new FileReader(inputFile));
                while((line = br.readLine()) != null) {
                    temp = line.split("\t");
                    densities.add(new Point2D.Double(Double.parseDouble(temp[LOC]), Double.parseDouble(temp[DENS])));
                }
                br.close();
            } 
            catch (FileNotFoundException e) {
                    System.out.println("File " + inputFile + " not found.");
                    e.printStackTrace();
            }
            catch(IOException g) {
                    System.out.println("IO Exception in writeFile");
                    g.printStackTrace();
            }
	}
	
        // Outputs the current list of NDRs to file, format: startPosition TAB endPosition
        public void writeOutput(File outFile) {
            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
                for(Point2D.Double ndr : metLimit) {
                    bw.write(ndr.getX() + "\t" + ndr.getY() + "\n");
                }
                bw.flush();
                bw.close();
            }
            catch (FileNotFoundException e) {
                System.out.println("File " + outFile + " not found.");
                e.printStackTrace();
            }
            catch(IOException e) {
                System.out.println("IO Exception in writeFile");
                e.printStackTrace();
            }
        }
        
        public void gapFinder(double densityLimit, double lengthLimit, double gapTolerance) {
            //int size = densities.size();
            ArrayList<Point2D.Double> lowOccupancy = new ArrayList<Point2D.Double>();
            
            // first, identify locations with occupancy under densityLimit
            // while doing this, combine across gaps
            for(Point2D.Double bp : densities) {
                if(bp.y < densityLimit) {
                    // check if it should be combined with the previous ndr
                    if(lowOccupancy.size() > 0 && bp.x - lowOccupancy.get(lowOccupancy.size() - 1).y <= gapTolerance) {
                        lowOccupancy.set(lowOccupancy.size() - 1, new Point2D.Double(lowOccupancy.get(lowOccupancy.size() - 1).x, bp.x));
                    }
                    else {
                        lowOccupancy.add(new Point2D.Double(bp.x, bp.x));
                    }
                }
            }
            // next, check lengths
            for(Point2D.Double ndr : lowOccupancy) {
                if(ndr.y - ndr.x >= lengthLimit) {
                    metLimit.add(ndr);
                }
            }
        }
        
	/*public void gapFinder(double densityLimit, double lengthLimit, double gapTolerance){
            int size = densities.size();
            int i = 0;
            double startPos = 0;
            while(i < size) {
                // continue incrementing i until no longer contained within a region of low nucleosome density
                // at this point, check length of segment
                if(densities.get(i).getY() > densityLimit) {
                    if(i >= 1) {
                        double endPos = densities.get(i-1).getX();
                        // check to see whether it was long enough
                        if((endPos - startPos) > lengthLimit) {
                            // check if it should be combined with the previous ndr
                            if(metLimit.size() > 0) {
                                Point2D.Double prev = metLimit.get(metLimit.size() - 1);
                                if(startPos - prev.getY() < gapTolerance) {
                                    double prevStart = prev.getX();
                                    metLimit.set(metLimit.size() - 1, new Point2D.Double(prevStart, endPos));
                                }
                                else {
                                    metLimit.add(new Point2D.Double(startPos, endPos));
                                }
                            }
                            else {
                                metLimit.add(new Point2D.Double(startPos, endPos));
                            }
                        }
                    }
                    if(i < densities.size() - 1) {
                        startPos = densities.get(i+1).getX();
                    }
                }
                i++;
            }
	}*/
	/* returns the number of currently identified NDRs
	*/
	public int getNDRCount() {
		return metLimit.size();
	}
}
