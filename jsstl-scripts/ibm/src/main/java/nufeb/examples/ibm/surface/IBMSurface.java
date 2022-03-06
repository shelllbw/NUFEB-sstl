package nufeb.examples.ibm.surface;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;


//import eu.quanticol.jsstl.core.dsl.ScriptLoader;
import eu.quanticol.jsstl.core.formula.Signal;
import eu.quanticol.jsstl.core.formula.jSSTLScript;
import eu.quanticol.jsstl.core.space.GraphModel;
import eu.quanticol.jsstl.core.io.TraGraphModelReader;

public class IBMSurface {

	public static void main(String[] args) throws Exception {
		//// /// %%%%%%% GRAPH %%%%% /////////
		TraGraphModelReader graphread = new TraGraphModelReader();
		GraphModel graph = graphread.read("models/sstl_surface_model.txt");
		long start = System.currentTimeMillis();

		System.out.println("Compute shortest path...");
		graph.dMNUFEBcomputation(30, 12, 24);
		double dmtime = (System.currentTimeMillis()-start)/1000.0;
		System.out.println("DM Computation: "+dmtime);

		//// /// %%%%%%% PROPERTY %%%%% //////////
		jSSTLScript script = new IBMSurfaceScript();
		// /// %%%%%%% DATA import %%%%%%%%%%%%/////////

		check(script, graph, getStratiParams());

		System.out.println("Done");
		
		dmtime = (System.currentTimeMillis()-start)/1000.0;
		System.out.println("MC Computation: "+dmtime);

	}

	private static ArrayList<HashMap<String, Double>> getStratiParams() {
		ArrayList<HashMap<String, Double>> toReturn = new ArrayList<>();
		for (double i = 10000; i < 760000; i = i + 10000) {
			// surface
			HashMap<String,Double> map1 = new HashMap<>();
			map1.put("minT", i-1);
			map1.put("maxT", i);
			map1.put("minD", 0.0);
			map1.put("maxD", 1.0);
			toReturn.add(map1);
			
		}

		return toReturn;
	}

	
	public static Signal readMeshTraj(GraphModel graph, String filename) throws IOException{
		FileReader fr = new FileReader(filename);
		BufferedReader br = new BufferedReader(fr);
		ArrayList<String> strings = new ArrayList<String>();
		String line = null;
		int timeSize = 0;
		while ((line = br.readLine()) != null) {
			strings.add(line);
			String[] key = line.split(",");

			if (key[0].equals("Time"))
				timeSize ++;
		}
		System.out.println("timeSize: "+timeSize);
		int locSize = graph.getNumberOfLocations();
		System.out.println("locSize: "+locSize);
		double[][][] data = new double[locSize][timeSize][1];
		double[] time = new double[timeSize];
		
		for (int t = 0; t < timeSize; t++) {
			String[] splitted = strings.get(t*locSize+t).split(",");
			time[t] = Double.valueOf(splitted[1]);
		}

		for (int t = 0; t < timeSize; t++) {
			for (int i = 0; i < locSize; i++) {	
				String[] splitted = strings.get(t*locSize+i+t+1).split(",");
				int loc = Integer.valueOf(splitted[0])-1;
				data[loc][t][0] = Double.valueOf(splitted[1]);
				//System.out.println("loc: "+loc + " value: "+  data[loc][t][0]);
			}
		}

		br.close();
		fr.close();
		Signal s = new Signal(graph, time, data);
		return s;
	}

	public static void check(jSSTLScript script, GraphModel graph, ArrayList<HashMap<String,Double>> params) throws IOException {

		Signal s = readMeshTraj(graph, "models/trajectories/surface/volf_all.txt");

		for (int i = 0; i < params.size(); i++) {
			int total = 0;
			HashMap<String,Double> parValues = params.get(i);
			
			PrintWriter stprinter = new PrintWriter("data/surface" + i +".txt");
			
			//Signal s = readMeshTraj(graph, "models/trajectories/test/test.txt");

			double[] boolSat1 = script.booleanSat("surface-find", parValues, graph, s);
			
			System.out.print("checking : " + parValues.get("maxT") + "-x" + + parValues.get("maxD"));
			
			for (int j = 0; j < graph.getNumberOfLocations(); j++) {
				if(boolSat1[j] > 0) {
					stprinter.println(1.0);
					total ++;
				} else {
					stprinter.println(0.0);
				}
			}
			
			stprinter.close();
			System.out.println("  total : " + total);
		}
	}
	
}
