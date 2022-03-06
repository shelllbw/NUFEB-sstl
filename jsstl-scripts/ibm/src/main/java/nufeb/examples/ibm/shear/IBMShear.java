package nufeb.examples.ibm.shear;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;




import nufeb.examples.ibm.shear.IBMShearScript;
//import eu.quanticol.jsstl.core.dsl.ScriptLoader;
import eu.quanticol.jsstl.core.formula.Signal;
import eu.quanticol.jsstl.core.formula.jSSTLScript;
import eu.quanticol.jsstl.core.space.GraphModel;
import eu.quanticol.jsstl.core.io.TraGraphModelReader;

public class IBMShear {

	public static void main(String[] args) throws Exception {
		//// /// %%%%%%% GRAPH %%%%% /////////
		TraGraphModelReader graphread = new TraGraphModelReader();
		GraphModel graph = graphread.read("models/sstl_shear_model.txt");
		long start = System.currentTimeMillis();

		System.out.println("Compute shortest path...");
		//graph.dMcomputation();
		double dmtime = (System.currentTimeMillis()-start)/1000.0;
		System.out.println("DM Computation: "+dmtime);

		//// /// %%%%%%% PROPERTY %%%%% //////////
		jSSTLScript script = new IBMShearScript();
		// /// %%%%%%% DATA import %%%%%%%%%%%%/////////

		check(script, graph, getShearParams());

		System.out.println("Done");
		
		dmtime = (System.currentTimeMillis()-start)/1000.0;
		System.out.println("MC Computation: "+dmtime);

	}

	private static ArrayList<HashMap<String, Double>> getShearParams() {
		ArrayList<HashMap<String, Double>> toReturn = new ArrayList<>();
		double pre = 0;
		for (double i = 0; i < 511200; i = i + 1800) {
			// define sstl region
			HashMap<String,Double> map1 = new HashMap<>();
			map1.put("minT", pre);
			map1.put("maxT", i);
			map1.put("fromX", 0.8e-4);
			map1.put("toX", 2e-4);
			toReturn.add(map1);
			
			if(i != 0) {
				pre = i - 1800;
			}
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

		Signal s = readMeshTraj(graph, "models/trajectories/shear-0.15/volf_all.txt");
		Signal s1 = readMeshTraj(graph, "models/trajectories/shear-0.15/gridx.txt");
		
		for (int i = 0; i < params.size(); i++) {
			
			double max = 0.0;
			
			HashMap<String,Double> parValues = params.get(i);
			
			PrintWriter stprinter = new PrintWriter("results/shear-0.15/"+"t" +i+".txt");
			
			//Signal s = readMeshTraj(graph, "models/trajectories/test/test.txt");

			double[] boolSat1 = script.booleanSat("streamer-formation", parValues, graph, s);
			double[] boolSat3 = script.booleanSat("grid-bound", parValues, graph, s1);
			
			System.out.print("checking : " + parValues.get("minT") + "-x" + + parValues.get("fromX"));
			
			for (int j = 0; j < graph.getNumberOfLocations(); j++) {
				if(boolSat1[j] > 0&& boolSat3[j] > 0) {
					if (max < s1.getValue(0,j,0))
						max = s1.getValue(0,j,0);
					stprinter.println(1.0);
				}
				else
					stprinter.println(0);
			}
			
			stprinter.close();
			System.out.println("  max, "+max);
		}
	}
}
