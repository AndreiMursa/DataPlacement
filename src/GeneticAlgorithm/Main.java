package GeneticAlgorithm;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import buildTimeStageAlgorithm.BuildTimeObject;

public class Main {
	// public static final double distributionWeight = 0.7 //prioritize data
	// distribution minimization
	// public static final double performanceWeight = 0.3 //prioritize performance
//	Genes to mutate 6
//	Time spent: 73 seconds
//	Time spent: 1 minutes
//	Best Fitness 1.9785458817596878E-4
//	Execution time 3661.2168833133937
//	Data scheduling times 1393
	private static List<Chromosome> individuals = new ArrayList<Chromosome>();
	private static List<Dataset> datasets = new ArrayList<Dataset>();
	private static List<DataCenter> dataCenters = new ArrayList<DataCenter>();

	public static void main(String args[]) {
		// generateDatasets();
		Utils.latencyMatrix = generateLatencyMatrix();
		generateDataCenters();
		generateDatasetsObjects();
		// displayMatrix(Utils.latencyMatrix);
		generatePopulation();			
		//writeDatasetsPlacementMatrix(individuals.get(0).getPlacementMatrix());
		Request request = new Request(Utils.numberOfComputations);
		request.generateComputations();
		request.createAssociationMatrix();
		int[][] associationMatrix = request.getAssociationMatrix();
		
		
//		BuildTimeObject btAlg = new BuildTimeObject(datasets, dataCenters, request.getComputations());
//		Chromosome cm = new Chromosome(Utils.datasetsCount, Utils.dataCentersCount);
//		cm.setPlacementMatrix(btAlg.getPlacementMatrix());
//		evaluate(cm, request, associationMatrix);
//		System.out.println(cm.getFitnessValue());
//		System.out.println("Data scheduling " + cm.getTotalNrDataScheduling());
//		System.out.println(cm.getTotalEPT());
//		System.out.println();
//		//displayMatrix(btAlg.getDependencyMatrix());
//		System.out.println("TERMINAT ALGORITHM BUILD TIME");
		

//	//	System.out.println("Placement Matrix");
//		int[][] placementMatrix = readDatasetsPlacementMatrix();
//	//	displayMatrix(placementMatrix);
////		for(Dataset d : datasets) {
////			System.out.println(d.getFileSize());
////		}

		long timeStart = System.nanoTime();
		int generation = 0;
		Chromosome bestSolution = null;
		double bestFitness = 0.0;
		int bestNrScheduling = 9999999;
		double bestExecTime = Double.MAX_VALUE;
		Chromosome best1 = null;
		Chromosome best2 = null;
		while (generation < Utils.maxGeneration) {
			System.out.println("Evaluation generation " + generation);
			evaluatePopulation(request, associationMatrix);
			double average = 0.0;
			for(Chromosome c: individuals) {
				average+=c.getFitnessValue();
			}
			System.out.println("Generation average: " + average/Utils.population);
			individuals = rouletteWheelSelection(bestSolution);
//			for (Chromosome c : individuals) {
//				System.out.println(c.getFitnessValue());
//			}

			for (Chromosome c : individuals) {
				if (c.getFitnessValue() > bestFitness) {
					bestFitness = c.getFitnessValue();
					bestSolution = c;
				}
				if (c.getTotalEPT() < bestExecTime) {
					best1 = c;
					bestExecTime = c.getTotalEPT();
				}
				if (c.getTotalNrDataScheduling() < bestNrScheduling) {
					best2 = c;
					bestNrScheduling = c.getTotalNrDataScheduling();
				}
			}
			if (generation % 10 == 0) {
				System.out.println("Best fitness at gen " + generation + " : " + bestFitness);
				System.out.println(
						"Best nr scheduling at gen " + generation + " : " + bestSolution.getTotalNrDataScheduling());
				System.out.println("Best time at gen " + generation + " : " + bestSolution.getTotalEPT());
			}
			crossoverOperation();
			mutationOperation();
			generation++;

		}
		long timeEnd = System.nanoTime();
		long timeElapsed = timeEnd - timeStart;
		System.out.println("Time spent: " + timeElapsed / 1000000000 + " seconds");
		System.out.println("Time spent: " + timeElapsed / 60000000000L + " minutes");
		System.out.println("Best Fitness " + bestFitness);
		System.out.println("Execution time " + bestSolution.getTotalEPT());
		System.out.println("Data scheduling times " + bestSolution.getTotalNrDataScheduling());
		System.out.println("BEST 1 Execution time chromosome: " + best1.getTotalEPT());
		System.out.println("BEST 1 Sched time chromosome: " + best1.getTotalNrDataScheduling());
		System.out.println();
		System.out.println("BEST 2 Execution time chromosome: " + best2.getTotalEPT());
		System.out.println("BEST 2 Sched time chromosome: " + best2.getTotalNrDataScheduling());

		// System.out.println("Best Result " + bestFiteness);
		// assignDatasetsToDatacenters(placementMatrix);
//
//		//System.out.println("Association Matrix");
//		Request request = new Request(Utils.numberOfComputations);
//		request.generateComputations();
//		request.createAssociationMatrix();
//		int[][] associationMatrix = request.getAssociationMatrix();
//		//displayMatrix(associationMatrix);
//		//System.out.println("Z Matrix");
//		int[][] matrixZ = matrixMultiplication(associationMatrix, placementMatrix);
//		//displayMatrix(matrixZ);
//		// System.out.println();
//		int[][] uMatrixZ = applyFunctionOnMatrix(matrixZ);
//		//displayMatrix(uMatrixZ);
//		List<Computation> computations = request.getComputations();
//		double totalEPT = 0.0;
//		for (Computation c : computations) {
//			double ept = getTotalEPTForComputation(c, placementMatrix, uMatrixZ);
//			totalEPT+= ept;
//		    //System.out.println("EPT for computation " + c.getIndex() + " is: "+ept + " seconds ??????, hope so");
//		}
//		System.out.println("Total execution time for all computations = " + totalEPT);
//		
//		int totalNrDataScheduling = 0;
//		for(int i=0; i< uMatrixZ.length;i++) {
//			int nrDataScheduling = 0;
//			for(int j=0; j< uMatrixZ[1].length;j++) {
//				nrDataScheduling+= uMatrixZ[i][j];
//			}
//			if(nrDataScheduling >= 1) {
//				nrDataScheduling--;
//				System.out.println("For computation " + i + " there are " + nrDataScheduling + " data scheduling");
//				totalNrDataScheduling+= nrDataScheduling;
//			}
//		}
//		
//		System.out.println("Total nr data scheduling " + totalNrDataScheduling);
//		System.out.println("Placement matrix");
//		displayMatrix(placementMatrix);
//		System.out.println("Association matrix");
//		displayMatrix(associationMatrix);
//		System.out.println("uMatrixZ matrix");
//		displayMatrix(uMatrixZ);
	}

	public static void evaluatePopulation(Request request, int[][] associationMatrix) {
		double bestFiteness = 0.0;
		individuals.stream().parallel().forEach(c -> evaluate(c, request, associationMatrix));
//		for (Chromosome chromosome : individuals) {
//			evaluate(chromosome, request, associationMatrix);
//
//		}
	}

	public static void evaluate(Chromosome chromosome, Request request, int[][] associationMatrix) {
		assignDatasetsToDatacenters(chromosome.getPlacementMatrix());
		int[][] matrixZ = matrixMultiplication(associationMatrix, chromosome.getPlacementMatrix());
		int[][] uMatrixZ = applyFunctionOnMatrix(matrixZ);
		List<Computation> computations = request.getComputations();
		double totalEPT = 0.0; // CHANGE BACK TO ZERO
		for (Computation c : computations) {
			double ept = getTotalEPTForComputation(c, chromosome.getPlacementMatrix(), uMatrixZ);
			totalEPT += ept;
			// System.out.println("EPT for computation " + c.getIndex() + " is: "+ept + "
			// seconds ??????, hope so");
		}
		// System.out.println("Total execution time for all computations = " +
		// totalEPT);

		int totalNrDataScheduling = 0;
		for (int i = 0; i < uMatrixZ.length; i++) {
			int nrDataScheduling = 0;
			for (int j = 0; j < uMatrixZ[1].length; j++) {
				nrDataScheduling += uMatrixZ[i][j];
			}
			if (nrDataScheduling >= 1) {
				nrDataScheduling--;
				// System.out.println("For computation " + i + " there are " + nrDataScheduling
				// + " data scheduling");
				totalNrDataScheduling += nrDataScheduling;
			}
		}

		// totalEPT =3*totalEPT/100;
		chromosome.setTotalNrDataScheduling(totalNrDataScheduling);
		chromosome.setTotalEPT(totalEPT);
		chromosome.setFitnessValue(1 / (totalEPT + totalNrDataScheduling));

//		System.out.println("Total nr data scheduling " + totalNrDataScheduling);
//		System.out.println("Total nr ept " + totalEPT);
//		System.out.println(1/(totalEPT + totalNrDataScheduling));
//		if(1/(totalEPT + totalNrDataScheduling) >= bestFiteness) {
//			bestFiteness = 1/(totalEPT + totalNrDataScheduling);
//		}
	}

	public static List<Chromosome> rouletteWheelSelection(Chromosome bestSolution) {
		List<Chromosome> newPopulation = new ArrayList<Chromosome>();
		List<Double> probabilities = new ArrayList<Double>();
		double totalFitnessSum = 0.0;
		for (Chromosome c : individuals) {
			totalFitnessSum += c.getFitnessValue();
		}
		double nextCheckpoint = 0.0;
		for (Chromosome c : individuals) {
			double prob = c.getFitnessValue() / totalFitnessSum;
			nextCheckpoint += prob;
			probabilities.add(nextCheckpoint);
		}
		int k = 0;
		if (bestSolution != null) {
			for (int i = 0; i < Utils.population/10; i++) {
				newPopulation.add(bestSolution);
				k++;
			}
		}
		Random randNum = new Random();
		for (int i = k; i < Utils.population; i++) {
			double selectionProbability = randNum.nextDouble();
			for (int j = 0; j < probabilities.size(); j++) {
				if (selectionProbability <= probabilities.get(j)) {
					newPopulation.add(individuals.get(j));
					break;
				}
			}
		}

		return newPopulation;

	}

	private static void crossoverOperation() {
		Random randNum = new Random();
		List<Chromosome> newPopulation = new ArrayList<Chromosome>();
		for (int i = 0; i < individuals.size(); i += 2) {
			
			Chromosome child1 = new Chromosome(Utils.datasetsCount, Utils.dataCentersCount);
			Chromosome child2 = new Chromosome(Utils.datasetsCount, Utils.dataCentersCount);
			Chromosome parent1 = individuals.get(i);
			Chromosome parent2 = individuals.get(i + 1);
			int[][] matrix1 = parent1.getPlacementMatrix();
			int[][] matrix2 = parent2.getPlacementMatrix();
			int crossoverPoint =(int) (matrix1.length * Utils.crossoverRate); //1 + randNum.nextInt(Utils.datasetsCount - 3);
			for (int k = crossoverPoint; k < matrix1.length; k++) {
				int[] temp = matrix1[k];
				matrix1[k] = matrix2[k];
				matrix2[k] = temp;
			}

			child1.setPlacementMatrix(matrix1);
			child2.setPlacementMatrix(matrix2);
			newPopulation.add(child1);
			newPopulation.add(child2);
		}

		individuals = newPopulation;
	}

	private static void mutationOperation() {
		int numberOfGenesToMutate = (int) (Utils.datasetsCount * 70 / 100);
		// System.out.println("Genes to mutate " + numberOfGenesToMutate);
		individuals.stream().parallel().forEach(c -> mutate(c, numberOfGenesToMutate));
//		for(Chromosome chromosome : individuals) {
//			
//		}
	}

	private static void mutate(Chromosome chromosome, int numberOfGenesToMutate) {
		Random randNum = new Random();
		double mutationChance = randNum.nextDouble();
		if (mutationChance <= Utils.mutationRate) {
			int nrOfGenesMutated = 0;
			List<Integer> mutatedGenes = new ArrayList<Integer>();
			while (nrOfGenesMutated < numberOfGenesToMutate) {
				int gene = randNum.nextInt(Utils.datasetsCount);
				if (!mutatedGenes.contains(gene)) {
					mutatedGenes.add(gene);

					int[] array = chromosome.getPlacementMatrix()[gene];
					// shuffle the datasets in datacenters
					for (int i = 0; i < array.length; i++) {
						int randomIndexToSwap = randNum.nextInt(array.length);
						int temp = array[randomIndexToSwap];
						array[randomIndexToSwap] = array[i];
						array[i] = temp;
					}

					nrOfGenesMutated++;
				}
			}
		}
	}

	private static void generatePopulation() {
		int discardedIndividuals = 0;
		for (int i = 0; i < Utils.population; i++) {
			Chromosome newIndividual = new Chromosome(Utils.datasetsCount, Utils.dataCentersCount);
			newIndividual.generateInitialPlacement(Utils.datasetsCount, Utils.dataCentersCount,
					Utils.replicationFactor);
			// newIndividual.displayIndividual();
			if (verifyIndividual(newIndividual.getPlacementMatrix())) {
				individuals.add(newIndividual);
			} else {
				discardedIndividuals++;
				i--;
			}

		}

		System.out.println("Weak bloodlines = " + discardedIndividuals);
	}

	public static boolean verifyIndividual(int[][] placementMatrix) {
		for (int j = 0; j < placementMatrix[1].length; j++) {
			long sumDatasetsSize = 0;
			for (int i = 0; i < placementMatrix.length; i++) {
				if (placementMatrix[i][j] == 1) {
					sumDatasetsSize += getDatasetAtIndex(i).getFileSize();
				}
			}
			if (sumDatasetsSize >= getDataCenterAtIndex(j).getCapacity()) {
				return false;
			}
		}

		return true;
	}

	private static void generateDataCenters() {
		Random randNum = new Random();
		for (int i = 0; i < Utils.dataCentersCount; i++) {
			double mbps = 8.0 + randNum.nextDouble() * 10;
			System.out.println("Mbps = " + mbps);
			double GHz = 1.0 + randNum.nextDouble() * 2;
			System.out.println("GHz = " + GHz);
			dataCenters.add(new DataCenter(i, Utils.dataCenterCapacity, mbps, GHz, Utils.latencyMatrix[i]));
		}
	}

	private static void generateDatasetsObjects() {
		for (int i = 0; i < Utils.datasetsCount; i++) {
			File fout = new File("out" + i + ".txt");
			datasets.add(new Dataset(fout, i));
		}
	}

	private static void assignDatasetsToDatacenters(int[][] placementMatrix) {
		for (int j = 0; j < Utils.dataCentersCount; j++) {
			List<Dataset> datasetsForDatacenter = new ArrayList<Dataset>();
			for (int i = 0; i < Utils.datasetsCount; i++) {
				if (placementMatrix[i][j] == 1) {
					datasetsForDatacenter.add(getDatasetAtIndex(i));
				}
			}
			getDataCenterAtIndex(j).setDatasetsToProcess(datasetsForDatacenter);
		}
	}

	private static void generateDatasets() {
		int createdDatasets = 0;
		BufferedReader br = null;
		try {
			// C:\\Users\\Andrei\\Desktop\\Files\\Test-Data.txt
			br = new BufferedReader(new FileReader("full_log1.csv"));
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		String line = "";
		try {

			Random randNum = new Random();
			for(int i=0;i < 10;i++) {
				long readLines = 0;
				br = new BufferedReader(new FileReader("full_log1.csv"));
				File fout = new File("out" + createdDatasets + ".txt");
				FileOutputStream fos = new FileOutputStream(fout);
				BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fos));
				int numberOfLines = (1+randNum.nextInt(5))*Utils.linesPerFile;
				while ((line = br.readLine()) != null && readLines <= numberOfLines) {
					if (readLines == numberOfLines) {
						readLines = 0;
						numberOfLines = (1+randNum.nextInt(5))*Utils.linesPerFile;
						createdDatasets++;
						bw.close();
						fout = new File("out" + createdDatasets + ".txt");
						fos = new FileOutputStream(fout);
						bw = new BufferedWriter(new OutputStreamWriter(fos));
					}
					bw.write(line);
					bw.newLine();
					readLines++;
				}
				bw.close();
				br.close();
			}


		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static double calculateFitnessValue(Chromosome individual, int[][] associationMatrix) {
		double value = 0.0;
		int[][] placementMatrix = individual.getPlacementMatrix();
		int[][] matrixZ = matrixMultiplication(associationMatrix, placementMatrix);

		return value;
	}

	private static void writeDatasetsPlacementMatrix(int[][] bestSolution) {
		BufferedWriter bw = null;
		try {
			File fout = new File("initial_placement.txt");
			FileOutputStream fos = new FileOutputStream(fout);
			bw = new BufferedWriter(new OutputStreamWriter(fos));

			for (int i = 0; i < Utils.datasetsCount; i++) {
				for (int j = 0; j < Utils.dataCentersCount; j++) {
					bw.write(bestSolution[i][j] + " ");
				}
				bw.newLine();
			}

			bw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static int[][] readDatasetsPlacementMatrix() {
		int[][] initialPlacement = new int[Utils.datasetsCount][Utils.dataCentersCount];
		BufferedReader br = null;
		try {
			// C:\\Users\\Andrei\\Desktop\\Files\\Test-Data.txt
			br = new BufferedReader(new FileReader("initial_placement.txt"));

			String line = "";
			int j = 0;
			while ((line = br.readLine()) != null) {
				String[] values = line.split(" ");
				for (int i = 0; i < Utils.dataCentersCount; i++) {
					initialPlacement[j][i] = Integer.parseInt(values[i]);
				}
				j++;
			}
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return initialPlacement;
	}

	public static void displayMatrix(int[][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[1].length; j++) {
				System.out.print(matrix[i][j] + " ");
			}
			System.out.println();
		}
	}

	public static int[][] matrixMultiplication(int[][] A, int[][] B) {
		int[][] c = new int[A.length][B[1].length];

		for (int i = 0; i < A.length; i++) {
			for (int j = 0; j < B[1].length; j++) {
				for (int k = 0; k < A[1].length; k++) {
					c[i][j] += A[i][k] * B[k][j];
				}
			}
		}

		return c;
	}

	public static int[][] applyFunctionOnMatrix(int[][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[1].length; j++) {
				if (matrix[i][j] != 0) {
					matrix[i][j] = 1;
				}
			}
		}

		return matrix;
	}

	public static double getTotalEPTForComputation(Computation c, int[][] placementMatrix, int[][] uMatrixZ) {
		double totalEPT = 0.0;
		int index = c.getIndex();
		double[][] edtsMatrix = new double[Utils.dataCentersCount][Utils.dataCentersCount];
		double taskSize = 0.0;
		for (Integer idx : c.getProcessedDatasets()) {
			taskSize += getDatasetAtIndex(idx).getFileSize();
		}
		for (int i = 0; i < Utils.latencyMatrix.length; i++) {
			for (int j = 0; j < Utils.latencyMatrix[1].length; j++) {
				DataCenter currentDataCenter = getDataCenterAtIndex(i);
				if (i == j) {
					edtsMatrix[i][j] = 0;
				} else {
					edtsMatrix[i][j] = taskSize / currentDataCenter.getEstimatedBandwidth()
							+ (double) Utils.latencyMatrix[i][j] / 1000;
				}
			}
		}
		List<Integer> usedDatacenters = new ArrayList<Integer>();

		for (int i = 0; i < uMatrixZ[1].length; i++) {
			if (uMatrixZ[index][i] != 0) {
				usedDatacenters.add(i);
			}
		}

		List<Integer> visitedDataCenters = new ArrayList<Integer>();
		int lastVisited = -1;
		for (int i = 0; i < uMatrixZ[1].length; i++) {
			if (uMatrixZ[index][i] != 0) {
				// usedDatacenters.add(i);
				visitedDataCenters.add(i);
				// first should always be the first found datacenter, and the edt is 0( no
				// scheduling , this is where is starts)
				if (visitedDataCenters.size() < usedDatacenters.size()) {
					boolean ok = true;
					while (ok) {
						lastVisited = getMinFromMAtrix(edtsMatrix, "row", i);
						if (!visitedDataCenters.contains(lastVisited)) {
							visitedDataCenters.add(lastVisited);
							totalEPT += edtsMatrix[i][lastVisited];
							ok = false;
						} else {
							// if datacenter already visited, add max value so we dont choose it again
							edtsMatrix[i][lastVisited] = Double.MAX_VALUE;
						}
					}

				}
				break;
			}
		}

		while (!(visitedDataCenters.size() == usedDatacenters.size())) {
			int idx = getMinFromMAtrix(edtsMatrix, "row", lastVisited);
			if (!visitedDataCenters.contains(idx)) {
				// edts.add(edtsMatrix[lastVisited][idx]);
				totalEPT += edtsMatrix[lastVisited][idx];
				lastVisited = idx;
				visitedDataCenters.add(idx);
			} else {
				// if datacenter already visited, add max value so we dont choose it again
				edtsMatrix[lastVisited][idx] = Double.MAX_VALUE;
			}
		}

		for (Integer idx : usedDatacenters) {
			taskSize = 0.0;
			DataCenter dataCenter = getDataCenterAtIndex(idx);
			List<Integer> computationDatasets = c.getProcessedDatasets();
			for (Dataset d : dataCenter.getDatasetsToProcess()) {
				if (computationDatasets.contains(d.getIndex())) {
					taskSize += d.getFileSize();
				}
			}

			totalEPT += taskSize / dataCenter.getEstimatedProcessingCapability();

		}

//		for (int j = 0; j < uMatrixZ[1].length; j++) {
//			for (int k = 0; k < uMatrixZ[1].length; k++) {
//				if (j != k && uMatrixZ[index][j] != 0 && uMatrixZ[index][k] != 0) {
//					double edt = 0.0;
//					DataCenter currentDataCenter = getDataCenterAtIndex(k);
//					double bdw = currentDataCenter.getEstimatedBandwidth();
//					double latency = Utils.latencyMatrix[j][k];
//					double taskSize = 0.0;
//					for (Integer idx : c.getProcessedDatasets()) {
//						taskSize += getDatasetAtIndex(idx).getFileSize();
//					}
//					taskSize = taskSize * 10 / 100;
//					// System.out.println("For computation " + index + " task size = " + taskSize);
//					edt = taskSize / bdw + latency;
//					edts.put("" + j + "->" + k, edt);
//				}
//			}
//		}
//
//		double min = Double.MAX_VALUE;
//		String executionOrder = "";
//		for (String edt : edts.keySet()) {
//			if (edts.get(edt) < min) {
//				min = edts.get(edt);
//				executionOrder = edt;
//			}
//		}
//
//		List<String> previouslyVisitedDatacenters = new ArrayList<String>();
//		List<Double> edtsInOrder = new ArrayList<Double>();
//		edtsInOrder.add(0.0); // first edt is 0, because here we start the computation
//		int counter = 1;
//		if (!executionOrder.isEmpty()) {
//			String currentDataset = executionOrder.split("->")[1];
//			previouslyVisitedDatacenters.add(executionOrder.split("->")[0]);
//			// Find the best order to get multiply the EDT with EPT for every data
//			// scheduling
//			boolean allDatacentersVisited = false;
//			while (!allDatacentersVisited) {
//				min = Double.MAX_VALUE;
//				for (String edt : edts.keySet()) {
//					if (edts.get(edt) < min && currentDataset.equals(edt.split("->")[0])
//							&& !previouslyVisitedDatacenters.contains(edt.split("->")[1])) {
//						min = edts.get(edt);
//						executionOrder = edt;
//					}
//				}
//				counter++;
//				currentDataset = executionOrder.split("->")[1];
//				previouslyVisitedDatacenters.add(executionOrder.split("->")[0]);
//				edtsInOrder.add(edts.get(executionOrder));
//				if (counter >= nrTimesDataCentersAreAccessed - 1) {
//					allDatacentersVisited = true;
//					previouslyVisitedDatacenters.add(currentDataset);
//				}
//			}
//
//			// MAYBE HERE END THE EDTS CALCULATION
//			
//			// Just printing stuff
////			System.out.println("Computation " + index);
////			for (String s : previouslyVisitedDatacenters) {
////				System.out.print(s + "->");
////			}
////			System.out.println();
////			for (Double d : edtsInOrder) {
////				System.out.print(d.doubleValue() + "->");
////			}
////			System.out.println();
//			//
//
//			// Calculate EPT for each data center
//			// get from uMatrixZ all datacenters where there is computation going on
//			List<Integer> usedDatacenters = new ArrayList<Integer>();
//			for (int i = 0; i < uMatrixZ[1].length; i++) {
//				if (uMatrixZ[index][i] == 1) {
//					usedDatacenters.add(i);
//				}
//			}
//			// get datacenter by id, get datasets from computation and see which match
//			List<Double> estimatedProcessingTimes = new ArrayList<Double>();
//			for (Integer idx : usedDatacenters) {
//				double taskSize = 0.0;
//				DataCenter dataCenter = getDataCenterAtIndex(idx);
//				List<Integer> computationDatasets = c.getProcessedDatasets();
//				for (Dataset d : dataCenter.getDatasetsToProcess()) {
//					if (computationDatasets.contains(d.getIndex())) {
//						taskSize += d.getFileSize();
//					}
//				}
//
//				estimatedProcessingTimes.add(taskSize / dataCenter.getEstimatedProcessingCapability());
//
//			}
//
		// calculate totalEPT, just add all edts and all epts
//			for (Double d : estimatedProcessingTimes) {
//				totalEPT += d;
//			}
//
//		}
		return totalEPT;
	}

	public static Dataset getDatasetAtIndex(int index) {
		for (Dataset d : datasets) {
			if (d.getIndex() == index) {
				return d;
			}
		}

		return null;
	}

	public static DataCenter getDataCenterAtIndex(int index) {
		for (DataCenter d : dataCenters) {
			if (d.getIndex() == index) {
				return d;
			}
		}
		return null;
	}

	public static double getInfoFromDatacenter(int index) {
		return getDataCenterAtIndex(index).getEstimatedBandwidth();
	}

	public static double getPrcCapFromDatacenter(int index) {
		return getDataCenterAtIndex(index).getEstimatedProcessingCapability();
	}

	public static int getLatencyBetweenDatacenter(int index1, int index2) {
		return Utils.latencyMatrix[index1][index2];
	}

	public static int[][] generateLatencyMatrix() {
		Random randNum = new Random();
		int[][] latencyMatrix = new int[Utils.dataCentersCount][Utils.dataCentersCount];
		for (int i = 0; i < Utils.dataCentersCount; i++) {
			for (int j = i; j < Utils.dataCentersCount; j++) {
				if (j != i) {
					latencyMatrix[i][j] = 40 + randNum.nextInt(40);
					latencyMatrix[j][i] = latencyMatrix[i][j];
				}
			}
		}
		// displayMatrix(latencyMatrix);
		return latencyMatrix;
	}

	public static int getMinFromMAtrix(double[][] matrix, String s, int idx) {
		double min = Double.MAX_VALUE;
		int index = 0;

		if (s.equals("row")) {
			for (int i = 0; i < matrix[1].length; i++) {
				if (matrix[idx][i] < min) {
					index = i;
					min = matrix[idx][i];
				}
			}
		} else if (s.equals("col")) {
			for (int i = 0; i < matrix.length; i++) {
				if (matrix[i][idx] < min) {
					index = i;
					min = matrix[i][idx];
				}
			}
		}

		return index;
	}
	
	public static List<Dataset> getDatasets(){
		return datasets;
	}

	public static void methodToTime() {
		BufferedReader br = null;
		try {
			// C:\\Users\\Andrei\\Desktop\\Files\\Test-Data.txt
			br = new BufferedReader(new FileReader("out1.txt"));

			String line = "";
			int j = 0;
			while ((line = br.readLine()) != null) {
				String[] values = line.split(" ");
				values[0] = values[0] + values[0];
				values[0] = "sdfsdfsdf";
				values[0] = "gfdgdf";
				values[0] = "sdfhgdfgsdfsdf";
				values[0] = "sdfsdfsdf";
				values[0] = "gfdgdf";
				values[0] = "sdfhgdfgsdfsdf";
				values[0] = "sdfsdfsdf";
				values[0] = "gfdgdf";
				values[0] = "sdfhgdfgsdfsdf";
				values[0] = "sdfsdfsdf";
				values[0] = "gfdgdf";
				values[0] = "sdfhgdfgsdfsdf";
				j += j + j + j - 10000;
				j++;
			}
			System.out.println(j);
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
